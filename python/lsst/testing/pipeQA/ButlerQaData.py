import sys, os, re, copy

import lsst.daf.persistence             as dafPersist
import lsst.afw.detection               as afwDet
import lsst.afw.image                   as afwImage
import lsst.meas.astrom                 as measAstrom
import lsst.afw.geom                    as afwGeom
import lsst.afw.cameraGeom              as cameraGeom
import lsst.meas.algorithms             as measAlg

# we need these for meas_extensions
# ... they are never explicitly used
try: import lsst.meas.extensions.shapeHSM
except: pass
try: import lsst.meas.extensions.rotAngle
except: pass
try: import lsst.meas.extensions.photometryKron
except: pass


import numpy
import math
import pyfits

import Manifest   as manifest
import CameraInfo as qaCamInfo

import QaDataUtils as qaDataUtils
import simRefObject as simRefObj
import source       as pqaSource

from QaData import QaData

#######################################################################
#
#
#
#######################################################################
class ButlerQaData(QaData):
    """ """

    #######################################################################
    #
    #######################################################################
    def __init__(self, label, rerun, cameraInfo, dataDir, **kwargs):
        """
        @param label A label to refer to the data
        @param rerun The rerun to retrieve
        @param cameraInfo A cameraInfo object describing the camera for these data
        @param dataDir The full path to the directory containing the data registry file.
        
        @param haveManifest verify files in dataDir are present according to manifest
        @param verifyChecksum verify files in dataDir have correct checksum as listed in manifest
        """
        
        QaData.__init__(self, label, rerun, cameraInfo)
        self.rerun = rerun
        self.dataDir = dataDir


        ###############################################
        # handle keyword args
        ###############################################
        self.kwargs      = kwargs
        self.dataId         = self.kwargs.get('dataId', {})
        self.haveManifest   = self.kwargs.get('haveManifest', False)
        self.verifyChecksum = self.kwargs.get('verifyChecksum', False)
        self.shapeAlg       = self.kwargs.get('shapeAlg', 'HSM_REGAUSS')

        knownAlgs = ["HSM_REGAUSS", "HSM_BJ", "HSM_LINEAR", "HSM_SHAPELET", "HSM_KSB"]
        if not self.shapeAlg in set(knownAlgs):
            knownStr = "\n".join(knownAlgs)
            raise Exception("Unknown shape algorithm: %s.  Please choose: \n%s\n" % (self.shapeAlg, knownStr))

        ###############################################
        # check the manifest, if requested
        # haveManifest = True is a bit slowish
        # verifyChecksum = True is quite slow
        manifest.verifyManifest(self.dataDir, verifyExists=self.haveManifest,
                                verifyChecksum=self.verifyChecksum)


        # This (dataId fetching) needs a better design, but will require butler/mapper change, I think.
        #
        # these obscure things refer to the names assigned to levels in the data hierarchy
        # eg. for lsstSim:   dataInfo  = [['visit',1], ['snap', 0], ['raft',0], ['sensor',0]]
        # a level is considered a discriminator if it represents different pictures of the same thing
        # ... so the same object may appear in multiple 'visits', but not on multiple 'sensors'
        # dataInfo is passed in from the derived class as it's specific to each mapper
        
        dataIdRegexDict = {}
        for array in self.dataInfo:
            dataIdName, dataIdDiscrim = array

            # if the user requested eg. visit=1234.*
            # pull that out of kwargs and put it in dataIdRegexDict
            if self.dataId.has_key(dataIdName):
                dataIdRegexDict[dataIdName] = self.dataId[dataIdName]
                

        #######################################
        # get butler
        self.outMapper = self.cameraInfo.getMapper(self.dataDir, rerun=self.rerun)
        self.outButler = dafPersist.ButlerFactory(mapper=self.outMapper).create()

        
        ####################################################
        # make a list of the frames we're asked to care about

        # get all the available raw inputs
        self.availableDataTuples = self.outButler.queryMetadata(cameraInfo.rawName, self.dataIdNames,
                                                                format=self.dataIdNames)

        # availableDataTuples may be a *very* *large* list.  Be sure to call reduceAvailableDataTupleList
        self.dataTuples = self.availableDataTuples

    
    def reduceAvailableDataTupleList(self, dataIdRegexDict):
        """Reduce availableDataTupleList by keeping only dataIds that match the input regex."""
        self.dataTuples = self._regexMatchDataIds(dataIdRegexDict, self.availableDataTuples)
        
    def initCache(self):

        QaData.initCache(self)
        # need to intialize these differently than base class
        # ... Db has 'object' and 'source' matching to be cached
        self.matchListCache = { 'obj': {}, 'src': {} }
        self.matchQueryCache = { 'obj' : {}, 'src': {} }

        
    def getDataName(self):
        """Get a string representation describing this data set. """
        return os.path.realpath(self.dataDir) + " rerun="+str(self.rerun)

    def getVisits(self, dataIdRegex):
        """ Return explicit visits matching for a dataIdRegex.

        @param dataIdRegex dataId dict containing regular expressions of data to retrieve.
        """
        visits = []

        dataTuplesToFetch = self._regexMatchDataIds(dataIdRegex, self.dataTuples, verbose=False)

        for dataTuple in dataTuplesToFetch:
            dataId = self._dataTupleToDataId(dataTuple)
            visits.append(self.cameraInfo.dataIdCameraToStandard(dataId)['visit'])
        return sorted(set(visits))
    


    def breakDataId(self, dataIdRegex, breakBy):
        """Take a dataId with regexes and return a list of dataId regexes
        which break the dataId by raft, or ccd.

        @param dataId    ... to be broken
        @param breakBy   'visit', 'raft', or 'ccd'
        """

        if not re.search("(visit|raft|ccd)", breakBy):
            raise Exception("breakBy must be 'visit','raft', or 'ccd'")

        if breakBy == 'visit':
            return [dataIdRegex]


        dataTuplesToFetch = self._regexMatchDataIds(dataIdRegex, self.dataTuples)


        dataIdDict = {}
        # handle lsst/hsc different naming conventions
        ccdConvention = 'ccd'
        if not dataIdRegex.has_key('ccd'):
            ccdConvention = 'sensor'

        for dataTuple in dataTuplesToFetch:
            thisDataId = self._dataTupleToDataId(dataTuple)
            visit, raft, ccd = thisDataId['visit'], thisDataId['raft'], thisDataId[ccdConvention]
            
            if breakBy == 'raft':
                ccd = dataIdRegex[ccdConvention]

            key = str(visit) + str(raft) + str(ccd)
            dataIdDict[key] = {
                'visit': str(visit),
                'raft' : raft,
                ccdConvention : ccd,
                'snap': '0'
                }
            
        # store the list of broken dataIds 
        self.brokenDataIdList = []
        for key in sorted(dataIdDict.keys()):
            self.brokenDataIdList.append(dataIdDict[key])
        
        return copy.copy(self.brokenDataIdList)
    

    def getMatchListBySensor(self, dataIdRegex, useRef=None):
        """Get a dict of all SourceMatches matching dataId, with sensor name as dict keys.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        flookup = {
            "u":"u", "g": "g", "r":"r", "i":"i", "z":"z", "y":"z",
            "B":"g", 'V':"r", 'R':"r", 'I':"i",
            }
        

        
        dataTuplesToFetch = self._regexMatchDataIds(dataIdRegex, self.dataTuples)
        
        # get the datasets corresponding to the request
        matchListDict = {}
        typeDict = {}
        for dataTuple in dataTuplesToFetch:
            dataId = self._dataTupleToDataId(dataTuple)
            dataKey = self._dataTupleToString(dataTuple)
            
            if self.matchListCache[useRef].has_key(dataKey):
                typeDict[dataKey] = copy.copy(self.matchListCache[useRef][dataKey])
                continue

            
            filterObj = self.getFilterBySensor(dataId)
            filterName = "unknown"
            if filterObj.has_key(dataKey) and hasattr(filterObj[dataKey], 'getName'):
                filterName = filterObj[dataKey].getName()
                filterName = flookup[filterName]
                
            # make sure we actually have the output file
            isWritten = self.outButler.datasetExists('icMatch', dataId) and \
                self.outButler.datasetExists('calexp', dataId)
            multiplicity = {}
            matchList = []
            
            if not isWritten:
                print str(dataTuple) + " output file missing.  Skipping."
                continue
            
            else:

                self.printStartLoad("Loading MatchList for: " + dataKey + "...")
                
                matches = measAstrom.astrom.readMatches(self.outButler, dataId)
                
                sourcesDict    = self.getSourceSetBySensor(dataId)
                refObjectsDict = self.getRefObjectSetBySensor(dataId)

                calibDict = self.getCalibBySensor(dataId)
                calib = calibDict[dataKey]

                
                fmag0, fmag0err = calib.getFluxMag0()
                for m in matches:
                    srefIn, sIn, dist = m
                    if ((not srefIn is None) and (not sIn is None)):

                        if not matchListDict.has_key(dataKey):
                            refCatObj = pqaSource.RefCatalog()
                            refCat    = refCatObj.catalog
                            catObj    = pqaSource.Catalog()
                            cat       = catObj.catalog

                            matchListDict[dataKey] = []

                            refRaKey   = refCatObj.keyDict['Ra']
                            refDecKey  = refCatObj.keyDict['Dec']
                            refPsfKey  = refCatObj.keyDict['PsfFlux']
                            refApKey   = refCatObj.keyDict['ApFlux']
                            refModKey  = refCatObj.keyDict['ModelFlux']
                            refInstKey = refCatObj.keyDict['InstFlux']

                            psfKey     = catObj.keyDict['PsfFlux']
                            apKey      = catObj.keyDict['ApFlux']
                            modKey     = catObj.keyDict['ModelFlux']
                            instKey    = catObj.keyDict['InstFlux']

                            psfErrKey  = catObj.keyDict['PsfFluxErr']
                            apErrKey   = catObj.keyDict['ApFluxErr']
                            modErrKey  = catObj.keyDict['ModelFluxErr']
                            instErrKey = catObj.keyDict['InstFluxErr']

                            intCenKey  = catObj.keyDict['FlagPixInterpCen']     
                            negKey     = catObj.keyDict['FlagNegative']    
                            edgeKey    = catObj.keyDict['FlagPixEdge']     
                            badCenKey  = catObj.keyDict['FlagBadCentroid'] 
                            satCenKey  = catObj.keyDict['FlagPixSaturCen']     
                            extKey     = catObj.keyDict['Extendedness']


                        matchList = matchListDict[dataKey]

                        # reference objects
                        sref = refCat.addNew()

                        sref.setId(srefIn.getId()) # this should be refobjId
                        sref.setD(refRaKey, srefIn.getRa().asDegrees())
                        sref.setD(refDecKey, srefIn.getDec().asDegrees())
                        flux = srefIn.get('flux')

                        sref.setD(refPsfKey, flux)
                        sref.setD(refApKey, flux)
                        sref.setD(refModKey, flux)
                        sref.setD(refInstKey, flux)

                        # sources
                        s = cat.addNew()
                        s.setId(sIn.getId())
                        isStar = srefIn.get('stargal')+0.0
                        s.setD(catObj.keyDict['Extendedness'], isStar)

                        s.setD(catObj.keyDict['XAstrom'], sIn.getX())
                        s.setD(catObj.keyDict['YAstrom'], sIn.getY())
                        s.setD(catObj.keyDict['Ra'], sIn.getRa().asDegrees())
                        s.setD(catObj.keyDict['Dec'], sIn.getDec().asDegrees())
                        s.setD(psfKey,  sIn.getPsfFlux())
                        s.setD(apKey,   sIn.getApFlux())
                        s.setD(modKey,  sIn.getModelFlux())
                        s.setD(instKey, sIn.getInstFlux())

                        s.setD(intCenKey, sIn.get('flags.pixel.interpolated.center')+0.0)
                        s.setD(negKey,    sIn.get('flags.negative')+0.0)
                        s.setD(edgeKey,   sIn.get('flags.pixel.edge')+0.0)
                        s.setD(badCenKey, sIn.get('flags.badcentroid')+0.0)
                        s.setD(satCenKey, sIn.get('flags.pixel.saturated.center')+0.0)
                        s.setD(extKey,    sIn.get('classification.extendedness')+0.0)

                        fmag0, fmag0Err = calib.getFluxMag0()

                        # fluxes
                        s.setD(psfKey,   s.getD(psfKey)/fmag0)
                        s.setD(apKey,    s.getD(apKey)/fmag0)
                        s.setD(modKey,   s.getD(modKey)/fmag0)
                        s.setD(instKey,  s.getD(instKey)/fmag0)

                        # flux errors
                        psfFluxErr  = qaDataUtils.calibFluxError(sIn.getPsfFlux(), sIn.getPsfFluxErr(),
                                                                 fmag0, fmag0Err)
                        s.setD(psfErrKey, psfFluxErr)

                        apFluxErr   = qaDataUtils.calibFluxError(sIn.getApFlux(),  sIn.getApFluxErr(),
                                                                 fmag0, fmag0Err)
                        s.setD(apErrKey, apFluxErr)

                        modFluxErr  = qaDataUtils.calibFluxError(sIn.getModelFlux(), sIn.getModelFluxErr(),
                                                                 fmag0, fmag0Err)
                        s.setD(modErrKey, modFluxErr)

                        instFluxErr = qaDataUtils.calibFluxError(sIn.getInstFlux(),  sIn.getInstFluxErr(),
                                                                 fmag0, fmag0Err)
                        s.setD(instErrKey, instFluxErr)

                        if multiplicity.has_key(s.getId()):
                            multiplicity[s.getId()] += 1
                        else:
                            multiplicity[s.getId()] = 1

                        matchList.append([sref, s, dist])

                self.dataIdLookup[dataKey] = dataId
                

                ######
                ######
                ######
                # Determine which are orphans, blends, straight matches, and non-detections

                sources    = sourcesDict[dataKey]
                if refObjectsDict.has_key(dataKey):
                    refObjects = refObjectsDict[dataKey]
                else:
                    refObjects = simRefObj.SimRefObjectSet() # an empty set


                typeDict[dataKey] = {}

                refIds     = []
                for ro in refObjects:
                    refIds.append(ro.getId())

                srcIds     = []
                for so in sources:
                    srcIds.append(so.getId())

                matRef     = []
                matSrc     = []
                for ma in matchList:
                    matRef.append(ma[0].getId())
                    matSrc.append(ma[1].getId())

                refIds = set(refIds)
                srcIds = set(srcIds)
                matRef = set(matRef)
                matSrc = set(matSrc)

                undetectedIds = refIds - matRef
                orphanIds     = srcIds - matSrc
                matchedIds    = srcIds & matSrc   # Does not know about duplicates
                #print 'Undet, orphan, matched:', len(undetectedIds), len(orphanIds), len(matchedIds)

                undetected = []
                orphans    = []
                matched    = []
                blended    = []
                for ro in refObjects:
                    if ro.getId() in undetectedIds:
                        undetected.append(ro)
                matchListById = dict([(m[1].getId(), m) for m in matchList])
                matchIdsSet = set(matchListById.keys())
                for so in sources:
                    soid = so.getId()
                    if soid in orphanIds:
                        orphans.append(so)
                    if soid in matchedIds and soid in matchIdsSet:
                        if multiplicity[soid] == 1:
                            matched.append(matchListById[soid])
                        else:
                            blended.append(matchListById[soid])



                typeDict[dataKey]['orphan']     = orphans
                typeDict[dataKey]['matched']    = matchList # a hack b/c src and icSrc Ids are different
                typeDict[dataKey]['blended']    = blended
                typeDict[dataKey]['undetected'] = undetected

                # cache it
                self.matchListCache[useRef][dataKey] = typeDict[dataKey]


                self.printMidLoad('\n        %s: Undet, orphan, matched, blended = %d %d %d %d' % (
                    dataKey,
                    len(typeDict[dataKey]['undetected']),
                    len(typeDict[dataKey]['orphan']),
                    len(typeDict[dataKey]['matched']),
                    len(typeDict[dataKey]['blended']))
                                  )

                self.matchQueryCache[useRef][dataKey] = True

                # Determine which are orphans, blends, straight matches, and non-detections
                ######
                ######
                ######

                self.printStopLoad()

        return copy.copy(typeDict)



    #######################################################################
    #
    #######################################################################
    def getSourceSetBySensor(self, dataIdRegex):
        """Get a dict of all Sources matching dataId, with sensor name as dict keys.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """
        
        dataTuplesToFetch = self._regexMatchDataIds(dataIdRegex, self.dataTuples)

        # get the datasets corresponding to the request
        ssDict = {}
        for dataTuple in dataTuplesToFetch:
            dataId = self._dataTupleToDataId(dataTuple)
            dataKey = self._dataTupleToString(dataTuple)
            
            if self.sourceSetCache.has_key(dataKey):
                ssDict[dataKey] = copy.copy(self.sourceSetCache[dataKey])
                continue

            self.printStartLoad("Loading SourceSets for: " + dataKey + "...")
            
            # make sure we actually have the output file
            isWritten = self.outButler.datasetExists('src', dataId)
            if isWritten:
                sourceCatalog = self.outButler.get('src', dataId)

                calibDict = self.getCalibBySensor(dataId)
                calib = calibDict[dataKey]

                if not calib is None:
                    fmag0, fmag0Err = calib.getFluxMag0()
                else:
                    print "Warning: no calib available, fluxes uncalibrated."
                    fmag0, fmag0Err = 1.0, 1.0


                catObj = pqaSource.Catalog()
                cat  = catObj.catalog

                raKey  = catObj.keyDict['Ra']     
                decKey = catObj.keyDict['Dec']    
                xKey   = catObj.keyDict['XAstrom']
                yKey   = catObj.keyDict['YAstrom']

                ixxKey = catObj.keyDict['Ixx']
                iyyKey = catObj.keyDict['Iyy']
                ixyKey = catObj.keyDict['Ixy']

                intCenKey  = catObj.keyDict['FlagPixInterpCen']     
                negKey     = catObj.keyDict['FlagNegative']    
                edgeKey    = catObj.keyDict['FlagPixEdge']     
                badCenKey  = catObj.keyDict['FlagBadCentroid'] 
                satCenKey  = catObj.keyDict['FlagPixSaturCen']     
                extKey     = catObj.keyDict['Extendedness']


                psfKey = catObj.keyDict['PsfFlux']
                apKey  = catObj.keyDict['ApFlux']
                modKey = catObj.keyDict['ModelFlux']
                instKey = catObj.keyDict['InstFlux']
            
                psfErrKey = catObj.keyDict['PsfFluxErr']
                apErrKey  = catObj.keyDict['ApFluxErr']
                modErrKey = catObj.keyDict['ModelFluxErr']
                instErrKey = catObj.keyDict['InstFluxErr']
                    
                for s in sourceCatalog:

                    rec = cat.addNew()
                    rec.setId(s.getId())

                    rec.setD(raKey,    float(s.getRa()))
                    rec.setD(decKey,   float(s.getDec()))
                    rec.setD(xKey,     float(s.getX())) #Astrom()))
                    rec.setD(yKey,     float(s.getY())) #Astrom()))
                    
                    # fluxes
                    rec.setD(psfKey,   float(s.getPsfFlux())/fmag0)
                    rec.setD(apKey,    float(s.getApFlux())/fmag0)
                    rec.setD(modKey,   float(s.getModelFlux())/fmag0)
                    rec.setD(instKey,  float(s.getInstFlux())/fmag0)

                    # shapes
                    rec.setD(ixxKey,   float(s.getIxx()))
                    rec.setD(iyyKey,   float(s.getIyy()))
                    rec.setD(ixyKey,   float(s.getIxy()))
                    
                    # flags
                    rec.setD(intCenKey, s.get('flags.pixel.interpolated.center')+0.0)
                    rec.setD(negKey,    s.get('flags.negative')+0.0)
                    rec.setD(edgeKey,   s.get('flags.pixel.edge')+0.0)
                    rec.setD(badCenKey, s.get('flags.badcentroid')+0.0)
                    rec.setD(satCenKey, s.get('flags.pixel.saturated.center')+0.0)
                    rec.setD(extKey,    s.get('classification.extendedness')+0.0)
                    
                    
                    # flux errors
                    psfFluxErr  = qaDataUtils.calibFluxError(float(s.getPsfFlux()), float(s.getPsfFluxErr()),
                                                             fmag0, fmag0Err)
                    rec.setD(psfErrKey, psfFluxErr)

                    apFluxErr   = qaDataUtils.calibFluxError(float(s.getApFlux()),  float(s.getApFluxErr()),
                                                             fmag0, fmag0Err)
                    rec.setD(apErrKey, apFluxErr)
                    
                    modFluxErr  = qaDataUtils.calibFluxError(float(s.getModelFlux()), float(s.getModelFluxErr()),
                                                             fmag0, fmag0Err)
                    rec.setD(modErrKey, modFluxErr)
                    
                    instFluxErr = qaDataUtils.calibFluxError(float(s.getInstFlux()),  float(s.getInstFluxErr()),
                                                             fmag0, fmag0Err)
                    rec.setD(instErrKey, instFluxErr)


                self.sourceSetCache[dataKey] = catObj.catalog
                ssDict[dataKey] = copy.copy(catObj.catalog)
                self.dataIdLookup[dataKey] = dataId

 
            else:
                print str(dataTuple) + " output file missing.  Skipping."
                
            self.printStopLoad()
                
        return ssDict

    def getSourceSet(self, dataIdRegex):
        """Get a SourceSet of all Sources matching dataId.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        ssDict = self.getSourceSetBySensor(dataIdRegex)
        ssReturn = []
        for key, ss in ssDict.items():
            ssReturn += ss
        return ssReturn




    def getRefObjectSetBySensor(self, dataIdRegex):
        """Get a dict of all Catalog Sources matching dataId, with sensor name as dict keys.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """
        
        # if the dataIdRegex is identical to an earlier query, we must already have all the data
        dataIdStr = self._dataIdToString(dataIdRegex)
        if self.refObjectQueryCache.has_key(dataIdStr):
            sroDict = {}
            # get only the ones that match the request
            for key, sro in self.refObjectCache.items():
                if re.search(dataIdStr, key):
                    sroDict[key] = sro
            return sroDict

        self.printStartLoad("Loading RefObjects for: " + dataIdStr + "...")
        
        # parse results and put them in a sourceSet
        mastConfig = measAstrom.astrom.MeasAstromConfig()
        astrom = measAstrom.astrom.Astrometry(mastConfig)
        
        dataTuplesToFetch = self._regexMatchDataIds(dataIdRegex, self.dataTuples)

        # get the datasets corresponding to the request
        sroDict = {}
        for dataTuple in dataTuplesToFetch:
            dataId = self._dataTupleToDataId(dataTuple)
            dataKey = self._dataTupleToString(dataTuple)
            
            wcs = self.getWcsBySensor(dataId)[dataKey]
            filterName = self.getFilterBySensor(dataId)[dataKey].getName()
            imageSize = self.calexpCache[dataKey]['NAXIS1'], self.calexpCache[dataKey]['NAXIS2']
            pixelMargin = 0.0
            refCat = astrom.getReferenceSourcesForWcs(wcs, imageSize, filterName, pixelMargin)
        
            if not sroDict.has_key(dataKey):
                sroDict[dataKey] = []
            sros = sroDict[dataKey]
            
            for rec in refCat:
                sro = simRefObj.SimRefObject()
                sro.refObjectId = rec.getId()
                sro.isStar = rec.get('stargal') + 0

                coo = rec.get('coord')
                sro.setRa(coo.getRa().asDegrees())
                sro.setDecl(coo.getDec().asDegrees())
                sro.setFlux(rec.get('flux'), filterName)

                sros.append(sro)

        self.printStopLoad()
        
        # cache it
        self.refObjectQueryCache[dataIdStr] = True
        for k, sro in sroDict.items():
            self.refObjectCache[k] = sroDict[k]

        return sroDict





    def loadCalexp(self, dataIdRegex):
        """Load the calexp data for data matching dataIdRegex.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """
        
        dataTuplesToFetch = self._regexMatchDataIds(dataIdRegex, self.dataTuples)

        # get the datasets corresponding to the request
        for dataTuple in dataTuplesToFetch:
            dataId = self._dataTupleToDataId(dataTuple)
            dataKey = self._dataTupleToString(dataTuple)

            if self.calexpCache.has_key(dataKey):
                continue

            self.printStartLoad("Loading Calexp for: " + dataKey + "...")

            if self.outButler.datasetExists('calexp_md', dataId):
                calexp_md = self.outButler.get('calexp_md', dataId)
                
                self.wcsCache[dataKey]      = afwImage.makeWcs(calexp_md)

                ccdName = calexp_md.getAsString('DETNAME').strip()
                names = ccdName.split()

                if len(names) > 1:
                    raftName = names[0]
                else:
                    raftName = ""
                raftName = raftName.strip()

                self.detectorCache[dataKey] = self.cameraInfo.detectors[ccdName] #ccdDetector
                if len(raftName) > 0:
                    self.raftDetectorCache[dataKey] = self.cameraInfo.detectors[raftName]

                
                self.filterCache[dataKey]   = afwImage.Filter(calexp_md)
                self.calibCache[dataKey]    = afwImage.Calib(calexp_md)

                # store the calexp as a dict
                if not self.calexpCache.has_key(dataKey):
                    self.calexpCache[dataKey] = {}

                nameLookup = qaDataUtils.getCalexpNameLookup()
                for n in calexp_md.names():
                    val = calexp_md.get(n)
                    self.calexpCache[dataKey][n] = val

                    # assign an alias to provide the same name as the database version uses.
                    if nameLookup.has_key(n):
                        n2 = nameLookup[n]
                        self.calexpCache[dataKey][n2] = val

                # if we're missing anything in nameLookup ... put in a NaN
                for calexpName,qaName in nameLookup.items():
                    if not self.calexpCache[dataKey].has_key(qaName):
                        self.calexpCache[dataKey][qaName] = numpy.NaN


                # check the fwhm ... we know we need it
                # NOTE that we actually try to load this from the SEEING
                #  keyword in the calexp_md.  So a fwhm=-1 here, doesn't mean
                #  it wasn't already set by SEEING
                sigmaToFwhm = 2.0*math.sqrt(2.0*math.log(2.0))
                width = calexp_md.get('NAXIS1')
                height = calexp_md.get('NAXIS2')
                try:
                    psf = self.outButler.get("psf", visit=dataId['visit'],
                                             raft=dataId['raft'], sensor=dataId['sensor'])
                    attr = measAlg.PsfAttributes(psf, width // 2, height // 2)
                    fwhm = attr.computeGaussianWidth() * self.wcsCache[dataKey].pixelScale().asArcseconds() * sigmaToFwhm
                except Exception, e:
                    fwhm = -1.0

                if (self.calexpCache[dataKey].has_key('fwhm') and
                    numpy.isnan(self.calexpCache[dataKey]['fwhm'])):
                    self.calexpCache[dataKey]['fwhm'] = fwhm
                
                self.calexpQueryCache[dataKey] = True
                self.dataIdLookup[dataKey] = dataId


                
            else:
                calibFilename = self.outButler.get('calexp_filename', dataId)
                print "\nSkipping " + str(dataTuple) + ". Calib output file missing:"
                print "   "+str(calibFilename)
                

            self.printStopLoad()
            


            
    def getCalexpEntryBySensor(self, cache, dataIdRegex):
        """Fill and return the dict for a specified calexp cache.

        @param cache The cache dictionary to return
        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        dataTuplesToFetch = self._regexMatchDataIds(dataIdRegex, self.dataTuples)

        # get the datasets corresponding to the request
        entryDict = {}
        for dataTuple in dataTuplesToFetch:
            dataId = self._dataTupleToDataId(dataTuple)
            dataKey = self._dataTupleToString(dataTuple)

            self.loadCalexp(dataId)
            if cache.has_key(dataKey):
                entryDict[dataKey] = cache[dataKey]
            else:
                entryDict[dataKey] = None
                
        return entryDict



    #######################################################################
    # utility to go through a list of data Tuples and return
    #  the ones which match regexes for the corresponding data type
    # so user can say eg. raft='0,\d', visit='855.*', etc
    #######################################################################
    def _regexMatchDataIds(self, dataIdRegexDict, availableDataTuples, verbose=False):
        """Match available data with regexes in a dataId dictionary        
        
        @param dataIdRegexDict dataId dict of regular expressions for data to be handled.
        @param availableDataTuples data sets available to be retrieved.
        """

        # go through the list of what's available, and compare to what we're asked for
        # Put matches in a list of tuples, eg. [(vis1,sna1,raf1,sen1),(vis2,sna2,raf2,sen2)] 
        dataTuples = []
        for dataTuple in availableDataTuples:
            if verbose:
                print dataTuple
            # start true and fail if any dataId keys fail ... eg. 'visit' doesn't match
            match = True
            for i in range(len(self.dataIdNames)):
                dataIdName = self.dataIdNames[i]   # eg. 'visit', 'sensor', etc
                regexForThisId = dataIdRegexDict.get(dataIdName, '.*') # default to '.*' or 'anything'
                dataId = dataTuple[i]

                # if it doesn't match, this frame isn't to be run.
                if not re.search(str(regexForThisId),  str(dataId)):
                    match = False
                    break
                
                # ignore the guiding ccds on the hsc camera
                if re.search('^hsc.*', self.cameraInfo.name) and dataIdName == 'ccd' and dataId > 99:
                    match = False
                    break

            if match:
                dataTuples.append(dataTuple)
                
        return dataTuples
                

    



#######################################################################
#
#
#
#######################################################################
def makeButlerQaData(label, rerun=None, camera=None, **kwargs):
    """Factory for a ButlerQaData object.

    @param label The basename directory in a TESTBED_PATH directory where the registry file exists.
    @param rerun The rerun of the data to retrieve

    Keyword args passed to ButlerQaData constructor:
    @param haveManifest verify files in dataDir are present according to manifest
    @param verifyChecksum verify files in dataDir have correct checksum as listed in manifest
    """
        
    testbedDir, testdataDir = qaDataUtils.findDataInTestbed(label)

    # make sure LsstSim is last in the list (its 'verifyRegistries()' will pass for all cameras)
    cameraKeys = ["hsc", "suprimecam", "suprimecam-old", "sdss", "lsstsim"]
    cameraInfos = {
#       "cfht": qaCamInfo.CfhtCameraInfo(), # XXX CFHT camera geometry is currently broken following #1767
        "hsc" : qaCamInfo.HscCameraInfo(),
        "suprimecam": qaCamInfo.SuprimecamCameraInfo(),
        "suprimecam-old": qaCamInfo.SuprimecamCameraInfo(True),
        "sdss" : qaCamInfo.SdssCameraInfo(),
        "lsstsim": qaCamInfo.LsstSimCameraInfo(),
        }

    # try each camera ***in-order***
    # note that LsstSim will look like all other as it uses the same registry for data and calib
    # ... must test it last
    cameraToUse = None
    if not camera is None:
        cameraToUse = cameraInfos[camera]
    else:

        for cameraKey in cameraKeys:
            cameraInfo = cameraInfos[cameraKey]
            # if the mapper couldn't be found, we can't use this camera
            hasMapper = not cameraInfo.mapperClass is None
            validReg = cameraInfo.verifyRegistries(testdataDir)
            #print cameraInfo.name, "valid: ", validReg
            if hasMapper and validReg:
                cameraToUse = cameraInfo
                break

    if cameraToUse is None:
        raise Exception("Can't find registries usable with any mappers.")
    else:
        if rerun is None:
            rerun = cameraToUse.getDefaultRerun()
        print "label:       ", label
        print "rerun:       ", rerun
        print "camera:      ", cameraToUse.name
        print "testdataDir: ", testdataDir
        return ButlerQaData(label, rerun, cameraToUse, testdataDir, **kwargs)


