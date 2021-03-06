#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import sys, os, re, copy
import datetime

import lsst.daf.persistence             as dafPersist
import lsst.afw.detection               as afwDet
import lsst.afw.image                   as afwImage
import lsst.meas.astrom                 as measAstrom
import lsst.afw.geom                    as afwGeom
import lsst.afw.cameraGeom              as cameraGeom
import lsst.meas.algorithms             as measAlg

import lsst.meas.photocal               as measPhotCal
import lsst.obs.hsc.colorterms       as osubColorTerms
# we need these for meas_extensions
# ... they are never explicitly used
try: import lsst.meas.extensions.shapeHSM
except: pass
try: import lsst.meas.extensions.rotAngle
except: pass
try: import lsst.meas.extensions.photometryKron
except: pass

haveMosaic=True
try: import lsst.meas.mosaic as measMos
except: haveMosaic=False

import numpy
import math
import pyfits

import CameraInfo as qaCamInfo

from QaDataUtils import QaDataUtils
qaDataUtils = QaDataUtils()

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
        """
        
        QaData.__init__(self, label, rerun, cameraInfo, qaDataUtils, **kwargs)
        self.rerun = rerun
        self.dataDir = dataDir


        ###############################################
        # handle keyword args
        ###############################################
        self.kwargs      = kwargs
        self.dataId         = self.kwargs.get('dataId', {})
        self.shapeAlg       = self.kwargs.get('shapeAlg', 'HSM_REGAUSS')

        knownAlgs = ["HSM_REGAUSS", "HSM_BJ", "HSM_LINEAR", "HSM_SHAPELET", "HSM_KSB"]
        if not self.shapeAlg in set(knownAlgs):
            knownStr = "\n".join(knownAlgs)
            raise Exception("Unknown shape algorithm: %s.  Please choose: \n%s\n" % (self.shapeAlg, knownStr))


        # these obscure things refer to the names assigned to levels in the data hierarchy
        # eg. for lsstSim:   dataInfo  = [['visit',1], ['snap', 0], ['raft',0], ['sensor',0]]
        # a level is considered a discriminator if it represents different pictures of the same thing
        # ... so the same object may appear in multiple 'visits', but not on multiple 'sensors'
        # dataInfo is passed in from the derived class as it's specific to each camera
        
        dataIdRegexDict = {}
        for array in self.dataInfo:
            dataIdName, dataIdDiscrim = array

            # if the user requested eg. visit=1234.*
            # pull that out of kwargs and put it in dataIdRegexDict
            if self.dataId.has_key(dataIdName):
                dataIdRegexDict[dataIdName] = self.dataId[dataIdName]
                

        #######################################
        # get butler
        self.butler = dafPersist.Butler(self.dataDir)

        
        ####################################################
        # make a list of the frames we're asked to care about

        # get all the available raw inputs
        self.availableDataTuples = self.butler.queryMetadata(cameraInfo.rawName, self.dataIdNames,
                                                                format=self.dataIdNames)

        # availableDataTuples may be a *very* *large* list.  Be sure to call reduceAvailableDataTupleList
        self.dataTuples = self.availableDataTuples

        self.alreadyTriedCalexp = set()
        
        
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

        dataTuplesToFetch = self._regexMatchDataIds(dataIdRegex, self.dataTuples, exact=False, verbose=False)

        for dataTuple in dataTuplesToFetch:
            dataId = self._dataTupleToDataId(dataTuple)
            visits.append(self.cameraInfo.dataIdCameraToStandard(dataId)['visit'])
        visits_sort = sorted(set(visits))
        return visits_sort


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

        dataIdDict = {}
        # handle lsst/hsc different naming conventions
        ccdConvention = 'ccd'
        if not dataIdRegex.has_key('ccd'):
            ccdConvention = 'sensor'

        
        exact = False
        if re.search("^\d+$", dataIdRegex[ccdConvention]):
            exact = True
        dataTuplesToFetch = self._regexMatchDataIds(dataIdRegex, self.dataTuples, exact=exact)


        for dataTuple in dataTuplesToFetch:
            thisDataId = self._dataTupleToDataId(dataTuple)
            visit = thisDataId['visit']
            raft = "NoRaft"
            if 'raft' in thisDataId:
                raft = thisDataId['raft']
            ccd = thisDataId[ccdConvention]
            
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
        
        dataIdListCopy = copy.copy(self.brokenDataIdList)
        return dataIdListCopy

    
    def verify(self, dataId):
        # just load the calexp, you'll need it anyway
        self.loadCalexp(dataId)
        key = self._dataIdToString(dataId, defineFully=True)
        haveIt = True if key in self.calexpQueryCache else False
        return haveIt
    
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
            isWritten = self.butler.datasetExists('icMatch', dataId) and \
                self.butler.datasetExists('calexp', dataId)
            multiplicity = {}
            matchList = []
            
            if not isWritten:
                self.log.log(self.log.WARN, str(dataTuple) + " output file missing.  Skipping.")
                continue
            
            else:

                self.printStartLoad("Loading MatchList for: " + dataKey + "...")
                
                matches = measAstrom.astrom.readMatches(self.butler, dataId)

                # use the ref fluxes to get colors and color correct
                astrom = measAstrom.astrom.Astrometry(measAstrom.config.MeasAstromConfig())
                _rad = 0.1*afwGeom.arcseconds
                sourcesDict    = self.getSourceSetBySensor(dataId)
                refObjectsDict = self.getRefObjectSetBySensor(dataId)

                calibDict = self.getCalibBySensor(dataId)
                calib = calibDict[dataKey]

                measPhotCal.colorterms.Colorterm.setColorterms(osubColorTerms.colortermsData)
                measPhotCal.colorterms.Colorterm.setActiveDevice("Hamamatsu") # ahghgh ... hard code
        
                cterm = osubColorTerms.colortermsData['Hamamatsu'][filterName]
                
                fmag0, fmag0err = calib.getFluxMag0()
                fmag0err = 0.0
                for m in matches:
                    srefIn, sIn, dist = m
                    if ((srefIn is not None) and (sIn is not None)):

                        if not matchListDict.has_key(dataKey):
                            refCatObj = pqaSource.RefCatalog()
                            refCat    = refCatObj.catalog
                            catObj    = pqaSource.Catalog()
                            cat       = catObj.catalog

                            matchListDict[dataKey] = []


                        matchList = matchListDict[dataKey]

                        # reference objects
                        sref = refCat.addNew()

                        sref.setId(srefIn.getId()) # this should be refobjId

                        fmag0, fmag0Err = calib.getFluxMag0()
                        
                        _ra = srefIn.getRa()
                        _dec = srefIn.getDec()
                        fullRefCat = astrom.getReferenceSources(_ra, _dec, _rad, filterName, allFluxes=True)
                        mPrimary   = -2.5*numpy.log10(fullRefCat.get(cterm.primary))
                        mSecondary = -2.5*numpy.log10(fullRefCat.get(cterm.secondary))
                        refMag = cterm.transformMags(filterName, mPrimary, mSecondary)
                        calmag = -2.5*numpy.log10(sIn.getApFlux()/fmag0)
                        #print mPrimary, refMag, refMag2, calmag, mPrimary-mSecondary, mPrimary-calmag, refMag-calmag, refMag2-calmag


                        
                        sref.setD(self.k_rRa, _ra.asDegrees())
                        sref.setD(self.k_rDec, _dec.asDegrees())
                        #flux = srefIn.get('flux')
                        flux = 10**(-refMag[0]/2.5)
                        refMerr = fullRefCat.get(cterm.primary+".err")
                        ferr = (flux*numpy.log(10.0)*0.4*refMerr)[0]
                        
                        sref.setD(self.k_rPsf, flux)
                        sref.setD(self.k_rAp, flux)
                        sref.setD(self.k_rMod, flux)
                        sref.setD(self.k_rInst, flux)
                        sref.setD(self.k_rPsfE, ferr)
                        sref.setD(self.k_rApE, ferr)
                        sref.setD(self.k_rModE, ferr)
                        sref.setD(self.k_rInstE, ferr)

                        # sources
                        s = cat.addNew()
                        s.setId(sIn.getId())
                        isStar = 0
                        if 'stargal' in srefIn.getSchema().getNames():
                            isStar = srefIn.get('stargal')
                        s.setD(self.k_ext, isStar)

                        s.setD(self.k_x,    sIn.getX())
                        s.setD(self.k_y,    sIn.getY())
                        s.setD(self.k_Ra,   sIn.getRa().asDegrees())
                        s.setD(self.k_Dec,  sIn.getDec().asDegrees())
                        s.setD(self.k_Psf,  sIn.getPsfFlux())
                        s.setD(self.k_Ap,   sIn.getApFlux())
                        s.setD(self.k_Mod,  sIn.getModelFlux())
                        s.setD(self.k_Inst, sIn.getInstFlux())
                        s.setI(self.k_intc, sIn.get('flags.pixel.interpolated.center'))
                        s.setI(self.k_neg,  sIn.get('flags.negative'))
                        s.setI(self.k_edg,  sIn.get('flags.pixel.edge'))
                        s.setI(self.k_bad,  sIn.get('flags.badcentroid'))
                        s.setI(self.k_satc, sIn.get('flags.pixel.saturated.center'))
                        s.setD(self.k_ext,  sIn.get('classification.extendedness'))

                        # fluxes
                        s.setD(self.k_Psf,   s.getD(self.k_Psf)/fmag0)
                        s.setD(self.k_Ap,    s.getD(self.k_Ap)/fmag0)
                        s.setD(self.k_Mod,   s.getD(self.k_Mod)/fmag0)
                        s.setD(self.k_Inst,  s.getD(self.k_Inst)/fmag0)

                        # flux errors
                        psfFluxErr  = qaDataUtils.calibFluxError(sIn.getPsfFlux(), sIn.getPsfFluxErr(),
                                                                 fmag0, fmag0Err)
                        s.setD(self.k_PsfE, psfFluxErr)

                        apFluxErr   = qaDataUtils.calibFluxError(sIn.getApFlux(),  sIn.getApFluxErr(),
                                                                 fmag0, fmag0Err)
                        s.setD(self.k_ApE, apFluxErr)

                        modFluxErr  = qaDataUtils.calibFluxError(sIn.getModelFlux(), sIn.getModelFluxErr(),
                                                                 fmag0, fmag0Err)
                        s.setD(self.k_ModE, modFluxErr)

                        instFluxErr = qaDataUtils.calibFluxError(sIn.getInstFlux(),  sIn.getInstFluxErr(),
                                                                 fmag0, fmag0Err)
                        s.setD(self.k_InstE, instFluxErr)

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


                if len(matchList) == 0:
                    self.log.log(self.log.WARN, '%s: NO MATCHED OBJECTS!  Undet, orphan, matched, blended = %d %d %d %d' % (
                            dataKey, len(undetected), len(orphans), len(matchList), len(blended))
                                 )
                else:
                    self.log.log(self.log.INFO, '%s: Undet, orphan, matched, blended = %d %d %d %d' % (
                            dataKey, len(undetected), len(orphans), len(matchList), len(blended))
                                 )


                self.matchQueryCache[useRef][dataKey] = True

                self.printStopLoad("MatchList load for: " + dataKey)

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
            isWritten = self.butler.datasetExists('src', dataId)
            if isWritten:
                sourceCatalog = self.butler.get('src', dataId)

                calibDict = self.getCalibBySensor(dataId)
                calib = calibDict[dataKey]

                if calib is not None:
                    fmag0, fmag0Err = calib.getFluxMag0()
                else:
                    self.log.log(self.log.WARN, "Warning: no calib available, fluxes uncalibrated.")
                    fmag0, fmag0Err = 1.0, 1.0

                #fmag0Err = 0.0
                #print fmag0, fmag0Err
                catObj = pqaSource.Catalog()
                cat  = catObj.catalog


                    
                for s in sourceCatalog:

                    rec = cat.addNew()
                    rec.setId(s.getId())

                    rec.setD(self.k_Ra,    float(s.getRa().asDegrees()))
                    rec.setD(self.k_Dec,   float(s.getDec().asDegrees()))
                    rec.setD(self.k_x,     float(s.getX()))
                    rec.setD(self.k_y,     float(s.getY()))
                    
                    # fluxes
                    rec.setD(self.k_Psf,   float(s.getPsfFlux())/fmag0)
                    rec.setD(self.k_Ap,    float(s.getApFlux())/fmag0)
                    rec.setD(self.k_Mod,   float(s.getModelFlux())/fmag0)
                    rec.setD(self.k_Inst,  float(s.getInstFlux())/fmag0)

                    # shapes
                    rec.setD(self.k_ixx,   float(s.getIxx()))
                    rec.setD(self.k_iyy,   float(s.getIyy()))
                    rec.setD(self.k_ixy,   float(s.getIxy()))
                    
                    # flags
                    rec.setI(self.k_intc, s.get('flags.pixel.interpolated.center'))
                    rec.setI(self.k_neg,  s.get('flags.negative'))
                    rec.setI(self.k_edg,  s.get('flags.pixel.edge'))
                    rec.setI(self.k_bad,  s.get('flags.badcentroid'))
                    rec.setI(self.k_satc, s.get('flags.pixel.saturated.center'))
                    rec.setD(self.k_ext,  s.get('classification.extendedness'))
                    rec.setI(self.k_nchild, 0) #s.get('deblend_nchild'))
                    
                    # flux errors
                    psfFluxErr  = qaDataUtils.calibFluxError(float(s.getPsfFlux()), float(s.getPsfFluxErr()),
                                                             fmag0, fmag0Err)
                    rec.setD(self.k_PsfE, psfFluxErr)

                    apFluxErr   = qaDataUtils.calibFluxError(float(s.getApFlux()),  float(s.getApFluxErr()),
                                                             fmag0, fmag0Err)
                    rec.setD(self.k_ApE, apFluxErr)
                    
                    modFluxErr  = qaDataUtils.calibFluxError(float(s.getModelFlux()), float(s.getModelFluxErr()),
                                                             fmag0, fmag0Err)
                    rec.setD(self.k_ModE, modFluxErr)
                    
                    instFluxErr = qaDataUtils.calibFluxError(float(s.getInstFlux()),  float(s.getInstFluxErr()),
                                                             fmag0, fmag0Err)
                    rec.setD(self.k_InstE, instFluxErr)


                self.sourceSetCache[dataKey] = catObj.catalog
                ssDict[dataKey] = copy.copy(catObj.catalog)
                self.dataIdLookup[dataKey] = dataId

 
            else:
                self.log.log(self.log.WARN, str(dataTuple) + " output file missing.  Skipping.")
                
            self.printStopLoad("SourceSets load for: " + dataKey)
                
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
                sro.isStar = 0
                if 'stargal' in rec.getSchema().getNames():
                    sro.isStar = rec.get('stargal') + 0

                coo = rec.get('coord')
                sro.setRa(coo.getRa().asDegrees())
                sro.setDecl(coo.getDec().asDegrees())
                sro.setFlux(rec.get('flux'), filterName)

                sros.append(sro)

        self.printStopLoad("RefObjects load for: " + dataIdStr)
        
        # cache it
        self.refObjectQueryCache[dataIdStr] = True
        for k, sro in sroDict.items():
            self.refObjectCache[k] = sroDict[k]

        return sroDict



    def getSummaryDataBySensor(self, dataIdRegex):
        """Get a dict of dict objects which contain specific summary data.
        
        @param dataIdRegex dataId dictionary with regular expressions to specify data to retrieve
        """
        calexp = self.getCalexpBySensor(dataIdRegex)
        summary = {}
        for dataId, ce in calexp.items():
            if not dataId in summary:
                summary[dataId] = {}

            if ce is None:
                ce = {}

            dobs = ce.get('DATE-OBS', '0000-00-00')
            summary[dataId]["DATE_OBS"] = datetime.datetime.strptime(dobs, "%Y-%m-%d")
            summary[dataId]["EXPTIME"]  = ce.get('EXPTIME', 0.0)
            summary[dataId]['RA']       = ce.get('RA', 0.0)
            summary[dataId]['DEC']      = ce.get('DEC', 0.0)

            summary[dataId]['ALT']          = ce.get('ALTITUDE', 0.0)
            summary[dataId]['AZ']           = ce.get('AZIMUTH', 0.0)
            summary[dataId]["SKYLEVEL"]     = ce.get('SKYLEVEL', 0.0)
            summary[dataId]["ELLIPT"]       = ce.get('ELLIPT', 0.0)
            summary[dataId]["ELL_PA"]       = ce.get('ELL_PA_MED', 0.0)
            summary[dataId]["AIRMASS"]      = ce.get('AIRMASS', 0.0)
            summary[dataId]["FLATNESS_RMS"] = ce.get('FLATNESS_RMS', 0.0)
            summary[dataId]["FLATNESS_PP"]  = ce.get('FLATNESS_PP', 0.0)
            summary[dataId]["SIGMA_SKY"]    = ce.get('SKYSIGMA', 0.0)
            summary[dataId]["SEEING"]       = ce.get('SEEING', 0.0)
            summary[dataId]['OBJECT']       = ce.get('OBJECT', 0.0)
            summary[dataId]["OSLEVEL1"]     = ce.get('OSLEVEL1', 0.0)
            summary[dataId]["OSLEVEL2"]     = ce.get('OSLEVEL2', 0.0)
            summary[dataId]["OSLEVEL3"]     = ce.get('OSLEVEL3', 0.0)
            summary[dataId]["OSLEVEL4"]     = ce.get('OSLEVEL4', 0.0)
            summary[dataId]["OSSIGMA1"]     = ce.get('OSSIGMA1', 0.0)
            summary[dataId]["OSSIGMA2"]     = ce.get('OSSIGMA2', 0.0)
            summary[dataId]["OSSIGMA3"]     = ce.get('OSSIGMA3', 0.0)
            summary[dataId]["OSSIGMA4"]     = ce.get('OSSIGMA4', 0.0)
            #hsc = cd.get('HST', '00:00:00')
            summary[dataId]["HST"]          = ce.get('HST', 0.0)
            summary[dataId]["INSROT"]       = ce.get('INR-STR', 0.0)
            summary[dataId]["PA"]           = ce.get('INST-PA', 0.0)            
            summary[dataId]["MJD"]          = ce.get('MJD', 0.0)
            summary[dataId]["FOCUSZ"]       = ce.get('FOC-VAL', 0.0)
            summary[dataId]["ADCPOS"]       = ce.get('ADC-STR', 0.0)

            summary[dataId]['GAIN1']        = ce.get("T_GAIN1", 0.0)
            summary[dataId]['GAIN2']        = ce.get("T_GAIN2", 0.0)
            summary[dataId]['GAIN3']        = ce.get("T_GAIN3", 0.0)
            summary[dataId]['GAIN4']        = ce.get("T_GAIN4", 0.0)
            summary[dataId]['CCDTEMP']      = ce.get("DET-TMED", -1.0)

        return summary
        

    def loadCalexp(self, dataIdRegex):
        """Load the calexp data for data matching dataIdRegex.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """
        
        dataTuplesToFetch = self._regexMatchDataIds(dataIdRegex, self.dataTuples)

        
        # get the datasets corresponding to the request
        for dataTuple in dataTuplesToFetch:
            dataId = self._dataTupleToDataId(dataTuple)
            dataKey = self._dataTupleToString(dataTuple)

            #dataRef = self.butler.getDataRef(butler, dataId)
            
            if self.calexpCache.has_key(dataKey) or (dataKey in self.alreadyTriedCalexp):
                continue

            self.printStartLoad("Loading Calexp for: " + dataKey + "...")

            if self.butler.datasetExists('calexp_md', dataId):
                calexp_md = self.butler.get('calexp_md', dataId)
                
                self.wcsCache[dataKey]      = afwImage.makeWcs(calexp_md)

                raftName, ccdName = self.cameraInfo.getRaftAndSensorNames(dataId)

                self.detectorCache[dataKey] = self.cameraInfo.detectors[ccdName] #ccdDetector
                if len(raftName) > 0:
                    self.raftDetectorCache[dataKey] = self.cameraInfo.detectors[raftName]

                
                self.filterCache[dataKey]   = afwImage.Filter(calexp_md)
                self.calibCache[dataKey]    = afwImage.Calib(calexp_md)

                # try the other MAGZERO values instead
                if 'MAGZERO' in calexp_md.names():
                    exptime = calexp_md.get("EXPTIME")
                    mag0 = calexp_md.get('MAGZERO')
                    fmag0 = 10**(mag0/2.5)*exptime
                    merr0 = calexp_md.get("MAGZERO_RMS")/numpy.sqrt(calexp_md.get("MAGZERO_NOBJ"))
                    fmerr0 = fmag0*numpy.log(10.0)*merr0*0.4
                    self.calibCache[dataKey].setFluxMag0(fmag0, fmerr0)

                
                # if we have a meas_mosaic value, use that for fmag0
                # need a try block since butler will raise an exception if registries don't include tract
                try:
                    if self.butler.datasetExists("fcr", dataId) and haveMosaic:
                        fcr_md = self.butler.get("fcr_md", dataId, immediate=True)
                        ffp    = measMos.FluxFitParams(fcr_md)
                        fmag0 = fcr_md.get("FLUXMAG0")
                        #print "FMAG0=", fmag0
                        self.calibCache[dataKey].setFluxMag0(fmag0)
                except:
                    pass
                    
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
                    psf = self.butler.get("psf", visit=dataId['visit'],
                                             raft=dataId['raft'], sensor=dataId['sensor'])
                    attr = measAlg.PsfAttributes(psf, width // 2, height // 2)
                    fwhm = attr.computeGaussianWidth() * self.wcsCache[dataKey].pixelScale().asArcseconds() * sigmaToFwhm
                except Exception, e:
                    fwhm = -1.0

                if (self.calexpCache[dataKey].has_key('fwhm') and
                    numpy.isnan(float(self.calexpCache[dataKey]['fwhm']))):
                    self.calexpCache[dataKey]['fwhm'] = fwhm
                
                self.calexpQueryCache[dataKey] = True
                self.dataIdLookup[dataKey] = dataId


                
            else:
                calibFilename = self.butler.get('calexp_filename', dataId)
                self.log.log(self.log.WARN, "Skipping " + str(dataTuple) + ". Calib output file missing:")
                self.log.log(self.log.WARN, "   "+str(calibFilename))
                self.alreadyTriedCalexp.add(dataKey)

            self.printStopLoad("Calexp load for: " + dataKey)
            

            
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
    def _regexMatchDataIds(self, dataIdRegexDict, availableDataTuples, exact=True, verbose=False):
        """Match available data with regexes in a dataId dictionary        
        
        @param dataIdRegexDict dataId dict of regular expressions for data to be handled.
        @param availableDataTuples data sets available to be retrieved.
        """

        # go through the list of what's available, and compare to what we're asked for
        # Put matches in a list of tuples, eg. [(vis1,sna1,raf1,sen1),(vis2,sna2,raf2,sen2)] 
        dataTuples = []
        for dataTuple in availableDataTuples:
            if verbose:
                self.log.log(self.log.INFO, str(dataTuple))
            # start true and fail if any dataId keys fail ... eg. 'visit' doesn't match
            match = True
            for i in range(len(self.dataIdNames)):
                dataIdName = self.dataIdNames[i]   # eg. 'visit', 'sensor', etc
                regexForThisId = dataIdRegexDict.get(dataIdName, '.*') # default to '.*' or 'anything'
                dataId = dataTuple[i]

                if exact:
                    if str(regexForThisId) != str(dataId):
                        match = False
                        break
                else:
                    # if it doesn't match, this frame isn't to be run.
                    if not re.search(str(regexForThisId),  str(dataId)):
                        match = False
                        break
                
                # ignore the guiding ccds on the hsc camera
                if re.search('^hsc.*', self.cameraInfo.name) and dataIdName == 'ccd' and dataId > 103:
                    match = False
                    break

            if match:
                dataTuples.append(dataTuple)
                
        return dataTuples
                

    
