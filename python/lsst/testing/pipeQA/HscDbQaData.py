import sys, os, re, copy
import time, numpy

import lsst.afw.detection               as afwDet
import lsst.afw.image                   as afwImage
import lsst.meas.algorithms             as measAlg
import lsst.afw.geom                    as afwGeom
import lsst.afw.coord                   as afwCoord
import lsst.afw.cameraGeom              as cameraGeom

import CameraInfo                       as qaCamInfo

from HscDatabaseQuery import DbInterface, DatabaseIdentity
from QaData        import QaData

from HscQaDataUtils import HscQaDataUtils
qaDataUtils = HscQaDataUtils()

import simRefObject as simRefObj
import source       as pqaSource


class Source(object):
    def __init__(self):
        self.sourceId       = None
        
class MatchedSource(object):
    def __init__(self):
        self.simRefObjectId = None
        self.u              = Source()
        self.g              = Source()
        self.r              = Source()
        self.i              = Source()
        self.z              = Source()
        self.y              = Source()

class Timer(object):

    def __init__(self):
        self.t0 = {}
        self.t = {}
        self.labels = []
        
    def start(self, label):
        if not self.t.has_key(label):
            self.labels.append(label)
            self.t[label] = 0.0
        self.t0[label] = time.time()
        
    def stop(self, label):
        self.t[label] += time.time() - self.t0[label]

    def write(self):
        for label in self.labels:
            print label, "%.2f" % (self.t[label])
            

#########################################################################
#
#
#
#########################################################################
class HscDbQaData(QaData):
    #Qa__init__(self, label, rerun, dataInfo):

    def __init__(self, database, rerun, cameraInfo):
        """
        @param database The name of the database to connect to
        @param rerun The data rerun to use
        @param cameraInfo A cameraInfo object describing the camera for these data
        """
        QaData.__init__(self, database, rerun, cameraInfo, qaDataUtils)
        self.dbId        = DatabaseIdentity(self.label)
        self.dbInterface = DbInterface(self.dbId)

        self.refStr = {'obj' : ('Obj', 'object'), 'src' : ('Src', 'source') }


        self.tableSuffix = ""
        if cameraInfo.name == 'suprimecam':
            self.tableSuffix = "_sup"


        # handle backward compatibility of database names
        keyList = []
        tname = 'frame_sourcelist' + self.tableSuffix
        sql = "select column_name from information_schema.columns where table_name = '"+tname+"';"
        results = self.dbInterface.execute(sql)
        for r in results:
            keyList.append(r[0])

        # default to new names
        self.dbAliases = {
            #"flux_Gaussian" : "instFlux",
            #"flux_ESG"      : "modelFlux",
            #'flux_gaussian' : 'flux_gaussian',
            }
        # reset to old names if new names not present
        for k,v in self.dbAliases.items():
            if k in keyList:
                self.dbAliases[k] = k
        

                
    def initCache(self):

        QaData.initCache(self)
        # need to intialize these differently than base class
        # ... Db has 'object' and 'source' matching to be cached
        self.matchListCache = { 'obj': {}, 'src': {} }
        self.matchQueryCache = { 'obj' : {}, 'src': {} }
        

    def calibFluxError(self, f, df, f0, df0):
        val = numpy.NaN 
        try:
            return (df/f + df0/f0)*f/f0
        except FloatingPointError:
            val = numpy.NaN 
        return val

    
    def getMatchListBySensor(self, dataIdRegex, useRef='src'):
        """Get a dict of all SourceMatches matching dataId, with sensor name as dict keys.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        calib          = self.getCalibBySensor(dataIdRegex)
        sourcesDict    = self.getSourceSetBySensor(dataIdRegex)
        refObjectsDict = self.getRefObjectSetBySensor(dataIdRegex)
        

        sourceLookupByRef = {}
        for k, sources in sourcesDict.items():
            for s in sources:
                sourceLookupByRef[s.get("RefId")] = s
            
        
        self.verifyDataIdKeys(dataIdRegex.keys(), raiseOnFailure=True)

        setMethods = [x for x in qaDataUtils.getSourceSetAccessors()]
        selectList = ["s."+x for x in qaDataUtils.getSourceSetDbNames(self.dbAliases)]
        selectStr  = ", ".join(selectList)

        tname = 'frame' + self.tableSuffix
        sql  = 'select sce.filter from '+tname+' as sce'
        sql += ' where '
        haveAllKeys = True


        # b/c of diff cameras, dataId keys and ScienceCcdExposure schema are have different names
        # eg. visit vs. run-field, raft vs. camcol ...
        sceNames = [
            [x[0], "sce."+x[1]]
            for x in self.cameraInfo.dataIdDbNames.items() if not re.search("snap", x[0])
            ]
        nDataId = len(sceNames)
        
        idWhereList = []
        for keyNames in sceNames:
            key, sqlName = keyNames
            if dataIdRegex.has_key(key):
                likeEqual = self._sqlLikeEqual(sqlName, dataIdRegex[key])
                idWhereList.append(likeEqual)
            else:
                haveAllKeys = False
        idWhere = " and ".join(idWhereList)

        # if there are no regexes (ie. actual wildcard expressions),
        #  we can check the cache, otherwise must run the query
        
        if not re.search("\%", idWhere) and haveAllKeys:
            dataIdCopy = copy.copy(dataIdRegex)
            dataIdCopy['snap'] = "0"
            key = self._dataIdToString(dataIdCopy, defineFully=True)
            if self.matchListCache[useRef].has_key(key):
                return {key : self.matchListCache[useRef][key]}

        # if the dataIdRegex is identical to an earlier query, we must already have all the data
        dataIdStr = self._dataIdToString(dataIdRegex, defineFully=True)
        if self.matchQueryCache[useRef].has_key(dataIdStr):
            matchListDict = {}
            # get only the ones that match the request
            for key, matchList in self.matchListCache[useRef].items():
                if re.search(dataIdStr, key):
                    matchListDict[key] = matchList
            return matchListDict

        
        sql += idWhere
        result = self.dbInterface.execute(sql)
        filterName = result[0]

        arcsecErr = 1.0/206265.0


        flagBad   = "s.flag_badctd"
        flagSat   = "s.flag_pixsttctr"
        flagIntrp = "s.flag_pixiplctr"
        flagEdge  = "s.flag_pixedg"
        flagNeg   = "s.flag_neg"


        slist = 'frame_sourcelist' + self.tableSuffix
        mlist = 'frame_matchlist' + self.tableSuffix
        ftab  = 'frame' + self.tableSuffix
        
        if True:
            sql  = 'select '+','.join(zip(*sceNames)[1])+', m.ref_flux, '
            sql += 'm.ref_ra2000, m.ref_dec2000, m.ra2000, m.dec2000, '
            sql += ", ".join([flagBad,flagSat,flagIntrp,flagEdge,flagNeg]) + ", "
            sql += 'm.classification_extendedness, m.ref_id, m.obj_id, '
            sql += selectStr
            sql += '  from '+slist+' as s, '+ftab+' as sce, '+mlist+' as m'
            #sql += '  where (sce.frame_id = m.frame_id) and (s.obj_id = m.obj_id) '
            sql += '  where (sce.frame_id = m.frame_id) and '
            sql += '        (sce.frame_id = s.frame_id) and '
            sql += '        (abs(s.ra2000 - m.ra2000) < '+str(arcsecErr)+") and "
            sql += '        (abs(s.dec2000 - m.dec2000) < '+str(arcsecErr)+") "
            sql += '    and '+idWhere
        else:
            # this will have to be updated for the different dataIdNames when non-lsst cameras get used.
            sql  = 'select '+ ", ".join(zip(*sceNames)[1])+', sro.%sMag, sro.ra, sro.decl, sro.isStar, sro.refObjectId, s.sourceId, '%(filterName)
            sql += ' rom.n%sMatches,' % (self.refStr[useRef][0])
            sql += selectStr
            sql += '  from Source as s, Science_Ccd_Exposure as sce,'
            sql += '    Ref%sMatch as rom, RefObject as sro' % (self.refStr[useRef][0])
            sql += '  where (s.scienceCcdExposureId = sce.scienceCcdExposureId)'
            sql += '    and (s.%sId = rom.%sId) and (rom.refObjectId = sro.refObjectId)' % \
                   (self.refStr[useRef][1], self.refStr[useRef][1])
            if useRef == 'obj':
                sql += '    and (s.objectID is not NULL) '
            sql += '    and '+idWhere

        self.printStartLoad("Loading MatchList ("+ self.refStr[useRef][1]  +") for: " + dataIdStr + "...")


        # run the query
        #print sql

        results  = self.dbInterface.execute(sql)

        self.sqlCache['match'][dataIdStr] = sql


        # parse results and put them in a sourceSet
        multiplicity = {}
        matchListDict = {}
        i_count = 0
        for row in results:


            nFields = 13 + nDataId
            
            refflux, ra, dec, srcRa, srcDec, \
                isBad, isSat, isIntrp, isEdge, isNeg, \
                isStar, refObjId, srcId = row[nDataId:nFields]
            mag = -2.5*numpy.log10(refflux)
            
            ra = numpy.degrees(ra)
            dec = numpy.degrees(dec)

            nMatches = 1
            dataIdTmp = {}
            for j in range(nDataId):
                idName = sceNames[j][0]
                dataIdTmp[idName] = row[j]


            key = self._dataIdToString(dataIdTmp, defineFully=True)            
            self.dataIdLookup[key] = dataIdTmp

            if not matchListDict.has_key(key):
                refCatObj = pqaSource.RefCatalog()
                refCat    = refCatObj.catalog
                catObj    = pqaSource.Catalog(qaDataUtils)
                cat       = catObj.catalog
                
                matchListDict[key] = []

                    
            matchList = matchListDict[key]

            # reference objects
            sref = refCat.addNew()

            sref.setId(refObjId)
            #sref.setId(srcId)
            sref.setD(self.k_rRa, ra)
            sref.setD(self.k_rDec, dec)

            # clip at -30
            if mag < -30:
                mag = -30
            flux = refflux #10**(-mag/2.5)

            sref.setD(self.k_rPsf, flux)
            sref.setD(self.k_rAp, flux)
            sref.setD(self.k_rMod, flux)
            sref.setD(self.k_rInst, flux)

            # sources
            s = cat.addNew()

            realId = srcId
            if sourceLookupByRef.has_key(refObjId):
                realId = sourceLookupByRef[refObjId].getId()
            s.setId(realId)
            s.setD(self.k_ext, 0.0 if isStar else 1.0)

            
            i = 0
            for value in row[nFields:]:
               if not value is None:
                    setKey = catObj.setKeys[i]
                    keyName = catObj.keyNames[i]
                    if isinstance(value, str):
                        #print ord(value)
                        value = 1 if ord(value) else 0
                    #if keyName == 'Ra':
                    #    value = numpy.degrees(value)
                    #if keyName == 'Dec':
                    #    value = numpy.degrees(value)
                    #print keyName, value, type(value)
                    if value is None:
                        value = numpy.nan
                    s.set(setKey, value)
               i += 1


            # overwrite the values we loaded into these (those assumed the sourcelist flags,
            # and we're using the icsource ones.
            s.setI(self.k_intc,  isIntrp)
            s.setI(self.k_neg,   isNeg)
            s.setI(self.k_edg,   isEdge)
            s.setI(self.k_bad,   isBad)
            s.setI(self.k_satc,  isSat)

            # calibrate it
            fmag0, fmag0Err = calib[key].getFluxMag0()

            # WARNGING - fmag0Err is wrong, it's the error on the per-second value (magzero_rms in db)
            
            # HSC database stores magzero and magzero_rms, convert to fmag0,fmag0Err
            #fmag0Err = fmag0Err*fmag0*numpy.log(10.0)/2.5
            #fmag0    = 10.0**(fmag0/2.5)


            # fluxes 
            zp   = -2.5*numpy.log10(fmag0)
            imag = -2.5*numpy.log10(s.getD(self.k_Psf))

            s.setD(self.k_Psf,   s.getD(self.k_Psf)/fmag0)
            s.setD(self.k_Ap,    s.getD(self.k_Ap)/fmag0)
            s.setD(self.k_Mod,   s.getD(self.k_Mod)/fmag0)
            s.setD(self.k_Inst,  s.getD(self.k_Inst)/fmag0)
            

            i_count += 1
                
            # flux errors
            psfFluxErr  = qaDataUtils.calibFluxError(s.getD(self.k_Psf), s.getD(self.k_PsfE),
                                                     fmag0, fmag0Err)
            s.setD(self.k_PsfE, psfFluxErr)

            apFluxErr   = qaDataUtils.calibFluxError(s.getD(self.k_Ap),  s.getD(self.k_ApE),
                                                     fmag0, fmag0Err)
            s.setD(self.k_ApE, apFluxErr)

            modFluxErr  = qaDataUtils.calibFluxError(s.getD(self.k_Mod), s.getD(self.k_ModE),
                                                     fmag0, fmag0Err)
            s.setD(self.k_ModE, modFluxErr)

            instFluxErr = qaDataUtils.calibFluxError(s.getD(self.k_Inst),  s.getD(self.k_InstE),
                                                     fmag0, fmag0Err)
            s.setD(self.k_InstE, instFluxErr)

            dist = 0.0

            matchList.append([sref, s, dist])
            
            multiplicity[s.getId()] = nMatches

        
        ######
        ######
        ######
        # Determine which are orphans, blends, straight matches, and non-detections
        typeDict = {}
        for key in matchListDict.keys():
            matchList = matchListDict[key]

            sources    = sourcesDict[key]
            
            if refObjectsDict.has_key(key) and len(refObjectsDict[key]) > 0:
                refObjects = refObjectsDict[key]
            else:
                # just use the matched objects ... doomed to succeed
                refObjects = simRefObj.SimRefObjectSet() # an empty set
                for m in matchList:
                    sref, s, dist = m
                    sro = simRefObj.SimRefObject()
                    sro.setId(sref.getId())
                    refObjects.append(sro)
                
            typeDict[key] = {}
                            
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
                        
            self.printMidLoad('\n        %s: Undet, orphan, matched, blended = %d %d %d %d' % (
                key, len(undetected), len(orphans), len(matched), len(blended))
                              )

            typeDict[key]['orphan']     = orphans
            typeDict[key]['matched']    = matched
            typeDict[key]['blended']    = blended
            typeDict[key]['undetected'] = undetected

            # cache it
            self.matchListCache[useRef][key] = typeDict[key]
            
            # Determine which are orphans, blends, straight matches, and non-detections
            ######
            ######
            ######
        
        # cache it
        self.matchQueryCache[useRef][dataIdStr] = True

        self.printStopLoad()
        
        return typeDict


    def getSourceSetBySensor(self, dataIdRegex):
        """Get a dict of all Sources matching dataId, with sensor name as dict keys.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        # verify that the dataId keys are valid
        self.verifyDataIdKeys(dataIdRegex.keys(), raiseOnFailure=True)

        setMethods = [x for x in qaDataUtils.getSourceSetAccessors()]
        selectList = ["s."+x for x in qaDataUtils.getSourceSetDbNames(self.dbAliases)]
        selectStr  = ", ".join(selectList)

        # b/c of diff cameras, dataId keys and ScienceCcdExposure schema are have different names
        # eg. visit vs. run-field, raft vs. camcol ...
        sceNames = [
            [x[0], "sce."+x[1]]
            for x in self.cameraInfo.dataIdDbNames.items() if not re.search("snap", x[0])
            ]
        
        # this will have to be updated for the different dataIdNames when non-lsst cameras get used.
        slist = 'frame_sourcelist' + self.tableSuffix
        ftab = 'frame' + self.tableSuffix
        
        sql  = 'select '+", ".join(zip(*sceNames)[1])+',s.obj_id,'+selectStr
        sql += '  from '+slist+' as s, '+ftab+' as sce'
        sql += '  where (s.frame_id = sce.frame_id)'
        haveAllKeys = True

        for keyNames in sceNames:
            key, sqlName = keyNames
            if dataIdRegex.has_key(key):
                sql += '    and '+self._sqlLikeEqual(sqlName, dataIdRegex[key])
            else:
                haveAllKeys = False
        
        # if there are no regexes (ie. actual wildcard expressions),
        #  we can check the cache, otherwise must run the query
        
        if not re.search("\%", sql) and haveAllKeys:
            dataIdCopy = copy.copy(dataIdRegex)
            dataIdCopy['snap'] = "0"
            key = self._dataIdToString(dataIdCopy, defineFully=True)
            if self.sourceSetCache.has_key(key):
                return {key : self.sourceSetCache[key]}

        # if the dataIdRegex is identical to an earlier query, we must already have all the data
        dataIdStr = self._dataIdToString(dataIdRegex, defineFully=True)
        if self.queryCache.has_key(dataIdStr):
            ssDict = {}
            # get only the ones that match the request
            for key, ss in self.sourceSetCache.items():
                if re.search(dataIdStr, key):
                    ssDict[key] = ss
            return ssDict

        self.queryCache[dataIdStr] = True
        
        self.printStartLoad("Loading SourceSets for: " + dataIdStr + "...")

        # run the query
        results  = self.dbInterface.execute(sql)
        self.sqlCache['src'][dataIdStr] = sql
        calib = self.getCalibBySensor(dataIdRegex)

        
        # parse results and put them in a sourceSet
        ssDict = {}
        for k in calib.keys():
            catObj = pqaSource.Catalog(qaDataUtils)
            ssDict[k] = catObj.catalog

            if False:
                psfKey = catObj.keyDict['PsfFlux']
                apKey  = catObj.keyDict['ApFlux']
                modKey = catObj.keyDict['ModelFlux']
                instKey = catObj.keyDict['InstFlux']

                psfErrKey = catObj.keyDict['PsfFluxErr']
                apErrKey  = catObj.keyDict['ApFluxErr']
                modErrKey = catObj.keyDict['ModelFluxErr']
                instErrKey = catObj.keyDict['InstFluxErr']

                fPixInterpCenKey = catObj.keyDict['FlagPixInterpCen']
                fNegativeKey     = catObj.keyDict['FlagNegative']
                fPixEdgeKey      = catObj.keyDict['FlagPixEdge']
                fBadCentroidKey  = catObj.keyDict['FlagBadCentroid']
                fPixSaturCenKey  = catObj.keyDict['FlagPixSaturCen']
            
        for row in results:

            # get the values for the dataId
            i = 0
            dataIdTmp = {}
            for idName, dbName in sceNames:
                dataIdTmp[idName] = row[i]
                i += 1
            sid = row[i]
            #refid = row[i+1]

            #if refid:
            #    sid = refid
            nIdKeys = i+1

            key = self._dataIdToString(dataIdTmp, defineFully=True)
            self.dataIdLookup[key] = dataIdTmp

            s = ssDict[key].addNew()
            
            s.setId(sid)
            
            i = 0
            for value in row[nIdKeys:]:
                #print value, type(value), catObj.keyNames[i]
                if not value is None:
                    setKey = catObj.setKeys[i]
                    keyName = catObj.keyNames[i]
                    if isinstance(value, str) and len(value) == 1:
                        value = 1 if ord(value) else 0
                    if keyName == 'Ra':
                        value = numpy.degrees(value)
                    if keyName == 'Dec':
                        value = numpy.degrees(value)
                    s.set(setKey, value)
                i += 1


            # calibrate it
            fmag0, fmag0Err = calib[key].getFluxMag0()

            # WARNGING - fmag0Err is wrong, it's the error on the per-second value (magzero_rms in db)
            
            # HSC database stores magzero and magzero_rms, convert to fmag0,fmag0Err
            #fmag0Err = fmag0Err*fmag0*numpy.log(10.0)/2.5
            #fmag0    = 10.0**(-fmag0/2.5)

            if (fmag0 == 0.0):
                continue

            # fluxes
            s.setD(self.k_Psf,   s.getD(self.k_Psf)/fmag0)
            s.setD(self.k_Ap,    s.getD(self.k_Ap)/fmag0)
            s.setD(self.k_Mod,   s.getD(self.k_Mod)/fmag0)
            s.setD(self.k_Inst,  s.getD(self.k_Inst)/fmag0)

            # flux errors
            psfFluxErr  = qaDataUtils.calibFluxError(s.getD(self.k_Psf), s.getD(self.k_PsfE),
                                                     fmag0, fmag0Err)
            s.setD(self.k_PsfE, psfFluxErr)

            apFluxErr   = qaDataUtils.calibFluxError(s.getD(self.k_Ap),  s.getD(self.k_ApE),
                                                     fmag0, fmag0Err)
            s.setD(self.k_ApE, apFluxErr)

            modFluxErr  = qaDataUtils.calibFluxError(s.getD(self.k_Mod), s.getD(self.k_ModE),
                                                     fmag0, fmag0Err)
            s.setD(self.k_ModE, modFluxErr)

            instFluxErr = qaDataUtils.calibFluxError(s.getD(self.k_Inst),  s.getD(self.k_InstE),
                                                     fmag0, fmag0Err)
            s.setD(self.k_InstE, instFluxErr)
                

        # cache it
        for k, ss in ssDict.items():
            self.sourceSetCache[k] = ssDict[k]
        
        self.printStopLoad()

        return ssDict


    def getDataIdsFromRegex(self, dataIdRegex):


        # b/c of diff cameras, dataId keys and ScienceCcdExposure schema are have different names
        # eg. visit vs. run-field, raft vs. camcol ...
        sceNames = [
            [x[0], "sce."+x[1]]
            for x in self.cameraInfo.dataIdDbNames.items() if not re.search("snap", x[0])
            ]
        
        haveAllKeys = True
        sqlDataId = []
        for keyNames in sceNames:
            key, sqlName = keyNames
            if dataIdRegex.has_key(key):
                sqlDataId.append(self._sqlLikeEqual(sqlName, dataIdRegex[key]))
            else:
                haveAllKeys = False
        sqlDataId = " and ".join(sqlDataId)

        ftab = 'frame' + self.tableSuffix
        
        sql  = "select "+", ".join(zip(*sceNames)[1])
        sql += "  from "+ftab+" as sce "
        sql += "  where " + sqlDataId

        dataIdList = []
        results  = self.dbInterface.execute(sql)
        nIds = len(sceNames)
        for r in results:
            dataId = {}
            i = 0
            for idName, dbName in sceNames:
                dataId[idName] = r[i]
                i += 1

            dataIdList.append(dataId)
            
        return dataIdList


    def getVisitMatchesBySensor(self, matchDatabase, matchVisit, dataIdRegex):
        """ Get a dict of all Catalog Sources matching dataId, but
        within another Science_Ccd_Exposure's polygon"""


        # If the dataIdEntry is identical to an earlier query, we must already have all the data
        dataIdStr = self._dataIdToString(dataIdRegex, defineFully=True)  # E.g. visit862826551-snap.*-raft.*-sensor.*

        if self.visitMatchQueryCache.has_key(matchDatabase):
            if self.visitMatchQueryCache[matchDatabase].has_key(matchVisit):
                vmqCache = self.visitMatchQueryCache[matchDatabase][matchVisit]

                if vmqCache.has_key(dataIdStr):
                    vmCache = self.visitMatchCache[matchDatabase][matchVisit]
                    vmDict = {}
                    for key, ss in vmCache.items():
                        if re.search(dataIdStr, key):
                            vmDict[key] = ss
                    return vmDict


        # Load each of the dataIds
        dataIdList = self.getDataIdsFromRegex(dataIdRegex)

        # Set up the outputs
        calib = self.getCalibBySensor(dataIdRegex)
        vmDict = {}
        for k in calib.keys():
            vmDict[k] = []

        for dataIdEntry in dataIdList:
            visit, raft, sensor = dataIdEntry['visit'], dataIdEntry['raft'], dataIdEntry['sensor']
            dataIdEntryStr = self._dataIdToString(dataIdEntry, defineFully=True) # E.g. visit862826551-snap0-raft30-sensor20

            haveAllKeys = True
            sqlDataId = []
            for keyNames in [['visit', 'sce.visit'], ['raft', 'sce.raftName'], ['sensor', 'sce.ccdName']]:
                key, sqlName = keyNames
                if dataIdEntry.has_key(key):
                    sqlDataId.append(self._sqlLikeEqual(sqlName, dataIdEntry[key]))
                else:
                    haveAllKeys = False
            sqlDataId = " and ".join(sqlDataId)

            # Poly comes from our own database
            sql1  = 'SELECT poly FROM Science_Ccd_Exposure as sce '
            sql1 += 'WHERE %s ' % (sqlDataId) 
            sql1 += 'INTO @poly;'

            sql2  = 'CALL scisql.scisql_s2CPolyRegion(@poly, 20);'
        
            # Selection of source matches from the comparison database
            self.verifyDataIdKeys(dataIdRegex.keys(), raiseOnFailure=True)
            setMethods = ["set"+x for x in qaDataUtils.getSourceSetAccessors()]
            selectList = ["s."+x for x in qaDataUtils.getSourceSetDbNames(self.dbAliases)]
            selectStr = ", ".join(selectList)
            sql3  = 'SELECT sce.visit, sce.raftName, sce.ccdName, sce.filterName, '                # 4 values
            sql3 += ' sce.fluxMag0, sce.fluxMag0Sigma,'                                            # 2 values
            sql3 += '   CASE WHEN sce.filterId = 0 THEN sro.uMag'
            sql3 += '        WHEN sce.filterId = 1 THEN sro.gMag'
            sql3 += '        WHEN sce.filterId = 2 THEN sro.rMag'
            sql3 += '        WHEN sce.filterId = 3 THEN sro.iMag'
            sql3 += '        WHEN sce.filterId = 4 THEN sro.zMag'
            sql3 += '        WHEN sce.filterId = 5 THEN sro.yMag'
            sql3 += '   END as mag,'                                                               # 1 value
            sql3 += ' sro.ra, sro.decl, sro.isStar, sro.refObjectId,'                              # 4 values
            sql3 += selectStr
            sql3 += ' FROM %s.Source AS s USE INDEX FOR JOIN(IDX_htmId20)' % (matchDatabase)
            sql3 += ' INNER JOIN %s.Science_Ccd_Exposure AS sce ' % (matchDatabase)
            sql3 += ' ON (s.scienceCcdExposureId = sce.scienceCcdExposureId) AND (sce.visit = %s)' % (matchVisit)
            sql3 += '   INNER JOIN %s.RefSrcMatch AS rsm ON (s.sourceId = rsm.sourceId)' % (matchDatabase)
            sql3 += '   INNER JOIN %s.RefObject AS sro ON (sro.refObjectId = rsm.refObjectId)'  % (matchDatabase) 
            sql3 += '   INNER JOIN scisql.Region AS reg ON (s.htmId20 BETWEEN reg.htmMin AND reg.htmMax) '
            sql3 += 'WHERE scisql_s2PtInCPoly(s.ra, s.decl, @poly) = 1;'

            #if not re.search("\%", sql1) and haveAllKeys:
            #    dataIdCopy = copy.copy(dataIdEntry)
            #    dataIdCopy['snap'] = "0"
            #    key = self._dataIdToString(dataIdCopy)
            #    if self.visitMatchCache.has_key(key):
            #        vmDict[key] = self.visitMatchCache[key]
            #        continue
            
            self.printStartLoad("Loading DatasetMatches for: " + dataIdEntryStr + "...")
            self.dbInterface.execute(sql1)
            self.dbInterface.execute(sql2)
            results = self.dbInterface.execute(sql3)

            self.printMidLoad("Found %d matches..." % (len(results)))

            for row in results:
                s = pqaSource.Source()
                qaDataUtils.setSourceBlobsNone(s)
                sref = pqaSource.RefSource()
                qaDataUtils.setSourceBlobsNone(sref)

                nValues = 11
                mvisit, mraft, mccd, mfilt, fmag0, fmag0Err, mag, ra, dec, isStar, refObjId = row[:nValues]
                filt = afwImage.Filter(mfilt, True)

                sref.setId(refObjId)
                sref.setRa(ra)
                sref.setDec(dec)
                flux = 10**(-mag/2.5)
                sref.setPsfFlux(flux)
                sref.setApFlux(flux)
                sref.setModelFlux(flux)
                sref.setInstFlux(flux)

                i = 0
                for value in row[nValues:]:
                    method = getattr(s, setMethods[i])
                    if not value is None:
                        method(value)
                    i += 1
           
                for sss in [s, sref]:
                    if isStar == 1:
                        sss.setFlagForDetection(sss.getFlagForDetection() | pqaSource.STAR)
                    else:
                        sss.setFlagForDetection(sss.getFlagForDetection() & ~pqaSource.STAR)
        
                # fluxes
                s.setPsfFlux(s.getPsfFlux()/fmag0)
                s.setApFlux(s.getApFlux()/fmag0)
                s.setModelFlux(s.getModelFlux()/fmag0)
                s.setInstFlux(s.getInstFlux()/fmag0)
    
                # flux errors
                psfFluxErr  = qaDataUtils.calibFluxError(s.getPsfFlux(),   s.getPsfFluxErr(),   fmag0, fmag0Err)
                s.setPsfFluxErr(psfFluxErr)
    
                apFluxErr   = qaDataUtils.calibFluxError(s.getApFlux(),    s.getApFluxErr(),    fmag0, fmag0Err)
                s.setApFluxErr(apFluxErr)
    
                modFluxErr  = qaDataUtils.calibFluxError(s.getModelFlux(), s.getModelFluxErr(), fmag0, fmag0Err)
                s.setModelFluxErr(modFluxErr)
    
                instFluxErr = qaDataUtils.calibFluxError(s.getInstFlux(),  s.getInstFluxErr(),  fmag0, fmag0Err)
                s.setInstFluxErr(instFluxErr)
    
                vm = vmDict[dataIdEntryStr]
                vm.append( [sref, s, filt] )
    
            self.printStopLoad()
    
        # cache it
        if not self.visitMatchQueryCache.has_key(matchDatabase):
            self.visitMatchQueryCache[matchDatabase] = {}
            self.visitMatchCache[matchDatabase] = {}

        if not self.visitMatchQueryCache[matchDatabase].has_key(matchVisit):
            self.visitMatchQueryCache[matchDatabase][matchVisit] = {}
            self.visitMatchCache[matchDatabase][matchVisit] = {}

        self.visitMatchQueryCache[matchDatabase][matchVisit][dataIdStr] = True
        for k, ss in vmDict.items():
            self.visitMatchCache[matchDatabase][matchVisit][k] = ss
        
        return vmDict


    def getRefObjectSetBySensor(self, dataIdRegex):
        """Get a dict of all Catalog Sources matching dataId, with sensor name as dict keys.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        # b/c of diff cameras, dataId keys and ScienceCcdExposure schema are have different names
        # eg. visit vs. run-field, raft vs. camcol ...
        sceNames = [
            [x[0], "sce."+x[1]]
            for x in self.cameraInfo.dataIdDbNames.items() if not re.search("snap", x[0])
            ]

        
        # verify that the dataId keys are valid
        self.verifyDataIdKeys(dataIdRegex.keys(), raiseOnFailure=True)

        # figure out if we have yMag
        keyList = []
        short_circuit = True
        haveYmag = False
        if not short_circuit:
            sql = "show columns from RefObject;"
            results = self.dbInterface.execute(sql)
            for r in results:
                keyList.append(r[0])
            haveYmag = 'yMag' in keyList
        
        sroFields = simRefObj.fields
        if not haveYmag:
            sroFields = [x for x in sroFields if x != 'yMag']
        sroFieldStr = ", ".join(["sro."+field for field in sroFields])

        # if the dataIdEntry is identical to an earlier query, we must already have all the data
        dataIdStr = self._dataIdToString(dataIdRegex, defineFully=True)
        if self.refObjectQueryCache.has_key(dataIdStr):
            # get only the ones that match the request
            sroDict = {}
            for key, sro in self.refObjectCache.items():
                if re.search(dataIdStr, key):
                    sroDict[key] = sro
            return sroDict


        # get a list of matching dataIds 
        dataIdList = self.getDataIdsFromRegex(dataIdRegex)
            

        # Load each of the dataIds
        sroDict = {}

        # if we have no ref object table,
        short_circuit = True
        if short_circuit:
            for dataIdEntry in dataIdList:
                key = self._dataIdToString(dataIdEntry, defineFully=True)
                if not sroDict.has_key(key):
                    sroDict[key] = simRefObj.SimRefObjectSet()
                continue
            # cache it
            for k, sro in sroDict.items():
                self.refObjectCache[k] = sroDict[k]
            return sroDict

        
        for dataIdEntry in dataIdList:

            dataIdEntryStr = self._dataIdToString(dataIdEntry, defineFully=True)            
            
            haveAllKeys = True
            sqlDataId = []
            for keyNames in sceNames:
                key, sqlName = keyNames
                if dataIdEntry.has_key(key):
                    sqlDataId.append(self._sqlLikeEqual(sqlName, dataIdEntry[key]))
                else:
                    haveAllKeys = False
            sqlDataId = " and ".join(sqlDataId)



            nStep = 2

            if nStep == 2:
                sql  = 'SELECT scisql_s2CPolyToBin('
                sql += '   sce.corner1Ra, sce.corner1Decl, '
                sql += '   sce.corner2Ra, sce.corner2Decl, '
                sql += '   sce.corner3Ra, sce.corner3Decl, '
                sql += '   sce.corner4Ra, sce.corner4Decl) '
                sql += 'FROM Science_Ccd_Exposure as sce '
                sql += 'WHERE %s ' % (sqlDataId)
                #sq += '   (sce.visit = 887252941) AND'
                #sq += '   (sce.raftName = \'2,2\') AND'
                #sq += '   (sce.ccdName = \'1,1\');'
                sql += 'INTO @poly; '

                sql2 = 'SELECT %s ' % (sroFieldStr)
                sql2 += 'FROM '
                sql2 += '    RefObject AS sro '
                sql2 += 'WHERE '
                sql2 += '    (scisql_s2PtInCPoly(sro.ra, sro.decl, @poly) = 1) '

            # use a 3 step query
            elif nStep == 3:
                sql  = 'SELECT poly FROM Science_Ccd_Exposure as sce '
                sql += 'WHERE %s ' % (sqlDataId) 
                sql += 'INTO @poly;'

                sql2 = 'CALL scisql.scisql_s2CPolyRegion(@poly, 20);'

                sql3  = 'SELECT %s ' % (sroFieldStr)
                sql3 += 'FROM RefObject AS sro INNER JOIN '
                sql3 += '   scisql.Region AS reg ON (sro.htmId20 BETWEEN reg.htmMin AND reg.htmMax) '
                sql3 += 'WHERE scisql_s2PtInCPoly(sro.ra, sro.decl, @poly) = 1;'

            
            # if there are no regexes (ie. actual wildcard expressions),
            #  we can check the cache, otherwise must run the query

            if not re.search("\%", sql) and haveAllKeys:
                dataIdCopy = copy.copy(dataIdEntry)

                key = self._dataIdToString(dataIdCopy, defineFully=True)
                if self.refObjectCache.has_key(key):
                    sroDict[key] = self.refObjectCache[key]
                    continue

                        
            self.printStartLoad("Loading RefObjects for: " + dataIdEntryStr + "...")

            # run the queries
            if nStep == 2:
                self.dbInterface.execute(sql)
                results = self.dbInterface.execute(sql2)
            elif nStep == 3:
                self.dbInterface.execute(sql)
                self.dbInterface.execute(sql2)
                results = self.dbInterface.execute(sql3)
                    

            # parse results and put them in a sourceSet
            raftName, ccdName = self.cameraInfo.getRaftAndSensorNames(dataIdEntry)
            #visit, raft, sensor = dataIdEntry['visit'], dataIdEntry['raft'], dataIdEntry['sensor']
            wcs = self.getWcsBySensor(dataIdEntry)[dataIdEntryStr]
            #raftName = self.cameraInfo.raftKeyToName(raft)
            #ccdName = self.cameraInfo.ccdKeyToName(sensor)
            bbox = self.cameraInfo.getBbox(raftName, ccdName)
            
            for row in results:
                sroStuff = list(row[:])
                if not haveYmag:
                    sroStuff.append(0.0) # dummy yMag

                # ignore things near the edge
                # ... they wouldn't be detected, and we should know about them
                ra, dec = sroStuff[2], sroStuff[3]
                x, y = wcs.skyToPixel(afwCoord.Coord(afwGeom.PointD(ra, dec)))
                if qaDataUtils.atEdge(bbox, x, y):
                   continue

                
                dataIdTmp = dataIdEntry #{'visit':str(visit), 'raft':raft, 'sensor':sensor, 'snap':'0'}
                key = self._dataIdToString(dataIdTmp, defineFully=True)
                self.dataIdLookup[key] = dataIdTmp

                if not sroDict.has_key(key):
                    sroDict[key] = simRefObj.SimRefObjectSet()
                sros = sroDict[key]
                sros.push_back(simRefObj.SimRefObject(*sroStuff))

            self.refObjectQueryCache[dataIdStr] = True
            
            self.printStopLoad()

            
        # cache it
        for k, sro in sroDict.items():
            self.refObjectCache[k] = sroDict[k]

        return sroDict



    def getVisits(self, dataIdRegex):
        """ Return explicit visits matching for a dataIdRegex.

        @param dataIdRegex dataId dict containing regular expressions of data to retrieve.
        """


        columnsTmp = zip(*self.cameraInfo.dataInfo)[0]
        visitLike = zip(*self.cameraInfo.dataInfo)[1]
        dbNames = []
        columns = []
        for i in range(len(columnsTmp)):
            c = columnsTmp[i]
            n = visitLike[i]
            if not re.search("snap", c) and n > 0:
                columns.append(self.cameraInfo.dataIdDbNames[c])
                dbNames.append([c, self.cameraInfo.dataIdDbNames[c]])

        ftab = 'frame' + self.tableSuffix
        
        sql = "select distinct "+", ".join(columns)+" from "+ftab
        sql += "   where "
        haveAllKeys = True

        whereList = []
        for keyNames in dbNames:
            key, sqlName = keyNames
            if re.search("snap", key):
                continue
            if dataIdRegex.has_key(key):
                whereList.append(self._sqlLikeEqual(sqlName, dataIdRegex[key]))
            else:
                haveAllKeys = False
        sql += " and ".join(whereList)

        results = self.dbInterface.execute(sql)

        visits = []
        for r in results:
            dataId = {}
            for i in range(len(columns)):
                dataId[columns[i]] = str(r[i])
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


        # b/c of diff cameras, dataId keys and ScienceCcdExposure schema are have different names
        # eg. visit vs. run-field, raft vs. camcol ...
        sceNames = [
            [x[0], "sce."+x[1]]
            for x in self.cameraInfo.dataIdDbNames.items() if not re.search("snap", x[0])
            ]
        

        ftab = 'frame' + self.tableSuffix
        sql = "select "+", ".join(zip(*sceNames)[1])+" from "+ftab+" as sce"
        sql += "   where "
        whereList = []
        for keyNames in sceNames:
            key, sqlName = keyNames
            if dataIdRegex.has_key(key):
                whereList.append(self._sqlLikeEqual(sqlName, dataIdRegex[key]))
        sql += " and ".join(whereList)

        results = self.dbInterface.execute(sql)

        dataIdDict = {}
        #ccdConvention = 'ccd'
        #if not dataIdRegex.has_key('ccd'):
        #    ccdConvention = 'sensor'
            
        for r in results:

            i = 0
            thisDataId = {}
            for idName, dbName in sceNames:
                thisDataId[idName] = r[i]
                i += 1
            
            if breakBy == 'raft':
                # handle lsst/hsc different naming conventions
                ccd = dataIdRegex[self.cameraInfo.dataIdTranslationMap['sensor']]

            key = self._dataIdToString(thisDataId, defineFully=True)
            dataIdDict[key] = thisDataId


        # store the list of broken dataIds 
        self.brokenDataIdList = []
        for key in sorted(dataIdDict.keys()):
            self.brokenDataIdList.append(dataIdDict[key])
        
        return copy.copy(self.brokenDataIdList)


    def getSummaryDataBySensor(self, dataIdRegex):
        """Get a dict of dict objects which contain specific summary data.
        
        @param dataIdRegex dataId dictionary with regular expressions to specify data to retrieve
        """
        calexp = self.getCalexpBySensor(dataIdRegex)
        summary = {}
        for dataId, ce in calexp.items():
            if not dataId in summary:
                summary[dataId] = {}
            summary[dataId]["DATE_OBS"]     = ce['date_obs']
            summary[dataId]["EXPTIME"]      = ce['exptime']
            summary[dataId]['RA']           = ce['ra']
            summary[dataId]['DEC']          = ce['decl']
            summary[dataId]['ALT']          = ce['elevation']
            summary[dataId]['AZ']           = ce['azimuth']
            summary[dataId]["SKYLEVEL"]     = ce['skylevel']
            summary[dataId]["ELLIPT"]       = ce['ellipt']
            summary[dataId]["ELL_PA"]       = ce['ell_pa']
            summary[dataId]["AIRMASS"]      = ce['airmass']
            summary[dataId]["FLATNESS_RMS"] = ce['flatness_rms']
            summary[dataId]["FLATNESS_PP"]  = ce['flatness_pp']
            summary[dataId]["SIGMA_SKY"]    = ce['sigma_sky']
            summary[dataId]["SEEING"]       = ce['seeing']
            summary[dataId]['OBJECT']       = ce['object']
            #summary[dataId]["OSLEVEL1"]     = ce['oslevel1']
            #summary[dataId]["OSLEVEL2"]     = ce['oslevel2']
            #summary[dataId]["OSLEVEL3"]     = ce['oslevel3']
            #summary[dataId]["OSLEVEL4"]     = ce['oslevel4']
            #summary[dataId]["OSSIGMA1"]     = ce['ossigma1']
            #summary[dataId]["OSSIGMA2"]     = ce['ossigma2']
            #summary[dataId]["OSSIGMA3"]     = ce['ossigma3']
            #summary[dataId]["OSSIGMA4"]     = ce['ossigma4']
            summary[dataId]["HST"]          = ce['hst']
            summary[dataId]["INSROT"]       = ce['insrot']
            summary[dataId]["PA"]           = ce['pa']            
            summary[dataId]["MJD"]          = ce['mjd']
            summary[dataId]["FOCUSZ"]       = ce['focusz']
            summary[dataId]["ADCPOS"]       = ce['adcpos']
            #summary[dataId][""] = ce['']
            #summary[dataId][""] = ce['']
            
        return summary


    def loadCalexp(self, dataIdRegex):
        """Load the calexp data for data matching dataIdRegex.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        # verify that the dataId keys are valid
        self.verifyDataIdKeys(dataIdRegex.keys(), raiseOnFailure=True)

        # b/c of diff cameras, dataId keys and ScienceCcdExposure schema are have different names
        # eg. visit vs. run-field, raft vs. camcol ...
        sceDataIdNames = [
            #[x[0], "sce."+x[1]]
            x for x in self.cameraInfo.dataIdDbNames.items() if not re.search("snap", x[0])
            ]

        
        selectList = ["sce."+x for x in qaDataUtils.getSceDbNames(sceDataIdNames)]
        selectStr = ", ".join(selectList)

        # NOTE: HSC places focusz and adcpos in the Exposure table (not Frame)
        # So, getSceFoo mechanism is augmented slightly here to get the extra two values
        
        ftab = 'frame' + self.tableSuffix
        etab = 'exposure' + self.tableSuffix
        sql  = 'select '+selectStr
        sql += ', exp.focusz, exp.adcpos'
        sql += '  from '+ftab+' as sce, '+etab+' as exp'
        sql += '  where '

        haveAllKeys = True

        whereList = []
        for keyNames in [[x[0], "sce."+x[1]] for x in sceDataIdNames]:
            key, sqlName = keyNames
            if dataIdRegex.has_key(key):
                # Hsc DB stores ccd as int ... need to cast
                if sqlName == 'sce.ccd' and re.search("%", dataIdRegex[key]):
                    sqlName = 'sce.ccd::character'
                whereList.append(self._sqlLikeEqual(sqlName, dataIdRegex[key]))
            else:
                haveAllKeys = False
        sql += " (sce.exp_id = exp.exp_id) and "
        sql += " and ".join(whereList)


        # if there are no regexes (ie. actual wildcard expressions),
        #  we can check the cache, otherwise must run the query
        if not re.search("\%", sql) and haveAllKeys:
            dataIdCopy = copy.copy(dataIdRegex)
            dataIdCopy['snap'] = "0"
            key = self._dataIdToString(dataIdCopy, defineFully=True)
            if self.calexpQueryCache.has_key(key):
                return

        # if the dataIdRegex is identical to an earlier query, we must already have all the data
        dataIdStr = self._dataIdToString(dataIdRegex, defineFully=True)
        if self.calexpQueryCache.has_key(dataIdStr) and self.calexpQueryCache[dataIdStr]:
            return

        self.printStartLoad("Loading Calexp for: " + dataIdStr + "...")

        # run the query
        results  = self.dbInterface.execute(sql)

        for row in results:

            rowNames = qaDataUtils.getSceDbNames(sceDataIdNames) + ('focusz', 'adcpos')
            rowDict = dict(zip(rowNames, row))

            dataIdTmp = {}
            for idName, dbName in sceDataIdNames:
                dataIdTmp[idName] = rowDict[dbName]
                
            #visit, raft, sensor = rowDict['visit'], rowDict['raftName'], rowDict['ccdName']
            #dataIdTmp = {'visit':visit, 'raft':raft, 'sensor':sensor, 'snap':'0'}
            key = self._dataIdToString(dataIdTmp, defineFully=True)
            self.dataIdLookup[key] = dataIdTmp
            
            #print rowDict
            if not self.wcsCache.has_key(key) and rowDict.has_key('crval1'):
                crval = afwCoord.Coord(afwGeom.PointD(rowDict['crval1'], rowDict['crval2']))
                crpix = afwGeom.PointD(rowDict['crpix1'], rowDict['crpix2'])
                cd11, cd12, cd21, cd22 = rowDict['cd1_1'], rowDict['cd1_2'], rowDict['cd2_1'], rowDict['cd2_2']
                wcs = afwImage.makeWcs(crval, crpix, cd11, cd12, cd21, cd22)
                self.wcsCache[key] = wcs

            else:
                cd11, cd12, cd21, cd22 = 1.0, 0.0, 0.0, 1.0
                crval = afwCoord.Coord(afwGeom.PointD(0.0, 0.0))
                crpix = afwGeom.PointD(0.0, 0.0)
                wcs = afwImage.makeWcs(crval, crpix, cd11, cd12, cd21, cd22)
                self.wcsCache[key] = wcs

            if not self.detectorCache.has_key(key):
                raftName, ccdName = self.cameraInfo.getRaftAndSensorNames(dataIdTmp)
                self.detectorCache[key] = self.cameraInfo.detectors[ccdName] #ccdDetector
                if not raftName is None and len(raftName) > 0:
                    self.raftDetectorCache[key] = self.cameraInfo.detectors[raftName]

            if not self.filterCache.has_key(key):
                filt = afwImage.Filter(rowDict['filter'], True)
                self.filterCache[key] = filt
            
            if not self.calibCache.has_key(key):
                calib = afwImage.Calib()
                useMagzero = False
                if useMagzero:
                    mag0 = rowDict['magzero']
                    exptime = rowDict['exptime']
                    fmag0 = exptime*10**(mag0/2.5)
                else:
                    zp = rowDict['zeropt']
                    fmag0 = 10**(zp/2.5)
                fmag0Err =  rowDict['magzero_rms']*fmag0*numpy.log(10.0)/2.5
                #print mag0, exptime, fmag0
                
                calib.setFluxMag0(fmag0, fmag0Err)
                self.calibCache[key] = calib
                
            self.calexpCache[key] = rowDict
            self.calexpQueryCache[key] = True

        self.calexpQueryCache[dataIdStr] = True
        
        self.printStopLoad()



    def getCalexpEntryBySensor(self, cache, dataIdRegex):
        """Fill and return the dict for a specified calexp cache.

        @param cache The cache dictionary to return
        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        # get the datasets corresponding to the request
        self.loadCalexp(dataIdRegex)
        dataIdStr = self._dataIdToString(dataIdRegex, defineFully=True)
        entryDict = {}
        for dataKey in cache.keys():
            if re.search(dataIdStr, dataKey):
                entryDict[dataKey] = cache[dataKey]
        return entryDict



    def getSourceSet(self, dataIdRegex):
        """Get a SourceSet of all Sources matching dataId.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """

        ssDict = self.getSourceSetBySensor(dataIdRegex)
        ssReturn = []
        for key, ss in ssDict.items():
            ssReturn += ss
            
        return ssReturn


    ########################################
    #
    # utility method to generate a string for an sql 'where' clause
    # given 'field' and 'regex':
    #    - if regex is just a value - return "field = regex"
    #    - if regex is a regex - sub '.*' and '?' to '%' and return "field like 'regex'"
    ########################################
    def _sqlLikeEqual(self, field, regex):
        """Utility to convert a dataId regex to an sql 'where' clause.
        """

        regex = str(regex).strip()

        
        clause = ""
        # just a number
        if re.search('^\d+$', regex):
            if re.search('ccd', field):
                clause += field + " = '%s'" % (regex)
            else:
                clause += field + " = %s" % (regex)
                
        # comma-sep numbers
        elif re.search('^[\d,]+$', regex):
            clause += field + " = '%s'" % (regex)
        # if it's a filter
        elif re.search('^[ugrizy]$', regex):
            clause += field + " = '%s'" % (regex)            
        # .*  ?  % followed/preceeding by comma-sep numbers
        elif re.search('^(\.\*|\?|\%)[\d,]*$', regex) or re.search('^[\d,]*(\.\*|\?|\%)$', regex):
            regexSql = re.sub("(\.\*|\?)", "%", regex)
            clause += field + " like '%s'" % (regexSql)
        else:
            raise Exception("Regex for SQL can only use '.*', '?', or '%' at beginning or end of string ("+
                            field+"="+regex+")")

        return "("+clause+")"



