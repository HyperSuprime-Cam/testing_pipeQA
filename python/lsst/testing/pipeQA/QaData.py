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

import sys, os, re, copy, time
import numpy

import source as pqaSource

import lsst.pex.logging as pexLog

#######################################################################
#
#
#
#######################################################################
class QaData(object):
    """Base class for QA data retrieval."""

    #######################################################################
    #
    #######################################################################
    def __init__(self, label, rerun, cameraInfo, qaDataUtils, **kwargs):
        """
        @param label The name of this data set
        @param rerun The rerun to retrieve
        @param cameraInfo A cameraInfo object containing specs on the camera
        """
        
        self.label = label
        self.rerun = rerun
        self.cameraInfo = cameraInfo
        self.dataInfo = self.cameraInfo.dataInfo

        self.log = kwargs.get('log', pexLog.getDefaultLog())
        
        self.dataIdNames   = []
        self.dataIdDiscrim = []
        for array in self.dataInfo:
            dataIdName, dataIdDiscrim = array
            self.dataIdNames.append(dataIdName)
            self.dataIdDiscrim.append(dataIdDiscrim)


        self.initCache()

        self.loadDepth = 0
        self.lastPrint = None
        self.t0 = []

        self.brokenDataIdList = []


        self.ccdConvention = 'ccd'
        if self.cameraInfo.name == 'lsstSim':
            self.ccdConvention = 'sensor'
        if self.cameraInfo.name == 'sdss':
            self.ccdConvention = 'camcol'

            
                
        refCatObj     = pqaSource.RefCatalog()
        self.k_rPsf   = refCatObj.keyDict['PsfFlux']
        self.k_rAp    = refCatObj.keyDict['ApFlux']
        self.k_rMod   = refCatObj.keyDict['ModelFlux']
        self.k_rInst  = refCatObj.keyDict['InstFlux']
        self.k_rPsfE  = refCatObj.keyDict['PsfFluxErr']
        self.k_rApE   = refCatObj.keyDict['ApFluxErr']
        self.k_rModE  = refCatObj.keyDict['ModelFluxErr']
        self.k_rInstE = refCatObj.keyDict['InstFluxErr']
        self.k_rRa    = refCatObj.keyDict['Ra']
        self.k_rDec   = refCatObj.keyDict['Dec']
        
        catObj        = pqaSource.Catalog(qaDataUtils=qaDataUtils)
        self.k_x      = catObj.keyDict['XAstrom']
        self.k_y      = catObj.keyDict['YAstrom']
        self.k_Ra     = catObj.keyDict['Ra']
        self.k_Dec    = catObj.keyDict['Dec']
        
        self.k_Psf    = catObj.keyDict['PsfFlux']
        self.k_Ap     = catObj.keyDict['ApFlux']
        self.k_Mod    = catObj.keyDict['ModelFlux']
        self.k_Inst   = catObj.keyDict['InstFlux']

        self.k_PsfE   = catObj.keyDict['PsfFluxErr']
        self.k_ApE    = catObj.keyDict['ApFluxErr']
        self.k_ModE   = catObj.keyDict['ModelFluxErr']
        self.k_InstE  = catObj.keyDict['InstFluxErr']

        self.k_ext    = catObj.keyDict['Extendedness']

        self.k_intc   = catObj.keyDict["FlagPixInterpCen"]
        self.k_satc   = catObj.keyDict["FlagPixSaturCen"]
        self.k_edg    = catObj.keyDict["FlagPixEdge"]
        self.k_neg    = catObj.keyDict['FlagNegative']
        self.k_bad    = catObj.keyDict["FlagBadCentroid"]

        self.k_ixx    = catObj.keyDict['Ixx']
        self.k_iyy    = catObj.keyDict['Iyy']
        self.k_ixy    = catObj.keyDict['Ixy']
                
        self.k_nchild  = catObj.keyDict['deblend_nchild']
            

    def isFlagged(self, src):
        intcen = src.getI(self.k_intc)
        satcen = src.getI(self.k_satc)
        edge   = src.getI(self.k_edg)
        nchild = src.getI(self.k_nchild) > 0
        return  intcen or satcen or edge or nchild
        
    def printStartLoad(self, message):

        self.loadStr = ""
        if self.loadDepth > 0:
            self.loadStr += " "*4*self.loadDepth
            self.loadStr += message
        else:
            self.loadStr += message

        self.loadDepth += 1
        t0 = time.time()
        self.t0.append(t0)
        self.log.log(self.log.INFO, self.loadStr)

        
    def printStopLoad(self, message=""):
        t0 = self.t0[-1]
        self.t0 = self.t0[:-1]
        t_final = time.time()
        t_elapsed = t_final - t0
        done =  message + " ... done (%.2fs)." % t_elapsed
        self.loadStr = ""
        if self.loadDepth > 1:
            self.loadStr += " "*4*(self.loadDepth-1)
        self.loadStr += done

        self.log.log(self.log.INFO, self.loadStr)

        self.loadDepth -= 1
        

    def initCache(self):
        """Initialize all internal cache attributes. """

        # cache the dataId requests
        # they may contain regexes ... we won't know if our cached sourceSets have
        #   all the entries that match unless we redo the query
        #   But, if we've already done the identical query, we know we have everything
        self.queryCache = {}
        self.columnQueryCache = {}

        # cache source sets to avoid reloading the same thing
        self.sourceSetCache = {}
        self.sourceSetColumnCache = {}

        self.matchQueryCache = {}
        self.matchListCache = {}

        self.refObjectQueryCache = {}
        self.refObjectCache = {}

        self.visitMatchQueryCache = {}
        self.visitMatchCache = {}

        # cache calexp info, but not the MaskedImage ... it's too big.
        self.calexpQueryCache = {}
        self.calexpCache = {}
        self.wcsCache = {}
        self.detectorCache = {}
        self.raftDetectorCache = {}
        self.filterCache = {}
        self.calibCache = {}
        self.sqlCache = {"match": {}, "src": {}}
        
        self.performCache = {}
        
        # store the explicit dataId (ie. no regexes) for each key used in a cache
        self.dataIdLookup = {}

        self.cacheList = {
            "query"          : self.queryCache,  
            "columnQuery"    : self.columnQueryCache,
            "sourceSet"      : self.sourceSetCache,
            "sourceSetColumn"  : self.sourceSetColumnCache,
            "matchQuery"     : self.matchQueryCache,
            "matchList"      : self.matchListCache,
            "refObjectQuery" : self.refObjectQueryCache,
            "refObject"      : self.refObjectCache,
            "visitMatchQuery": self.visitMatchQueryCache,
            "calexpQuery"    : self.calexpQueryCache,
            "calexp"         : self.calexpCache, 
            "wcs"            : self.wcsCache,    
            "detector"       : self.detectorCache,
            "raftDetector"   : self.raftDetectorCache,
            "filter"         : self.filterCache, 
            "calib"          : self.calibCache,  
            "dataIdLookup"   : self.dataIdLookup,
            "sql"            : self.sqlCache,
            }


    def cachePerformance(self, dataIdStr, test, label, value):
        if isinstance(dataIdStr, dict):
            dataIdStr = self._dataIdToString(dataIdStr, defineFully=True)
            
        if not self.performCache.has_key(dataIdStr):
            self.performCache[dataIdStr] = {}
        if not self.performCache[dataIdStr].has_key(test):
            self.performCache[dataIdStr][test] = {}
        if not self.performCache[dataIdStr].has_key('total'):
            self.performCache[dataIdStr]['total'] = {}
            
        if not self.performCache[dataIdStr][test].has_key(label):
            self.performCache[dataIdStr][test][label] = 0.0
        if not self.performCache[dataIdStr]['total'].has_key(label):
            self.performCache[dataIdStr]['total'][label] = 0.0
            
        self.performCache[dataIdStr][test][label]    += value
        self.performCache[dataIdStr]['total'][label] += value

        
    def getPerformance(self, dataIdStr, test, label):
        if isinstance(dataIdStr, dict):
            dataIdStr = self._dataIdToString(dataIdStr, defineFully=True)
        if self.performCache.has_key(dataIdStr):
            if self.performCache[dataIdStr].has_key(test):
                if self.performCache[dataIdStr][test].has_key(label):
                    return self.performCache[dataIdStr][test][label]
        return None

    def clearCache(self):
        """Reset all internal cache attributes."""
        for cache in self.cacheList.values():
            for key in cache.keys():
                if isinstance(cache[key], dict):
                    for key2 in cache[key].keys():
                        del cache[key][key2]
                else:
                    del cache[key]
        self.initCache()

    def printCache(self):
        for name, cache in self.__dict__.items():
            if re.search("^_", name):
                continue
            n = 0
            if isinstance(cache, dict):
                for key in cache:
                    if hasattr(cache[key], "__len__"):
                        n += len(cache[key])
            if isinstance(cache, list):
                n = len(cache)
            print name, n
                

    def getDataName(self):
        """Get a string representation of ourself."""
        # should this be __str__ or __repr__ ?
        return self.label+" rerun="+str(self.rerun)
        


    def getSourceSetColumnsBySensor(self, dataIdRegex, accessors):
        """Get a specified Source field from all sources in SourceSet as a numpy array.
        
        @param dataIdRegex dataId dict with regular expressions for data to retrieve
        @param accessors  List of accessor method names (as string without 'get' prepended)
        """

        dataIdStr = self._dataIdToString(dataIdRegex, defineFully=True)

        # if they want a specific dataId and we have it ....
        ssTDict = {}
        if self.sourceSetColumnCache.has_key(dataIdStr):
            ssTDict[dataIdStr] = {}
            haveIt = True
            for accessor in accessors:
                if self.sourceSetColumnCache[dataIdStr].has_key(accessor):
                    ssTDict[dataIdStr][accessor] = self.sourceSetColumnCache[dataIdStr][accessor]
                else:
                    haveIt = False
                    break
            if haveIt:
                return ssTDict

        # get it and put the transpose in the cache
        # return a copy of what we cached.
        ssDict = self.getSourceSetBySensor(dataIdRegex)
        ssTDict = {}
        for k, ss in ssDict.items():
            if not self.sourceSetColumnCache.has_key(k):
                self.sourceSetColumnCache[k] = {}
            ssTDict[k] = {}
            for accessor in accessors:
                tmp = numpy.array([])
                for s in ss:
                    method = getattr(s, "get"+accessor)
                    value = method()
                    tmp = numpy.append(tmp, value)
                self.sourceSetColumnCache[k][accessor] = tmp
                ssTDict[k][accessor] = tmp
        
        return ssTDict




    def getWcsBySensor(self, dataIdRegex):
        """Get a dict of Wcs objects with sensor ids as keys.
        
        @param dataIdRegex dataId dictionary with regular expressions to specify data to retrieve
        """
        return self.getCalexpEntryBySensor(self.wcsCache, dataIdRegex)
    def getDetectorBySensor(self, dataIdRegex):
        """Get a dict of Detector objects with sensor ids as keys.
        
        @param dataIdRegex dataId dictionary with regular expressions to specify data to retrieve
        """
        return self.getCalexpEntryBySensor(self.detectorCache, dataIdRegex)
    def getFilterBySensor(self, dataIdRegex):
        """Get a dict of Filter objects with sensor ids as keys.
        
        @param dataIdRegex dataId dictionary with regular expressions to specify data to retrieve
        """
        return self.getCalexpEntryBySensor(self.filterCache, dataIdRegex)
    def getCalibBySensor(self, dataIdRegex):
        """Get a dict of Calib objects with sensor ids as keys.
        
        @param dataIdRegex dataId dictionary with regular expressions to specify data to retrieve
        """
        return self.getCalexpEntryBySensor(self.calibCache, dataIdRegex)
    def getCalexpBySensor(self, dataIdRegex):
        """Get a dict of Calib objects with sensor ids as keys.
        
        @param dataIdRegex dataId dictionary with regular expressions to specify data to retrieve
        """
        return self.getCalexpEntryBySensor(self.calexpCache, dataIdRegex)
    

    def verifyDataIdKeys(self, dataIdKeys, raiseOnFailure=True):
        """Verify that what we're asked for makes sense for this camera (ie. ccd vs. sensor).

        @param dataIdRegex dataId dictionary with regular expressions to specify data to retrieve
        @param raiseOnFailure Raise an exception if verification fails.
        """
        missingKeys = []
        for key in dataIdKeys:
            if not key in self.dataIdNames:
                missingKeys.append(key)
        if len(missingKeys) > 0:
            keyStr = ",".join(missingKeys)
            if raiseOnFailure:
                raise Exception("Key(s): "+keyStr+" not valid for camera "+str(self.cameraInfo))
            return False
        return True



    def getSourceClusters(self, dataIdRegex,
                          epsilonArcsec = 0.5,
                          minNeighbors = 1,
                          pointsPerLeaf = 100,
                          leafExtentThresholdArcsec = 0.5):
        import lsst.ap.cluster as apCluster

        """Get apClusters for all Sources matching dataIdRegex.

        @param epsilonArcsec Matching distance
        @param minNeighbors Fewest neighbors to accept
        @param pointsPerLeaf who knows?
        @param leafExtentThresholdArcsec Drawing a blank here too.
        """
        import lsst.ap.cluster as apCluster
        
        policy = pexPolicy.Policy()
        policy.set("epsilonArcsec", epsilonArcsec)
        policy.set("minNeighbors", minNeighbors)
        policy.set("pointsPerLeaf", pointsPerLeaf)
        policy.set("leafExtentThresholdArcsec", leafExtentThresholdArcsec)

        sources = self.getSourceSet(dataIdRegex)
        sourceClusters = apCluster.cluster(sources, policy)
        
        return sourceClusters

    #######################################################################
    #
    #######################################################################
    def _dataTupleToString(self, dataTuple):
        """Represent a dataTuple as a string.

        @param dataTuple The dataTupe to be converted.
        """
        
        s = []
        for i in xrange(len(self.dataIdNames)):
            name = self.dataIdNames[i]
            value = re.sub("[,]", "", str(dataTuple[i]))
            s.append(name + value)
        return "-".join(s)

    #######################################################################
    # utility to convert a data tuple to a dictionary using dataId keys
    #######################################################################
    def _dataTupleToDataId(self, dataTuple):
        """Convert a dataTuple to a dataId dict.

        @param dataTuple The dataTuple to be converted.
        """
        
        dataId = {}
        for i in xrange(len(self.dataIdNames)):
            dataIdName = self.dataIdNames[i]
            dataId[dataIdName] = dataTuple[i]
        return dataId

    #######################################################################
    # utility to convert a dataId dictionary to a tuple
    #######################################################################
    def _dataIdToDataTuple(self, dataId):
        """Convert a dataId to a dataTuple

        @param dataId The dataId to be converted.
        """
        
        # if snap isn't specified, we'll add it.
        dataIdCopy = copy.copy(dataId)
        if not dataIdCopy.has_key('snap'):
            dataIdCopy['snap'] = '0'
            
        dataList = []
        for i in xrange(len(self.dataIdNames)):
            dataIdName = self.dataIdNames[i]
            if dataIdCopy.has_key(dataIdName):
                dataList.append(dataIdCopy[dataIdName])
            else:
                raise Exception("key: "+dataIdName+" not in dataId: "+ str(dataIdCopy))
        return tuple(dataList)


    def _dataIdToString(self, dataId, defineFully=False):
        """Convert a dataId dict to a string.

        @param dataId The dataId to be converted
        """

        s = []
        for i in xrange(len(self.dataIdNames)):
            dataIdName = self.dataIdNames[i]
            if dataId.has_key(dataIdName):
                s.append( dataIdName + re.sub("[,]", "", str(dataId[dataIdName])))
            elif defineFully:
                if dataIdName == 'snap':
                    s.append(dataIdName + "0")
                else:
                    s.append( dataIdName + ".*" )
        s = "-".join(s)

        return s

    
    def reduceAvailableDataTupleList(self, dataIdRegexDict):
        pass

    def keyStringsToVisits(self, keyList):
        """Extract a list of visits from a list of sensor keys.

        @param keyList List of keys to have visits extracted from.
        """
        
        visits = {}
        for key in keyList:
            visit = str(self.dataIdLookup[key]['visit'])
            visits[visit] = True
        return visits.keys()

    
    #########################################################
    # pure virtual methods
        
    def getSourceSet(self, dataIdRegex):
        """Get a SourceSet of all Sources matching dataId.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """
        raise NotImplementedError, "Must define getSourceSet in derived QaData class."
    def getSourceSetBySensor(self, dataIdRegex):
        """Get a dict of all Sources matching dataId, with sensor name as dict keys.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """
        raise NotImplementedError, "Must define getSourceSetBySensor in derived QaData class."

    def getMatchListBySensor(self, dataIdRegex):
        """Get a dict of all SourceMatches matching dataId, with sensor name as dict keys.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """
        raise NotImplementedError, "Must define getMatchListBySensor in derived QaData class."

    def getCalexpEntryBySensor(self, cache, dataIdRegex):
        """Fill and return the dict for a specified calexp cache.

        @param cache The cache dictionary to return
        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """
        raise NotImplementedError, "Must define getCalexpEntryBySensor in derived QaData class."

    def getVisits(self, dataIdRegex):
        """Get a list of all visits matching dataIdRegex

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """
        raise NotImplementedError, "Must define getVisits in derived QaData class."

        

    def getSummaryDataBySensor(self, dataIdRegex):
        """Get a dict of dict objects which contain specific summary data.
        
        @param dataIdRegex dataId dictionary with regular expressions to specify data to retrieve
        """
        raise NotImplementedError("Please define getSummaryDataBySensor() in the derived class.")
    
    
            

    def loadCalexp(self, dataIdRegex):
        """Load the calexp data for data matching dataIdRegex.

        @param dataIdRegex dataId dict of regular expressions for data to be retrieved
        """
        raise NotImplementedError, "Must define loadCalexp in derived QaData class."

    def breakDataId(self, dataId, breakBy):
        """Take a dataId with regexes and return a list of dataId regexes
        which break the dataId by raft, or ccd.

        @param breakBy   'visit', 'raft', or 'ccd'
        """
        raise NotImplementedError, "Must define breakDataId in derived QaData class."
