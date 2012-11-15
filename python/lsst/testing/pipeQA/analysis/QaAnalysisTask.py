import os
import numpy
import cPickle as pickle
import eups

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

import lsst.testing.pipeQA.TestCode as testCode
import lsst.testing.pipeQA.figures as qaFig


class QaAnalysisConfig(pexConfig.Config):
    pass


class QaAnalysisTask(pipeBase.Task):
    """Baseclass for analysis classes."""
    ConfigClass  = QaAnalysisConfig
    _DefaultName = "qaAnalysis"


    def __init__(self, testLabel=None, useCache=False, wwwCache=True, summaryProcessing="delay",
                 lazyPlot='sensor', *args, **kwargs):
        """
        @param testLabel   A name for this kind of analysis test.
        """
        pipeBase.Task.__init__(self, *args, **kwargs)

        self.testSets  = {}
        self.testLabel = testLabel

        # if we're not going to use the cached values
        # we'll have to clean the output directory on our first call
        self.useCache     = useCache
        self.clean        = not useCache
        self.wwwCache     = wwwCache
        self.summaryProcessing = summaryProcessing

        self.summOpt = {
            'delay'     : 'delay',  # make summary figures after the final CCD is processed
            'summOnly'  : 'summOnly',   # only make summary figures
            'none'      : 'none'    # don't make summary figures
            }
        
        #self.delaySummary = delaySummary
        #self.noSummary    = noSummary
        #self.onlySummary  = onlySummary

        options = ['none', 'sensor', 'all']
        if not lazyPlot in options:
            raise ValueError, "lazyPlot must be: "+ ",".join(options) + " You said: "+lazyPlot
        
        self.lazyPlot  = lazyPlot


    def combineOutputs(self, data, dataId, label=None):
        ts = self.getTestSet(data, dataId, label, noSuffix=True)
        ts.accrete()
        #pass

        
    def getTestSet(self, data, dataId, label=None, noSuffix=False):
        """Get a TestSet object in the correct group.

        @param data    a QaData object
        @param dataId  a dataId dictionary
        @param label   a label for particular TestSet
        """

        dataIdStd = data.cameraInfo.dataIdCameraToStandard(dataId)
        group = dataIdStd['visit']


        raftName, ccdName = data.cameraInfo.getRaftAndSensorNames(dataId)
            
        filter = data.getFilterBySensor(dataId)
        # all sensors have the same filter, so just grab one
        key = filter.keys()[0]
        filterName = '?'
        if not filter[key] is None:
            filterName = filter[key].getName()

        if not label is None:
            label = self.__class__.__name__ + "."+ label
        else:
            label = self.__class__.__name__

        tsIdLabel = "visit-filter"
        tsId = str(group) + '-' + filterName
        if data.cameraInfo.name == 'sdss':
            tsId = group

        if noSuffix:
            return testCode.TestSet(label, group=tsId, clean=self.clean, wwwCache=self.wwwCache, sqliteSuffix="")
        else:
            sqliteSuffix = ccdName

            
        if not self.testSets.has_key(tsId):
            self.testSets[tsId] = testCode.TestSet(label, group=tsId, clean=self.clean,
                                                   wwwCache=self.wwwCache, sqliteSuffix=sqliteSuffix)
            
            self.testSets[tsId].addMetadata('dataset', data.getDataName())
            self.testSets[tsId].addMetadata(tsIdLabel, tsId)
            pqaVersion = eups.getSetupVersion('testing_pipeQA')
            dqaVersion = eups.getSetupVersion('testing_displayQA')
            self.testSets[tsId].addMetadata('PipeQA', pqaVersion)
            self.testSets[tsId].addMetadata('DisplayQA', dqaVersion)
            
            if hasattr(data, 'coaddTable') and not data.coaddTable is None:
                self.testSets[tsId].addMetadata('coaddTable', data.coaddTable)
            if hasattr(data, 'useForced'):
                self.testSets[tsId].addMetadata('forced', "True" if data.useForced else "False")

            key = data._dataIdToString(dataId, defineFully=True)
            sqlCache = data.sqlCache['match'].get(key, "")
            if sqlCache:
                self.testSets[tsId].addMetadata("SQL match", sqlCache)
            sqlCache = data.sqlCache['src'].get(key, "")
            if sqlCache:
                self.testSets[tsId].addMetadata("SQL src" ,  sqlCache)
                
        return self.testSets[tsId]


    def __str__(self):
        testLabel = ""
        if not self.testLabel is None:
            testLabel = "."+self.testLabel
        return self.__class__.__name__ + testLabel
    

    ##########################################
    # pure virtual methods
    #########################################
    def free(self):
        """Method to free attributes to minimize memory consumption."""
        pass
    
    def test(self):
        """Method to perform tests. """
        return []

    def plot(self):
        """Method to make figures."""
        return []

