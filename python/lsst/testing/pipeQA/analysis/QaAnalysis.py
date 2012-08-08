import os
import lsst.testing.pipeQA.figures as qaFig
import numpy
import cPickle as pickle
import eups

import lsst.testing.pipeQA.TestCode as testCode

class QaAnalysis(object):
    """Baseclass for analysis classes."""

    def __init__(self, testLabel=None, useCache=False, wwwCache=True, delaySummary=False):
        """
        @param testLabel   A name for this kind of analysis test.
        """
        
        self.testSets = {}
        self.testLabel = testLabel

        # if we're not going to use the cached values
        # we'll have to clean the output directory on our first call
        self.useCache = useCache
        self.clean    = not useCache
        self.wwwCache = wwwCache
        self.delaySummary = delaySummary
        
    
    def getTestSet(self, data, dataId, label=None):
        """Get a TestSet object in the correct group.

        @param data    a QaData object
        @param dataId  a dataId dictionary
        @param label   a label for particular TestSet
        """

        dataIdStd = data.cameraInfo.dataIdCameraToStandard(dataId)
        group = dataIdStd['visit']
        
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
        tsId = str(group)+ '-' + filterName
        if data.cameraInfo.name == 'sdss':
            tsId = group
            
        if not self.testSets.has_key(tsId):
            self.testSets[tsId] = testCode.TestSet(label, group=tsId, clean=self.clean,
                                                   wwwCache=self.wwwCache)
            self.testSets[tsId].addMetadata('dataset', data.getDataName())
            self.testSets[tsId].addMetadata(tsIdLabel, tsId)
            pqaVersion = eups.getSetupVersion('testing_pipeQA')
            dqaVersion = eups.getSetupVersion('testing_displayQA')
            self.testSets[tsId].addMetadata('PipeQA', pqaVersion)
            self.testSets[tsId].addMetadata('DisplayQA', dqaVersion)
            
            
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

