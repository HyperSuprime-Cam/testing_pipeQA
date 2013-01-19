import os, re
import numpy
import cPickle as pickle
import eups
import datetime

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


    def __init__(self, testLabel=None, useCache=False, wwwCache=True, delaySummary=False,
                 lazyPlot='sensor', *args, **kwargs):
        """
        @param testLabel   A name for this kind of analysis test.
        """
        pipeBase.Task.__init__(self, *args, **kwargs)

        self.testSets  = {}
        self.testLabel = testLabel

        # if we're not going to use the cached values
        # we'll have to clean the output directory on our first call
        self.useCache = useCache
        self.clean    = not useCache
        self.wwwCache = wwwCache
        self.delaySummary = delaySummary

        options = ['none', 'sensor', 'all']
        if not lazyPlot in options:
            raise ValueError, "lazyPlot must be: "+ ",".join(options) + " You said: "+lazyPlot
        
        self.lazyPlot  = lazyPlot

        
    def getTestSet(self, data, dataId, label=None):
        """Get a TestSet object in the correct group.

        @param data    a QaData object
        @param dataId  a dataId dictionary
        @param label   a label for particular TestSet
        """

        dataIdStd = data.cameraInfo.dataIdCameraToStandard(dataId)
        group = dataIdStd['visit']

        raftName, ccdName = data.cameraInfo.getRaftAndSensorNames(dataId)
        ccdName = data.cameraInfo.getDetectorName(raftName, ccdName)
        
        filter = data.getFilterBySensor(dataId)
        # all sensors have the same filter, so just grab one
        key = filter.keys()[0]
        filterName = '?'
        if not filter[key] is None:
            filterName = filter[key].getName()

        summaryInfo = data.getSummaryDataBySensor(dataId)[key]

        dateObs = summaryInfo.get('DATE_OBS', "unknown")
        if dateObs:
            dateObs = dateObs.strftime("%Y-%m-%d")
        dateObs += ' (MJD: %.6f)' % (summaryInfo.get("MJD", 0.0))
        
        expTime = summaryInfo.get('EXPTIME', 'unknown')
        raDec   = summaryInfo.get('RA', 'unk') + " " + summaryInfo.get('DEC', 'unk')
        alt = summaryInfo.get('ALT', None)
        az  = str(summaryInfo.get('AZ', None))
        if alt:
            altAz   = str(alt) + " " + az + " (airmass: %.2f)" % (summaryInfo.get("AIRMASS", "unk"))
        awbject = summaryInfo.get('OBJECT', None)
        hst = summaryInfo.get("HST", None)
        if hst:
            hst = hst.strftime("%H:%M:%S")
        else:
            hst = 'unk'
        
        insrot = summaryInfo.get("INSROT", 'unk')
        #focusz = summaryInfo.get("FOCUSZ", 'unk')
        inspa = summaryInfo.get("PA", 'unk')

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

            self.testSets[tsId].addMetadata('DATE_OBS', dateObs)
            self.testSets[tsId].addMetadata('HST', hst)
            self.testSets[tsId].addMetadata('EXPTIME', expTime)
            self.testSets[tsId].addMetadata('RaDec', raDec)
            self.testSets[tsId].addMetadata('AltAz', altAz)
            if awbject:
                self.testSets[tsId].addMetadata('Object', awbject)
            if insrot:
                self.testSets[tsId].addMetadata('Inst.Rot.', insrot)
            if inspa:
                self.testSets[tsId].addMetadata('Inst.PA', inspa)
                

            # version info
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
                self.testSets[tsId].addMetadata("SQL_match-"+ccdName, sqlCache)
            sqlCache = data.sqlCache['src'].get(key, "")
            if sqlCache:
                self.testSets[tsId].addMetadata("SQL_src-"+ccdName ,  sqlCache)
                
                
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

