import sys, os, re
import numpy
import platform


import lsst.meas.algorithms         as measAlg
import lsst.testing.pipeQA.figures  as qaFig
import lsst.afw.math                as afwMath
import lsst.testing.pipeQA.TestCode as testCode

import RaftCcdData                  as raftCcdData
import QaAnalysisUtils              as qaAnaUtil

from   .QaAnalysisTask              import QaAnalysisTask
import lsst.pex.config              as pexConfig



class SummaryQaConfig(pexConfig.Config):
    cameras = pexConfig.ListField(dtype = str,
                                  doc = "Cameras to run SummaryQaTask",
                                  default = ("suprimecam", "hscSim"))

    
class SummaryQaTask(QaAnalysisTask):
    ConfigClass = SummaryQaConfig
    _DefaultName = "summaryQa"

    def __init__(self, **kwargs):
        QaAnalysisTask.__init__(self, **kwargs)

        self.overscan_limits = [0.0, 2.0]
        
        self.description = """
        The page summarizes observing parameters for a visit.
        """

    def free(self):

        # please free all large data structures here.
        pass

        
    def test(self, data, dataId):
        
        self.detector         = data.getDetectorBySensor(dataId)
        self.filter           = data.getFilterBySensor(dataId)
        self.summary          = data.getSummaryDataBySensor(dataId)

        # create containers for data we're interested in
        self.overscan = raftCcdData.RaftCcdData(self.detector)
        
        # create a testset
        testSet = self.getTestSet(data, dataId)

        # this normally gets set in the plot as that's where the caching happens,
        # here we're stashing the nDetection and nCcd values, so we need to set it early.
        testSet.setUseCache(self.useCache)
        testSet.addMetadata({"Description": self.description})

        
        summaryBase = "summaryShelf"
        summDat = testSet.unshelve(summaryBase)

            
        for raft, ccd in self.overscan.raftCcdKeys():
            
            # add tests for acceptible numpy of empty sectors
            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)

            oscan = 1.0
            label = "overscan"
            comment = "value of overscan"
            test = testCode.Test(label, oscan, self.overscan_limits, comment, areaLabel=areaLabel)
            testSet.addTest(test)
            summDat[ccd] = oscan
            self.overscan.set(raft, ccd, oscan)

        testSet.shelve(summaryBase, summDat)


    def plot(self, data, dataId, showUndefined=False):

        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)
        isFinalDataId = False
        if len(data.brokenDataIdList) == 0 or data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True

        #################################
        # memory
        ################################

        if True:
            # make fpa figures - for all detections, and for matched detections
            summBase = "summ"

            summData, summMap       = testSet.unpickle(summBase, [None, None])
            summFig    = qaFig.FpaQaFigure(data.cameraInfo, data=summData, map=summMap)

            for raft, ccdDict in summFig.data.items():
                for ccd, value in ccdDict.items():

                    # set values for data[raft][ccd] (color coding)
                    # set values for map[raft][ccd]  (tooltip text)
                    if not self.overscan.get(raft, ccd) is None:
                        summ = self.overscan.get(raft, ccd)
                        summFig.data[raft][ccd] = summ
                        summFig.map[raft][ccd] = "%.1f" % (summ)

            testSet.pickle(summBase, [summFig.data, summFig.map])


            # make the figures and add them to the testSet
            # sample colormaps at: http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps
            if not self.delaySummary or isFinalDataId:
                self.log.log(self.log.INFO, "plotting FPAs")
                summFig.makeFigure(showUndefined=showUndefined, cmap="gist_heat_r",
                                      vlimits=self.overscan_limits, 
                                      title="Overscan",
                                      failLimits=self.overscan_limits)
                testSet.addFigure(summFig, summBase+".png",
                                  "summary", navMap=True)
                del summFig

            else:
                del summFig



