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

import QaPlotUtils                  as qaPlotUtil



class SummaryQaConfig(pexConfig.Config):
    cameras = pexConfig.ListField(dtype = str,
                                  doc = "Cameras to run SummaryQaTask",
                                  default = ("suprimecam", "hscSim"))

    
class SummaryQaTask(QaAnalysisTask):
    ConfigClass = SummaryQaConfig
    _DefaultName = "summaryQa"

    def __init__(self, **kwargs):
        QaAnalysisTask.__init__(self, **kwargs)

        self.limits = {
            'oslevel1'    : [100.0, 250.0],
            'oslevel2'    : [100.0, 250.0],
            'oslevel3'    : [100.0, 250.0],
            'oslevel4'    : [100.0, 250.0],
            'ossigma1'    : [1.0, 2.5],
            'ossigma2'    : [1.0, 2.5],
            'ossigma3'    : [1.0, 2.5],
            'ossigma4'    : [1.0, 2.5],
            'gain1'       : [2.0, 4.0],
            'gain2'       : [2.0, 4.0],
            'gain3'       : [2.0, 4.0],
            'gain4'       : [2.0, 4.0],
            'seeing'      : [0.0, 10.0],
            'ellipt'      : [0.0, 0.12],
            'ell_pa'      : [0.0, 360], 
            'flatness_rms': [0.0, 1.0],
            'flatness_pp' : [0.0, 1.0],
            'skylevel'    : [1.0, 1.0e5],
            'sigma_sky'   : [1.0, 1.0e3],
            'ccdtemp'     : [100.0,200.0],
            }
        self.fromCalexp = ['ccdtemp']
        



        self.aggregate_limits = [0.0, len(self.limits.keys())]
        
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
        self.calexp           = data.getCalexpBySensor(dataId)

        # create containers for data we're interested in
        self.overscan     = raftCcdData.RaftCcdData(self.detector)
        self.seeing       = raftCcdData.RaftCcdData(self.detector)
        self.ellipt       = raftCcdData.RaftCcdData(self.detector)
        self.ell_pa       = raftCcdData.RaftCcdData(self.detector)
        self.flatness_rms = raftCcdData.RaftCcdData(self.detector)
        self.flatness_pp  = raftCcdData.RaftCcdData(self.detector)
        self.skylevel     = raftCcdData.RaftCcdData(self.detector)
        self.sigma_sky    = raftCcdData.RaftCcdData(self.detector)
        self.ccdtemp      = raftCcdData.RaftCcdData(self.detector)
        
        self.aggregate    = raftCcdData.RaftCcdData(self.detector)
        
        # create a testset
        testSet = self.getTestSet(data, dataId)

        # this normally gets set in the plot as that's where the caching happens,
        # here we're stashing the nDetection and nCcd values, so we need to set it early.
        testSet.setUseCache(self.useCache)
        testSet.addMetadata({"Description": self.description})

        
        summaryBase = "summaryShelf"
        summDat = testSet.unshelve(summaryBase)

        key = None
        for key in self.summary.keys():
            raft = self.detector[key].getParent().getId().getName()
            ccd  = self.detector[key].getId().getName()

            aggregate = 0            
            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)

            for s in self.limits.keys():
                if s.upper() in self.summary[key]:
                    value = self.summary[key][s.upper()]
                elif s in self.calexp[key]:
                    value = self.calexp[key][s]
                else:
                    value = None
                label = s.lower()
                comment = s
                test = testCode.Test(label, value, self.limits[s], comment, areaLabel=areaLabel)
                testSet.addTest(test)
                #summDat[ccd] = skylevel

                # count failures to color code the ccd
                if not test.evaluate():
                    aggregate += 1 

            self.aggregate.set(raft, ccd, aggregate)

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
                    if not self.aggregate.get(raft, ccd) is None:
                        summ = self.aggregate.get(raft, ccd)
                        summFig.data[raft][ccd] = summ
                        summFig.map[raft][ccd] = "%.1f" % (summ)

            testSet.pickle(summBase, [summFig.data, summFig.map])


            # make the figures and add them to the testSet
            # sample colormaps at: http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps
            if not self.delaySummary or isFinalDataId:
                self.log.log(self.log.INFO, "plotting FPAs")
                summFig.makeFigure(showUndefined=showUndefined, cmap="gist_heat_r",
                                      vlimits=self.aggregate_limits, 
                                      title="No. of Parameters outside Specified Range",
                                      failLimits=self.aggregate_limits)
                testSet.addFigure(summFig, summBase+".png",
                                  "summary", navMap=True)
                del summFig

            else:
                del summFig



        if False:




            cacheLabel = "metadata_sum" #cache
            shelfData = {}

            # make any individual (ie. per sensor) plots
            for raft, ccd in self.sigma_sky.raftCcdKeys():

                x = 2.0*numpy.pi*numpy.arange(100)/100
                y = numpy.sin(x) + numpy.random.uniform(0.0, 10.0)

                dataDict = {'x' : x, 'y': y,
                            'summary' : False,
                            }

                self.log.log(self.log.INFO, "plotting %s" % (ccd))
                import SummaryQaPlot as plotModule
                label = data.cameraInfo.getDetectorName(raft, ccd)
                caption = "Dummy metadata summary " + label
                pngFile = cacheLabel+".png"

                if self.lazyPlot.lower() in ['sensor', 'all']:
                    testSet.addLazyFigure(dataDict, pngFile, caption,
                                          plotModule, areaLabel=label, plotargs="")
                else:
                    testSet.cacheLazyData(dataDict, pngFile, areaLabel=label)
                    fig = plotModule.plot(dataDict)
                    testSet.addFigure(fig, pngFile, caption, areaLabel=label)
                    del fig


            if not self.delaySummary or isFinalDataId:
                self.log.log(self.log.INFO, "plotting Summary figure")


                import SummaryQaPlot as plotModule
                label = 'all'
                caption = "Dummy metadata summary "+label
                pngFile = cacheLabel + ".png"

                if self.lazyPlot in ['all']:
                    testSet.addLazyFigure({}, pngFile, caption,
                                          plotModule, areaLabel=label, plotargs="")
                else:
                    dataDict, isSummary = qaPlotUtil.unshelveGlob(cacheLabel+"-all.png", testSet=testSet)
                    dataDict['summary'] = True
                    fig = plotModule.plot(dataDict)                
                    testSet.addFigure(fig, pngFile, caption, areaLabel=label)
                    del fig


            
