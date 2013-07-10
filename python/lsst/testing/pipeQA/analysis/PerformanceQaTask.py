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



def getMemUsageThisPid(size="rss"):
    """Generalization; memory sizes: rss, rsz, vsz."""
    #return 1.0
    return int(os.popen('ps -p %d -o %s | tail -1' % (os.getpid(), size)).read())



class PerformanceQaConfig(pexConfig.Config):
    cameras = pexConfig.ListField(dtype = str,
                                  doc = "Cameras to run PerformanceQaTask",
                                  default = ("lsstSim", "cfht", "suprimecam", "hscSim", "sdss", "coadd"))

    
class PerformanceQaTask(QaAnalysisTask):
    ConfigClass = PerformanceQaConfig
    _DefaultName = "performanceQa"

    def __init__(self, **kwargs):
        QaAnalysisTask.__init__(self, **kwargs)

        self.node   = platform.node()
        self.dist   = platform.dist() # a tuple e.g., ('redhat', '6.2', 'Santiago')
        
        fp_meminfo  = open('/proc/meminfo')
        meminfoList = fp_meminfo.readlines()
        fp_meminfo.close()

        self.meminfo = {}
        for line in meminfoList:
            fields = line.split()
            self.meminfo[fields[0]] = int(fields[1])

        self.memtotal = self.meminfo['MemTotal:']/1024.0  # convert to MB
        self.limits = [0.0, 0.125*self.memtotal]
        
        self.description = """
        The page summarizes various performance parameters associated with the pipeQA run.
        """

    def free(self):

        # please free all large data structures here.
        pass
        
    def test(self, data, dataId):
        
        self.detector         = data.getDetectorBySensor(dataId)
        self.filter           = data.getFilterBySensor(dataId)

        # create containers for data we're interested in
        self.mem   = raftCcdData.RaftCcdData(self.detector)
        self.testRuntime = raftCcdData.RaftCcdData(self.detector)
        self.plotRuntime = raftCcdData.RaftCcdData(self.detector)
        
        # create a testset
        testSet = self.getTestSet(data, dataId)

        # this normally gets set in the plot as that's where the caching happens,
        # here we're stashing the nDetection and nCcd values, so we need to set it early.
        testSet.setUseCache(self.useCache)
        testSet.addMetadata({"Description": self.description})

        
        performBase = "performShelf"
        mb = testSet.unshelve(performBase)
        
        for raft, ccd in self.mem.raftCcdKeys():
            
            # add tests for acceptible numpy of empty sectors
            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            
            label = "memory usage"
            mem = getMemUsageThisPid()/1024 # MB
            comment = "memory usage in MB (%.2f%% of %.0fMB sys.)" % (100.0*mem/self.memtotal, self.memtotal)
            test = testCode.Test(label, mem, self.limits, comment, areaLabel=areaLabel)
            testSet.addTest(test)
            mb[ccd] = mem
            self.mem.set(raft, ccd, mem)

            testRuntime = data.getPerformance(dataId, 'total', 'test-runtime')
            if testRuntime is None:
                testRuntime = 0.0
            self.testRuntime.set(raft, ccd, testRuntime)

            plotRuntime = data.getPerformance(dataId, 'total', 'plot-runtime')
            # if running with plot() forked, plot times are not defined.
            if plotRuntime is None:
                plotRuntime = 0.0
            self.plotRuntime.set(raft, ccd, plotRuntime)

            tTest = testCode.Test("test-runtime", testRuntime, [0.0, 3600],
                                  "Runtime for test()[s]", areaLabel=areaLabel)
            testSet.addTest(tTest)
                
            pTest = testCode.Test("plot-runtime", plotRuntime, [0.0, 3600],
                                  "Runtime for plot()[s] (%d plots)" % (qaFig.QaFigure.count),
                                  areaLabel=areaLabel)
            testSet.addTest(pTest)

            info = self.node + " " + " ".join(self.dist)
            infoTest = testCode.Test("generalInfo", 0.5, [0.0, 1.0], info, areaLabel=areaLabel)
            testSet.addTest(infoTest)
            
        testSet.shelve(performBase, mb)


    def plot(self, data, dataId, showUndefined=False):

        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)
        isFinalDataId = False
        if len(data.brokenDataIdList) == 0 or data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True

        #################################
        # memory
        ################################


        # make fpa figures - for all detections, and for matched detections
        memBase = "mem"

        memFig    = qaFig.FpaQaFigure(data.cameraInfo, data=None, map=None)

        if self.summaryProcessing != self.summOpt['summOnly']:
            for raft, ccdDict in memFig.data.items():
                for ccd, value in ccdDict.items():
                    # set values for data[raft][ccd] (color coding)
                    # set values for map[raft][ccd]  (tooltip text)
                    if not self.mem.get(raft, ccd) is None:
                        mem = self.mem.get(raft, ccd)
                        memFig.data[raft][ccd] = mem
                        memFig.map[raft][ccd] = "%.1fMB" % (mem)

                        label = data.cameraInfo.getDetectorName(raft, ccd)
                        testSet.pickle(memBase+label, [memFig.data, memFig.map])


        # make the figures and add them to the testSet
        # sample colormaps at: http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps

        if (self.summaryProcessing in [self.summOpt['summOnly'], self.summOpt['delay']]) and isFinalDataId:

            for raft, ccdDict in memFig.data.items():
                for ccd, value in ccdDict.items():
                    label = data.cameraInfo.getDetectorName(raft, ccd)
                    memDataTmp, memMapTmp       = testSet.unpickle(memBase+label, [None, None])
                    memFig.mergeValues(memDataTmp, memMapTmp)

            self.log.log(self.log.INFO, "plotting FPAs")
            memFig.makeFigure(showUndefined=showUndefined, cmap="gist_heat_r",
                                  vlimits=self.limits, 
                                  title="Memory Usage [MB]",
                                  failLimits=self.limits)
            testSet.addFigure(memFig, memBase+".png",
                              "memory usage in MB", navMap=True)
            
        del memFig



        #################################
        # runtime
        ################################

        # make fpa figures - for all detections, and for matched detections
        runtimeBase = "runtime"

        runtimeFig    = qaFig.FpaQaFigure(data.cameraInfo, data=None, map=None)

        #print "min/max", mintime, maxtime
        if self.summaryProcessing != self.summOpt['summOnly']:
            for raft, ccdDict in runtimeFig.data.items():
                for ccd, value in ccdDict.items():
                    if not self.testRuntime.get(raft, ccd) is None:
                        testRuntime = self.testRuntime.get(raft, ccd)
                        plotRuntime = self.plotRuntime.get(raft, ccd)
                        runtimeFig.data[raft][ccd] = testRuntime + plotRuntime
                        runtimeFig.map[raft][ccd] = "t=%.1fsec+p=%.1fsec" % (testRuntime, plotRuntime)

                        label = data.cameraInfo.getDetectorName(raft, ccd)
                        testSet.pickle(runtimeBase+label, [runtimeFig.data, runtimeFig.map])



        # make the figures and add them to the testSet
        # sample colormaps at: http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps
        if (self.summaryProcessing in [self.summOpt['summOnly'], self.summOpt['delay']]) and isFinalDataId:
            runArray = runtimeFig.getArray()
            if len(runArray) == 0:
                runArray = numpy.zeros(1)
            mintime, maxtime = runArray.min()-1.0, runArray.max()+1.0

            for raft, ccdDict in runtimeFig.data.items():
                for ccd, value in ccdDict.items():
                    label = data.cameraInfo.getDetectorName(raft, ccd)
                    runtimeDataTmp, runtimeMapTmp       = testSet.unpickle(runtimeBase+label, [None, None])
                    runtimeFig.mergeValues(runtimeDataTmp, runtimeMapTmp)
            
            runtimeFig.makeFigure(showUndefined=showUndefined, cmap="gist_heat_r",
                                  vlimits=[mintime, maxtime], 
                                  title="Runtime [Sec]",
                                  failLimits=[mintime, maxtime])
            testSet.addFigure(runtimeFig, runtimeBase+".png",
                              "Runtime", navMap=True)
                
        del runtimeFig


        if (self.summaryProcessing in [self.summOpt['summOnly'], self.summOpt['delay']]) and isFinalDataId:
            self.combineOutputs(data, dataId)
