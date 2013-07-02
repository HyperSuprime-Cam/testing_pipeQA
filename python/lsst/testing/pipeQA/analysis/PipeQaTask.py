import sys, re, os
import time
import copy
import datetime
import argparse
import traceback
import numpy

import lsst.pex.config         as pexConfig
import lsst.pipe.base          as pipeBase
import lsst.testing.pipeQA     as pipeQA

from lsst.pex.logging          import Trace
from .ZeropointFitQaTask       import ZeropointFitQaTask
from .EmptySectorQaTask        import EmptySectorQaTask
from .AstrometricErrorQaTask   import AstrometricErrorQaTask
from .PhotCompareQaTask        import PhotCompareQaTask
from .PsfShapeQaTask           import PsfShapeQaTask
from .CompletenessQaTask       import CompletenessQaTask
from .VignettingQaTask         import VignettingQaTask
from .VisitToVisitPhotQaTask   import VisitToVisitPhotQaTask
from .VisitToVisitAstromQaTask import VisitToVisitAstromQaTask

from .SummaryQaTask            import SummaryQaTask


class PipeQaConfig(pexConfig.Config):
    
    doZptFitQa      = pexConfig.Field(dtype = bool,
                                      doc = "Photometric Zeropoint: qaAnalysis.ZeropointFitQaTask",
                                      default = True)
    doEmptySectorQa = pexConfig.Field(dtype = bool,
                                      doc = "Empty Sectors: qaAnalysis.EmptySectorQaTask", default = True)
    doAstromQa      = pexConfig.Field(dtype = bool,
                                      doc = "Astrometric Error: qaAnalysis.AstrometricErrorQaTask",
                                      default = True)
    doPhotCompareQa = pexConfig.Field(dtype = bool,
                                      doc = "Photometric Error: qaAnalysis.PhotCompareQaTask", default = True)
    doPsfShapeQa    = pexConfig.Field(dtype = bool,
                                      doc = "Psf Shape: qaAnalysis.PsfShapeQaTask", default = True)
    doCompleteQa    = pexConfig.Field(dtype = bool,
                                      doc = "Photometric Depth: qaAnalysis.CompletenessQaTask",
                                      default = True)
    doVignettingQa  = pexConfig.Field(dtype = bool,
                                      doc = "Vignetting testing: qaAnalysis.VignettingQaTask",
                                      default = True)
    doVisitQa       = pexConfig.Field(dtype = bool,
                                      doc = "Visit to visit: qaAnalysis.VisitToVisitPhotQaTask and "+
                                      "qaAnalysis.VisitToVisitAstromQaTask", default = False)

    doSummaryQa     = pexConfig.Field(dtype = bool,
                                      doc = "Include a Summary Qa Page.",
                                      default = True)

    
    zptFitQa        = pexConfig.ConfigurableField(target = ZeropointFitQaTask,
                                                  doc = "Quality of zeropoint fit")
    emptySectorQa   = pexConfig.ConfigurableField(target = EmptySectorQaTask,
                                                  doc = "Look for missing matches")
    astromQa        = pexConfig.ConfigurableField(target = AstrometricErrorQaTask,
                                                  doc = "Quality of astrometric fit")
    photCompareQa   = pexConfig.ConfigurableField(target = PhotCompareQaTask,
                                                  doc = "Quality of photometry")
    psfShapeQa      = pexConfig.ConfigurableField(target = PsfShapeQaTask,
                                                  doc = "Shape of Psf")
    completeQa      = pexConfig.ConfigurableField(target = CompletenessQaTask,
                                                  doc = "Completeness of detection")
    vignettingQa    = pexConfig.ConfigurableField(target = VignettingQaTask,
                                                  doc = "Look for residual vignetting features")
    vvPhotQa        = pexConfig.ConfigurableField(target = VisitToVisitPhotQaTask,
                                                  doc = "Visit to visit photometry")
    vvAstromQa      = pexConfig.ConfigurableField(target = VisitToVisitAstromQaTask,
                                                  doc = "Visit to visit astrometry")

    summaryQa       = pexConfig.ConfigurableField(target = SummaryQaTask,
                                                  doc = "Include Summary QA information")

    shapeAlgorithm = pexConfig.ChoiceField(
        dtype = str,
        doc = "Shape Algorithm to load",
        default = "HSM_REGAUSS",
        allowed = {
            "HSM_REGAUSS": "Hirata-Seljac-Mandelbaum Regaussianization",
            "HSM_BJ" : "Hirata-Seljac-Mandelbaum Bernstein&Jarvis",
            "HSM_LINEAR" : "Hirata-Seljac-Mandelbaum Linear",
            "HSM_KSB" : "Hirata-Seljac-Mandelbaum KSB",
            "HSM_SHAPELET" : "Hirata-Seljac-Mandelbaum SHAPELET"
            }
        )


    
class PipeQaTask(pipeBase.Task):
    
    ConfigClass = PipeQaConfig
    _DefaultName = "pipeQa"
    
    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        
    def _makeArgumentParser(self):
        parser = argparse.ArgumentParser(usage=__doc__,
                                         formatter_class = argparse.RawDescriptionHelpFormatter)
        
        parser.add_argument("dataset", help="Dataset to use")
        parser.add_argument("-b", "--breakBy", default="visit",
                            help="Break the run by 'visit','raft', or 'ccd' (default=%(default)s)")
        parser.add_argument("-C", "--camera", default="lsstSim",
                            help="Specify a camera and override auto-detection (default=%(default)s)")
        parser.add_argument("-c", "--ccd", default=".*",
                            help="Specify ccd as regex (default=%(default)s)")
        parser.add_argument("-d", "--dataSource", default="db",
                            help="Specify the source of data to load",
                            choices=('butler', 'db'))
        parser.add_argument("-e", "--exceptExit", default=False, action='store_true',
                            help="Don't capture exceptions, fail and exit (default=%(default)s)")
        parser.add_argument("-f", "--forkFigure", default=False, action='store_true',
                            help="Make figures in separate process (default=%(default)s)")
        parser.add_argument("-F", "--useForced", default=False, action='store_true',
                            help="Use forced photometry (default=%(default)s)")        
        parser.add_argument("-g", "--group", default=None,
                            help="Specify sub-group of visits 'groupSize:whichGroup' (default=%(default)s)")
        parser.add_argument("-k", "--keep", default=False, action="store_true",
                            help="Keep existing outputs (default=%(default)s)")
        parser.add_argument("-r", "--raft", default=".*",
                            help="Specify raft as regex (default=%(default)s)")
        parser.add_argument("-R", "--rerun", default=None,
                            help="Rerun to analyse - only valid for hsc/suprimecam (default=%(default)s)")
        parser.add_argument("-s", "--snap", default=".*",
                            help="Specify snap as regex (default=%(default)s)")
        parser.add_argument("-S", "--summaryProcessing", default="delay",
                            help="Processing summary figures ['delay', 'none', or 'summOnly'] (default=%(default)s)")
        parser.add_argument("-t", "--test", default=".*",
                            help="Regex specifying which QaAnalysis to run (default=%(default)s)")
        parser.add_argument("-T", "--coaddTable", default="goodSeeing",
                            help="Specify coadd table to use (default=%(default)s)")        
        parser.add_argument("-V", "--verbosity", default=1,
                            help="Trace level for lsst.testing.pipeQA")
        parser.add_argument("-v", "--visit", default=".*",
                            help="Specify visit as regex OR color separated list. (default=%(default)s)")
        parser.add_argument("-z", "--lazyPlot", default='sensor',
                            help="Figures to be generated dynamically online "+
                            "[options: none, sensor, all] (default=%(default)s)")

        
        # visit-to-visit
        parser.add_argument("--doVisitQa", default=False, action='store_true',
                            help="Do visit-to-visit pipeQA, overriding the policy default")
        parser.add_argument("--matchDataset", default=None,
                            help="Specify another dataset to compare analysis to")
        parser.add_argument("--matchVisits", default=None, action="append",
                            help="Visits within this dataset to compare analysis to; default is same visit")
        
        
        parser.add_argument("--noWwwCache", default=False, action="store_true",
                            help="Disable caching of pass/fail (needed to run in parallel) (default=%(default)s)")

        # and add in ability to override config
        config = self.ConfigClass()
        config.photCompareQa.compareTypes = ("psf ap", "psf inst", "psf cat", "ap cat")        
        parser.set_defaults(config = config)
        parser.add_argument("--config", nargs="*", action=pipeBase.argumentParser.ConfigValueAction,
                            help="config override(s), e.g. -c foo=newfoo bar.baz=3", metavar="NAME=VALUE")

        return parser


    
    @staticmethod
    def _getMemUsageThisPid(size = "rss"):
        """Generalization; memory sizes: rss, rsz, vsz."""
        return int(os.popen('ps -p %d -o %s | tail -1' % (os.getpid(), size)).read())
    


    
    def makeSubtask(self, taskStr, **kwargs):
        self.log.log(self.log.INFO,
                     "Initializing task %s: %s" % (taskStr, " ".join(map(str, kwargs.values()))))

        configurableField = getattr(self.config, taskStr, None)
        if configurableField is None:
            raise KeyError("%s's config does not have field %r" % (self.getFullName, taskStr)) 
        return configurableField.apply(name=taskStr, parentTask=self, **kwargs)


    
    def runSubtask(self, subtask, data, thisDataId, visit, test, testset, exceptExit):
    
        subtaskName = subtask.__name__

        if exceptExit:

            if subtaskName == 'free':
                subtask()
            else:
                subtask(data, thisDataId)
        else:
            
            if thisDataId.has_key('raft'):
                label = "v%s_r%s_s%s_%s_%s" % (visit, thisDataId['raft'],
                                               thisDataId[data.ccdConvention], test, subtaskName)
            else:
                label = "v%s_s%s_%s_%s" % (visit, thisDataId[data.ccdConvention], test, subtaskName)


            failed = False
            s = ""
            try:
                if subtaskName == 'free':
                    subtask()
                else:
                    subtask(data, thisDataId)
            except Exception, e:
                failed = True
                exc_type, exc_value, exc_traceback = sys.exc_info()
                s = traceback.format_exception(exc_type, exc_value,
                                               exc_traceback)

                self.log.log(self.log.WARN, "Warning: Exception in QA processing of %s: %s" % (label, str(e)))

            if failed:
                testset.addTest(label, 1, [0, 0], "QA exception thrown (%s)" % (str(thisDataId)),
                                backtrace="".join(s))


                
            
    @pipeBase.timeMethod
    def parseAndRun(self, args):
        self.log.log(self.log.INFO, "PipeQA Start")

        argumentParser = self._makeArgumentParser()
        parsedCmd = argumentParser.parse_args(args)
        self.config = parsedCmd.config # include any overrides
        self.config.validate()
        self.config.freeze()

        # split by :
        visits = parsedCmd.visit
        if re.search(":", parsedCmd.visit):
            visits = parsedCmd.visit.split(":")
        dataIdInput = {
            'visit' : visits,
            'sensor': parsedCmd.ccd,
            'raft'  : parsedCmd.raft,
            'snap'  : parsedCmd.snap,
            }
        rerun        = parsedCmd.rerun
        exceptExit   = parsedCmd.exceptExit
        testRegex    = parsedCmd.test
        camera       = parsedCmd.camera
        retrievalType = parsedCmd.dataSource
        keep         = parsedCmd.keep
        breakBy      = parsedCmd.breakBy
        groupInfo    = parsedCmd.group
        summaryProcessing = parsedCmd.summaryProcessing
        forkFigure   = parsedCmd.forkFigure
        wwwCache     = not parsedCmd.noWwwCache
        useForced    = parsedCmd.useForced
        coaddTable   = parsedCmd.coaddTable
        lazyPlot     = parsedCmd.lazyPlot
        verbosity    = parsedCmd.verbosity

        summOpts = ["delay", "none", "summOnly"]
        if not summaryProcessing in summOpts:
            raise ValueError("summaryProcessing must be: "+", ".join(summOpts))
            
        # optional visitQA info
        matchDset    = parsedCmd.matchDataset
        matchVisits  = parsedCmd.matchVisits

        # finally, the dataset to run!
        dataset      = parsedCmd.dataset

        if not re.search("^(visit|raft|ccd)$", breakBy):
            self.log.log(self.log.FATAL, "breakBy (-b) must be 'visit', 'raft', or 'ccd'")
            sys.exit()
    
        if breakBy in ['raft', 'ccd'] and not keep:
            self.log.log(self.log.WARN, "You've specified breakBy=%s, which requires 'keep' (-k). "+
                         "I'll set it for you.")
            keep = True
        
        # Is this deprecated?
        Trace.setVerbosity('lsst.testing.pipeQA', int(verbosity))
        
        visitList = []
        if isinstance(dataIdInput['visit'], list):
            visitList = dataIdInput['visit']
            dataIdInput['visit'] = ".*"
            
        #if exceptExit:
        #    numpy.seterr(all="raise")
    
        data = pipeQA.makeQaData(dataset, rerun=rerun, camera=camera,
                                 shapeAlg = self.config.shapeAlgorithm,
                                 retrievalType=retrievalType,
                                 useForced=useForced, coaddTable=coaddTable)
    
        if data.cameraInfo.name == 'lsstSim' and  dataIdInput.has_key('ccd'):
            dataIdInput['sensor'] = dataIdInput['ccd']
            del dataIdInput['ccd']

        # convert this input format visit,raft,ccd to the names used by the instrument
        dataIdOrig = copy.copy(dataIdInput) # good to have for debugging
        dataIdInput = data.cameraInfo.dataIdStandardToCamera(dataIdInput)    
            
        # take what we need for this camera, ignore the rest
        dataId = {}
        for name in data.dataIdNames:
            if dataIdInput.has_key(name):
                dataId[name] = dataIdInput[name]
    
        # if they requested a key that doesn't exist for this camera ... throw
        for k, v in dataIdInput.items():
            if (not k in data.dataIdNames) and (v != '.*'):
                raise Exception("Key "+k+" not available for this dataset (camera="+data.cameraInfo.name+")")
    
            
        taskList = []
        # Simple ones
        for doTask, taskStr in ( (self.config.doZptFitQa,      "zptFitQa"),
                                 (self.config.doEmptySectorQa, "emptySectorQa"),
                                 (self.config.doAstromQa,      "astromQa"),
                                 (self.config.doPsfShapeQa,    "psfShapeQa"),
                                 (self.config.doCompleteQa,    "completeQa"),
                                 (self.config.doVignettingQa,  "vignettingQa"),
                                 (self.config.doSummaryQa,     "summaryQa") ):
            
            if doTask and (data.cameraInfo.name in eval("self.config.%s.cameras" % (taskStr))):
                stask = self.makeSubtask(taskStr, useCache=keep, wwwCache=wwwCache,
                                         summaryProcessing=summaryProcessing, lazyPlot=lazyPlot)
                taskList.append(stask)

                
        # Multiple permutations
        if self.config.doPhotCompareQa and (data.cameraInfo.name in self.config.photCompareQa.cameras):
            for types in self.config.photCompareQa.compareTypes:
                mag1, mag2 = types.split()
                starGxyToggle = types in self.config.photCompareQa.starGalaxyToggle
                stask = self.makeSubtask("photCompareQa", magType1=mag1, magType2=mag2,
                                         starGalaxyToggle=starGxyToggle, useCache=keep, wwwCache=wwwCache,
                                         summaryProcessing=summaryProcessing, lazyPlot=lazyPlot)
                taskList.append(stask)


        # Additional dependencies needed
        if self.config.doVisitQa:
            if matchDset == None and matchVisits == None:
                # we can't do it!
                self.log.log(self.log.FATAL, "Unable to run visit to visit Qa; "+
                             "please request a comparison visit or database")
                sys.exit(1)

            elif matchDset == None:
                matchDset = dataset

            elif matchVisits == None:
                matchVisits = []

            if data.cameraInfo.name in self.config.vvPhotQa.cameras:
                for mType in self.config.vvPhotQa.magTypes:
                    stask = self.makeSubtask("vvPhotQa", matchDset=matchDset, matchVisits=matchVisits,
                                             mType=mType, useCache=keep, wwwCache=wwwCache,
                                             summaryProcessing=summaryProcessing, lazyPlot=lazyPlot)
                    taskList.append(stask)

            if data.cameraInfo.name in self.config.vvAstromQa.cameras:
                stask = self.makeSubtask("vvAstromQa", matchDset = matchDset, matchVisits = matchVisits, 
                                         useCache=keep, wwwCache=wwwCache, summaryProcessing=summaryProcessing,
                                         lazyPlot=lazyPlot)
                taskList.append(stask)
                
                
        # Split by visit, and handle specific requests
        visitsTmp = data.getVisits(dataId)
        visits = []
        if len(visitList) > 0:
            for v in visitsTmp:
                if v in visitList:
                    visits.append(v)
        else:
            visits = visitsTmp
    
        groupTag = ""
        if not groupInfo is None:
            groupSize, whichGroup = map(int, groupInfo.split(":"))
            lo, hi = whichGroup*groupSize, (whichGroup+1)*groupSize
    
            nvisit = len(visits)
            if lo >= nvisit:
                self.log.log(self.log.FATAL,
                             "Can't run visits %d to %d as there are only %d visits" % (lo, hi, nvisit))
                sys.exit()
            if hi > nvisit:
                hi = nvisit
            visits = visits[lo:hi]
    
            self.log.log(self.log.INFO,
                         "Total of %d visits grouped by %d.  Running group %d with visits %d - %d:\n%s\n" % \
                             (nvisit, groupSize, whichGroup, lo, hi-1, "\n".join(visits)))
            groupTag = "%02d-%02d" % (groupSize, whichGroup)
    
        

        for visit in visits:
            visit_t0 = time.time()
    
            dataIdVisit = copy.copy(dataId)
            dataIdVisit['visit'] = visit
    
            # now break up the run into eg. rafts or ccds
            #  ... if we only run one raft or ccd at a time, we use less memory
            brokenDownDataIdList = data.breakDataId(dataIdVisit, breakBy)

            #if this is a summary only run, shortcircuit to the last ccd
            if summaryProcessing in ['summOnly']:
                brokenDownDataIdList = [brokenDownDataIdList[-1]]

            for thisDataId in brokenDownDataIdList:

                haveIt = data.verify(thisDataId)
                if not haveIt:
                    self.log.log(self.log.WARN, "Missing dataId="+str(thisDataId))
                    continue
                
                raftName, ccdName = data.cameraInfo.getRaftAndSensorNames(thisDataId)
                ccdName = data.cameraInfo.getDetectorName(raftName, ccdName)
                testset = pipeQA.TestSet(group="", label="QA-failures"+groupTag, wwwCache=wwwCache, sqliteSuffix=ccdName)

                for task in taskList:
    
                    test_t0 = time.time()
                    test = str(task)
                    if not re.search(testRegex, test):
                        continue
    
                    date = datetime.datetime.now().strftime("%a %Y-%m-%d %H:%M:%S")
                    visitLog = str(visit)
                    if breakBy.lower() == 'ccd':
                        visitLog += " ccd:" + str(thisDataId['ccd'])
                    self.log.log(self.log.INFO, "Running " + test + "  visit:" + visitLog + "  ("+date+")")
                    sys.stdout.flush() # clear the buffer before the fork
    
    
                    # try the test() method
                    t0 = time.time()
                    if summaryProcessing in ['delay', 'none']:
                        self.runSubtask(task.test, data, thisDataId, visit, test, testset, exceptExit)


                    t0 = time.time()
                    if forkFigure:
                        pid = os.fork()
                        if pid == 0:
                            # try the plot() method
                            self.runSubtask(task.plot, data, thisDataId, visit, test, testset, exceptExit)
                            os._exit(os.EX_OK) # sys.exit(0)
                        else:
                            os.waitpid(pid, 0)
                            
                    else:
                        # try the plot() method
                        self.runSubtask(task.plot, data, thisDataId, visit, test, testset, exceptExit)

                    
                    # try the free() method
                    # test() method only ran for 'delay' and 'none'.  Only then is there stuff to free
                    if summaryProcessing in ['delay', 'none']:
                        self.runSubtask(task.free, data, thisDataId, visit, test, testset, exceptExit)
    

                        
                # we're now done this dataId ... can clear the cache            
                data.clearCache()
                raftName = ""
                if thisDataId.has_key('raft'):
                    raftName = str(thisDataId['raft'])+"-"
                ccdName = ""
                if thisDataId.has_key("ccd"):
                    ccdName = str(thisDataId[data.ccdConvention])

                
        if summaryProcessing in ['summOnly']:
            ts = pipeQA.TestSet(group="", label="QA-failures"+groupTag, wwwCache=wwwCache, sqliteSuffix="")
            ts.accrete()
            ts.updateCounts()


        self.log.log(self.log.INFO, "PipeQA End")
        return pipeBase.Struct()


SETUP = """
setup meas_algorithms 4.9.0.1+1
setup -k meas_astrom 4.9.0.0+1
setup -k -r ~/LSST/DMS/testing_displayQA/
setup -k -r ~/LSST/DMS/testing_pipeQA/
setup -k -r ~/LSST/DMS/pipe_base/
setup -k pex_config 5.0.0.1+1
cd ~/LSST/DMS/testing_pipeQA/
setenv WWW_RERUN ticket2040
python bin/pipeQa.py -v 885449191 -r "2,2" -c "1,1" -e abecker_PT1_2_u_becker_2012_0209_181253 --config doEmptySectorQa=False
"""
