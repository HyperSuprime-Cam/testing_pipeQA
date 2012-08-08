#!/usr/bin/env python
#
# Original filename: bin/pipeQa.py
#
# Author: 
# Email: 
# Date: Mon 2011-04-11 13:11:01
# 
# Summary: 
# 
"""
%prog [options] dataset
  dataset = Directory in TESTBED_PATH, full path to data, or database name
"""

import sys
import re
import optparse
import os
import datetime
import traceback
import copy
import time

import lsst.testing.pipeQA as pipeQA
import lsst.testing.pipeQA.analysis     as qaAnalysis
import lsst.pex.policy as pexPolicy
from lsst.pex.logging import Trace

import numpy

def getMemUsageThisPid(size="rss"):
    """Generalization; memory sizes: rss, rsz, vsz."""
    #return 1.0
    return int(os.popen('ps -p %d -o %s | tail -1' % (os.getpid(), size)).read())



def tryThis(func, data, thisDataId, visit, test, testset, exceptExit):

    funcName = func.__name__
    
    if exceptExit:
        if funcName == 'free':
            func()
        else:
            func(data, thisDataId)

            
    else:
    
        if thisDataId.has_key('raft'):
            label = "v%s_r%s_s%s_%s_%s" % (visit, thisDataId['raft'], thisDataId[data.ccdConvention],
                                           test, funcName)
        else:
            label = "v%s_s%s_%s_%s" % (visit, thisDataId[data.ccdConvention], test, funcName)


        failed = False
        s = ""
        try:
            if funcName == 'free':
                func()
            else:
                func(data, thisDataId)
        except Exception, e:
            failed = True
            exc_type, exc_value, exc_traceback = sys.exc_info()
            s = traceback.format_exception(exc_type, exc_value,
                                           exc_traceback)

            print "Warning: Exception in QA processing of %s: %s" % (label, str(e))
            #print "       :", "".join(s)

        if failed:
            testset.addTest(label, 1, [0, 0], "QA exception thrown (%s)" % (str(thisDataId)),
                            backtrace="".join(s))



#############################################################
#
# Main body of code
#
#############################################################

def main(dataset, dataIdInput, rerun=None, doVisitQa=False, matchDset=None, matchVisits=None, testRegex=".*", camera=None,
         exceptExit=False, keep=False, wwwCache=True, breakBy='visit',
         groupInfo=None, delaySummary=False, forkFigure=False):

    visitList = []
    if isinstance(dataIdInput['visit'], list):
        visitList = dataIdInput['visit']
        dataIdInput['visit'] = ".*"
        
    if exceptExit:
        numpy.seterr(all="raise")

    # Policy for which plots to make and value of quality metrics
    policyDictName = "PipeQa.paf"
    policyFile = pexPolicy.DefaultPolicyFile("testing_pipeQA", policyDictName, "policy")
    policy  = pexPolicy.Policy.createPolicy(policyFile, policyFile.getRepositoryPath(), True)


    data = pipeQA.makeQaData(dataset, rerun=rerun, camera=camera,
                             shapeAlg=policy.get('shapeAlgorithm'))

    
    # convert this input format visit,raft,ccd to the names used by the instrument
    dataIdOrig = copy.copy(dataIdInput) # good to have for debugging
    dataIdInput = data.cameraInfo.dataIdStandardToCamera(dataIdInput)    
    
    # take what we need for this camera, ignore the rest
    dataId = {}
    for name in data.dataIdNames:
        if dataIdInput.has_key(name):
            dataId[name] = dataIdInput[name]

    data.reduceAvailableDataTupleList(dataId)
    
    # if they requested a key that doesn't exist for this camera ... throw
    for k, v in dataIdInput.items():
        if (not k in data.dataIdNames) and (v != '.*'):
            raise Exception("Key "+k+" not available for this dataset (camera="+data.cameraInfo.name+")")

    
    analysisList = []
    if data.cameraInfo.name in policy.getStringArray("doZptQa"):
        zptMin = policy.get("zptQaMetricMin")
        zptMax = policy.get("zptQaMetricMax")
        analysisList.append(qaAnalysis.ZeropointQaAnalysis(zptMin, zptMax, useCache=keep, wwwCache=wwwCache,
                                                           delaySummary=delaySummary))
    if data.cameraInfo.name in policy.getStringArray("doZptFitQa"):
        offsetMin = policy.get("zptFitQaOffsetMin")
        offsetMax = policy.get("zptFitQaOffsetMax")
        analysisList.append(qaAnalysis.ZeropointFitQa(offsetMin, offsetMax, useCache=keep, wwwCache=wwwCache,
                                                      delaySummary=delaySummary))
    if data.cameraInfo.name in policy.getStringArray("doEmptySectorQa"):
        maxMissing = policy.get("emptySectorMaxMissing")
        analysisList.append(qaAnalysis.EmptySectorQaAnalysis(maxMissing, nx = 4, ny = 4, useCache=keep,
                                                             wwwCache=wwwCache, delaySummary=delaySummary))
    if data.cameraInfo.name in policy.getStringArray("doAstromQa"):
        analysisList.append(qaAnalysis.AstrometricErrorQaAnalysis(policy.get("astromQaMaxErr"),
                                                                  useCache=keep, wwwCache=wwwCache,
                                                                  delaySummary=delaySummary))
    if data.cameraInfo.name in policy.getStringArray("doPhotCompareQa"):
        magCut   = policy.get("photCompareMagCut")
        deltaMin = policy.get("photCompareDeltaMin")
        deltaMax = policy.get("photCompareDeltaMax")
        rmsMax   = policy.get("photCompareRmsMax")
        derrMax  = policy.get("photCompareDerrMax")
        slopeMin = policy.get("photCompareSlopeMinSigma")
        slopeMax = policy.get("photCompareSlopeMaxSigma")
        for types in policy.getStringArray("photCompareTypes"):
            cmp1, cmp2 = types.split()
            starGxyToggle = False
            if types in policy.getStringArray("photCompareStarGalaxyToggle"):
                starGxyToggle = True
            analysisList.append(qaAnalysis.PhotCompareQaAnalysis(cmp1, cmp2, magCut, deltaMin, deltaMax,
                                                                 rmsMax, derrMax, slopeMin, slopeMax, starGxyToggle,
                                                                 useCache=keep,
                                                                 wwwCache=wwwCache,
                                                                 delaySummary=delaySummary))
    if data.cameraInfo.name in policy.getStringArray("doPsfShapeQa"):
        analysisList.append(qaAnalysis.PsfShapeQaAnalysis(policy.get("psfEllipMax"),
                                                          policy.get("psfFwhmMax"), useCache=keep,
                                                          wwwCache=wwwCache, delaySummary=delaySummary))
    if data.cameraInfo.name in policy.getStringArray("doCompleteQa"):
        analysisList.append(qaAnalysis.CompletenessQa(policy.get("completeMinMag"),
                                                      policy.get("completeMaxMag"), useCache=keep,
                                                      wwwCache=wwwCache, delaySummary=delaySummary))
    if data.cameraInfo.name in policy.getStringArray("doVignettingQa"):
        analysisList.append(qaAnalysis.VignettingQa(policy.get("vigMaxMedian"),
                                                    policy.get("vigMagRms"),
                                                    policy.get("vigMaxMag"), useCache=keep,
                                                    wwwCache=wwwCache, delaySummary=delaySummary))

    # allow command line override
    if (data.cameraInfo.name in policy.getStringArray("doVisitQa")) or doVisitQa:
        if matchDset == None and matchVisits == None:
            # we can't do it!
            print "Unable to run visit to visit Qa; please request a comparison visit or database"
            sys.exit(1)
        elif matchDset == None:
            matchDset = dataset
        elif matchVisits == None:
            matchVisits = []

        for mType in policy.getStringArray("visitQaMagTypes"):
            analysisList.append(qaAnalysis.VisitToVisitPhotQaAnalysis(matchDset, matchVisits, mType, 
                                                                      policy.get("visitQaMagCut"),
                                                                      policy.get("visitQaPhotDeltaMin"),
                                                                      policy.get("visitQaPhotDeltaMax"),
                                                                      policy.get("visitQaPhotRmsMax"),
                                                                      useCache=keep,
                                                                      wwwCache=wwwCache, 
                                                                      delaySummary=delaySummary))

        analysisList.append(qaAnalysis.VisitToVisitAstromQaAnalysis(policy.get("astromQaMaxErr"),
                                                                    matchDset, matchVisits, 
                                                                    useCache=keep,
                                                                    wwwCache=wwwCache, 
                                                                    delaySummary=delaySummary))


    # always run performance summary
    analysisList.append(qaAnalysis.performanceQa(useCache=keep, wwwCache=wwwCache,
                                                 delaySummary=delaySummary))
        
    # split by visit, and handle specific requests
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
            print "Can't run visits %d to %d as there are only %d visits" % (lo, hi, nvisit)
            sys.exit()
        if hi > nvisit:
            hi = nvisit
        visits = visits[lo:hi]

        print "Total of %d visits grouped by %d.  Running group %d with visits %d - %d:\n%s\n" % \
              (nvisit, groupSize, whichGroup, lo, hi-1, "\n".join(visits))
        groupTag = "%02d-%02d" % (groupSize, whichGroup)


    useFp = open("runtimePerformance%s.dat" % (groupTag), 'w')
    useFp.write("#%-11s %-24s %-32s %10s  %16s\n" %
                ("timestamp", "dataId", "testname", "t-elapsed", "resident-memory-kb-Mb"))
    

    # create progress tests for all visits
    if len(groupTag) > 0:
        groupTag = "."+groupTag
    progset = pipeQA.TestSet(group="", label="ProgressReport"+groupTag, wwwCache=wwwCache)
    for visit in visits:
        progset.addTest(visit, 0, [0, 1], "Not started.  Last dataId: None")

    testset = pipeQA.TestSet(group="", label="QA-failures"+groupTag, wwwCache=wwwCache)
    for visit in visits:

        visit_t0 = time.time()
        
        dataIdVisit = copy.copy(dataId)
        data.cameraInfo.setDataId(dataIdVisit, 'visit', visit)

        # now break up the run into eg. rafts or ccds
        #  ... if we only run one raft or ccd at a time, we use less memory
        brokenDownDataIdList = data.breakDataId(dataIdVisit, breakBy)

        for thisDataId in brokenDownDataIdList:
            
            for a in analysisList:

                test_t0 = time.time()
                test = str(a)
                if not re.search(testRegex, test) and not re.search('performance', test):
                    continue

                date = datetime.datetime.now().strftime("%a %Y-%m-%d %H:%M:%S")
                print "Running " + test + "  visit:" + str(visit) + "  ("+date+")"
                sys.stdout.flush() # clear the buffer before the fork

                # try the test() method
                t0 = time.time()
                tryThis(a.test, data, thisDataId, visit, test, testset, exceptExit)
                data.cachePerformance(thisDataId, test, "test-runtime", time.time() - t0)
                
                t0 = time.time()
                if forkFigure:
                    pid = os.fork()
                    if pid == 0:
                        # try the plot() method
                        tryThis(a.plot, data, thisDataId, visit, test, testset, exceptExit)
                        sys.exit(0)
                    else:
                        os.waitpid(pid, 0)
                else:
                    # try the plot() method
                    tryThis(a.plot, data, thisDataId, visit, test, testset, exceptExit)
                data.cachePerformance(thisDataId, test, "plot-runtime", time.time() - t0)
                
                # try the free() method
                tryThis(a.free, data, thisDataId, visit, test, testset, exceptExit)

                    
                memory = getMemUsageThisPid()
                test_tf = time.time()
                tstamp = time.mktime(datetime.datetime.now().timetuple())
                idstamp = ""
                for k,v in thisDataId.items():
                    idstamp += k[0]+str(v)
                useFp.write("%-12.1f %-24s %-32s %9.2fs %7d %7.2f\n" %
                            (tstamp, idstamp, test, test_tf-test_t0, memory, memory/1024.0))
                useFp.flush()

            # we're now done this dataId ... can clear the cache
            data.clearCache()
            
            raftName = ""
            if thisDataId.has_key('raft'):
                raftName = thisDataId['raft']+"-"
            ccdName = thisDataId[data.ccdConvention]
            progset.addTest(visit, 0, [1, 1], "Processing. Done %s%s." % (raftName,ccdName))
        progset.addTest(visit, 1, [1, 1], "Done processing.")

    useFp.close()


#############################################################
# end
#############################################################




if __name__ == '__main__':
    
    ########################################################################
    # command line arguments and options
    ########################################################################
    
    parser = optparse.OptionParser(usage=__doc__)
    parser.add_option("-b", "--breakBy", default="visit",
                      help="Break the run by 'visit','raft', or 'ccd' (default=%default)")
    parser.add_option("-C", "--camera", default=None,
                      help="Specify a camera and override auto-detection (default=%default)")
    parser.add_option("-c", "--ccd", default=".*",
                      help="Specify ccd as regex (default=%default)")
    parser.add_option("-d", "--delaySummary", default=False, action="store_true",
                      help="Delay making summary figures until all ccds are finished (default=%default)")
    parser.add_option("-e", "--exceptExit", default=False, action='store_true',
                      help="Don't capture exceptions, fail and exit (default=%default)")
    parser.add_option("-f", "--forkFigure", default=False, action='store_true',
                      help="Make figures in separate process (default=%default)")
    parser.add_option("-g", "--group", default=None,
                      help="Specify sub-group of visits to run 'groupSize:whichGroup' (default=%default)")
    parser.add_option("-k", "--keep", default=False, action="store_true",
                      help="Keep existing outputs (default=%default)")
    parser.add_option("-r", "--raft", default=".*",
                      help="Specify raft as regex (default=%default)")
    parser.add_option("-R", "--rerun", default=None,
                      help="Rerun to analyse - only valid for hsc/suprimecam (default=%default)")
    parser.add_option("-s", "--snap", default=".*",
                      help="Specify snap as regex (default=%default)")
    parser.add_option("-t", "--test", default=".*",
                      help="Regex specifying which QaAnalysis to run (default=%default)")
    parser.add_option("-V", "--verbosity", default=1,
                      help="Trace level for lsst.testing.pipeQA")
    parser.add_option("-v", "--visit", default=".*",
                      help="Specify visit as regex OR color separated list. (default=%default)")

    # visit-to-visit
    parser.add_option("--doVisitQa", default=False, action='store_true',
                      help="Do visit-to-visit pipeQA, overriding the policy default")
    parser.add_option("--matchDataset", default=None,
                      help="Specify another dataset to compare analysis to")
    parser.add_option("--matchVisits", default=None, action="append",
                      help="Visits within this dataset to compare analysis to; default is same visit")


    parser.add_option("--noWwwCache", default=False, action="store_true",
                      help="Disable caching of pass/fail (needed to run in parallel) (default=%default)")
    
    opts, args = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        sys.exit(1)
        
    #### handle visits

    # split by :
    visits = opts.visit
    if re.search(":", opts.visit):
        visits = opts.visit.split(":")

    dataId = {
        'visit': visits,
        'sensor': opts.ccd,
        'raft': opts.raft,
        'snap': opts.snap,
        }
        
    rerun = opts.rerun
    dataset, = args

    Trace.setVerbosity('lsst.testing.pipeQA', int(opts.verbosity))

    if not re.search("^(visit|raft|ccd)$", opts.breakBy):
        print "breakBy (-b) must be 'visit', 'raft', or 'ccd'"
        sys.exit()


    wwwCache = not opts.noWwwCache

    if re.search("(raft|ccd)", opts.breakBy) and not opts.keep:
        print "You've specified breakBy=%s, which requires 'keep' (-k). I'll set it for you."
        opts.keep = True
        
    main(dataset, dataId, rerun=rerun,
         doVisitQa = opts.doVisitQa, matchDset=opts.matchDataset, matchVisits=opts.matchVisits,
         testRegex=opts.test,          camera=opts.camera,
         exceptExit=opts.exceptExit,   keep=opts.keep,      wwwCache=wwwCache,
         breakBy=opts.breakBy, groupInfo=opts.group, delaySummary=opts.delaySummary,
         forkFigure=opts.forkFigure)
        
