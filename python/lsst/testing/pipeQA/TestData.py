import sys, os, glob, re, stat
import traceback

import sqlite

import eups
import lsst.pex.policy                  as pexPolicy
import lsst.pex.logging                 as pexLog
import lsst.daf.persistence             as dafPersist
from lsst.testing.pipeQA.Checksum       import Checksum
from lsst.testing.pipeQA.Manifest       import Manifest

import QaDataUtils as qaDataUtils

#import lsst.meas.extensions.shapeHSM.hsmLib as shapeHSM

    
#######################################################################
#
#
#
#######################################################################
class TestData(object):
    """ """

    #######################################################################
    #
    #######################################################################
    def __init__(self, label, rerun, mapperClass, dataInfo, defaultConfig, kwargs):
        """
        keyword args:
        haveManifest = boolean, verify files in dataDir are present according to manifest
        verifyChecksum = boolean, verify files in dataDir have correct checksum as listed in manifest
        astrometryNetData = eups package name for astrometryNetData package to use
        """

        ###############################################
        # handle inputs
        ###############################################
        self.label         = label
        self.dataIdNames   = []
        self.dataIdDiscrim = []
        self.defaultConfig = defaultConfig
        roots = self.defaultConfig['roots']

        self.kwargs         = kwargs
        self.dataIdRegex    = self.kwargs.get('dataId', {})
        self.haveManifest   = self.kwargs.get('haveManifest', False)
        self.verifyChecksum = self.kwargs.get('verifyChecksum', False)
        self.astrometryNetData = self.kwargs.get('astrometryNetData', None)

        self.rerun = rerun
        
        ##################
        # output directory

        # if the user provided one, use it ... otherwise use the default
        self.outDir = kwargs.get('outDir', roots['output'])
        self.testdataDir = roots['data']
        self.calibDir    = roots['calib']

        # need a separate output dir for things we generate outside the pipe (logs, etc)
        self.localOutDir = os.path.join(os.getcwd(), self.label+"out")
        if not os.path.exists(self.localOutDir):
            os.mkdir(self.localOutDir)
        
        # allow a short hand for 'write outputs locally' ... use the word 'local'
        if re.search("(local|\.)", self.outDir):
            self.outDir = self.localOutDir
        if not os.path.exists(self.outDir):
            os.mkdir(self.outDir)
        roots['output'] = self.outDir


        # This (dataId fetching) needs a better design, but will require butler/mapper change, I think.
        #
        # these obscure things refer to the names assigned to levels in the data hierarchy
        # eg. for lsstSim:   dataInfo  = [['visit',1], ['snap', 0], ['raft',0], ['sensor',0]]
        # a level is considered a discriminator if it represents different pictures of the same thing
        # ... so the same object may appear in multiple 'visits', but not on multiple 'sensors'
        # dataInfo is passed in from the derived class as it's specific to each mapper
        for array in dataInfo:
            dataIdName, dataIdDiscrim = array
            self.dataIdNames.append(dataIdName)
            self.dataIdDiscrim.append(dataIdDiscrim)

        # keep a list of any:
        #  - logfiles we write,
        #  - eups setup list files
        # ...so we can provide it to be imported by a TestSet
        self.logDir = os.path.join(self.localOutDir, "log")
        self.logFiles = []
        self.eupsSetupFiles = []
        

        
        ##########################
        # load the manifest and verify the checksum (if we're asked to ... it's slower)
        # haveManifest = True is a bit slowish
        # verifyChecksum = True is quite slow
        if self.haveManifest:
            manifest = Manifest(self.testdataDir)
            manifest.read()
            missingInputs   = manifest.verifyExists()
            if self.verifyChecksum:
                failedChecksums = manifest.verifyChecksum()

            msg = ""
            if (len(missingInputs) > 0):
                msg = "Missing input files listed in manifest:\n"
                msg += "\n".join(missingInputs) + "\n"
            if self.verifyChecksum and (len(failedChecksums) > 0):
                msg += "Failed checksums:\n"
                msg += "\n".join(failedChecksums) + "\n"
            if len(msg) > 1:
                raise Exception(msg)

                    

        #########################
        # see if setup changed
        # we should rerun our data if the user has setup against different packages
        print "Warning: Setup change verification not yet implemented."
            

        #######################################
        # get i/o butlers
        registry = os.path.join(self.testdataDir, 'registry.sqlite3')
        # note: only suprime/hsc mappers accept rerun arg
        if self.rerun is None:
            self.inMapper  = mapperClass(root=self.testdataDir, calibRoot=self.calibDir)
            self.outMapper = mapperClass(root=self.outDir, calibRoot=self.calibDir, registry=registry)
        else:
            self.inMapper  = mapperClass(rerun=self.rerun, root=self.testdataDir, calibRoot=self.calibDir)
            self.outMapper = mapperClass(rerun=self.rerun,
                                         root=self.outDir, calibRoot=self.calibDir, registry=registry)
        self.inButler  = dafPersist.ButlerFactory(mapper=self.inMapper).create()
        self.outButler = dafPersist.ButlerFactory(mapper=self.outMapper).create()

        
        ####################################################
        # make a list of the frames we're asked to care about

        # get all the available raw inputs
        self.availableDataTuples = self.inButler.queryMetadata('raw', self.dataIdNames,
                                                               format=self.dataIdNames)

        # of the data available, get a list of the ones the user actually wants us
        #  to run.  A bit sketchy here ... kwargs contains non-idname info as well.
        self.dataTuples = self._regexMatchDataIds(self.dataIdRegex, self.availableDataTuples)


        # if/when we run, we'll store tracebacks for any failed runs
        # we don't want to stop outright, but we should report failures
        self.uncaughtExceptionDict = {}


        
    #######################################################################
    # setup our own astrometryNetData package, if we have one
    #######################################################################
    def setupAstrometryNetData(self):
        
        if not (self.astrometryNetData is None):
            ok, version, reason = eups.Eups().setup('astrometry_net_data', self.astrometryNetData)
        
        
    #######################################################################
    # Run our data through pipette
    #######################################################################
#    def run(self, kwargs):
#        """Run pipette on the data we know about."""
#        
#        force             = kwargs.get('force', False)
#        overrideConfigs   = kwargs.get('overrideConfig', None)  # array of paf filenames
#
#        # setup a specific astromentry.net data package, if one is provided
#        self.setupAstrometryNetData()
#        
#            
#        # keep a record of the eups setups for the data we're running
#        eupsSetupFile = os.path.join(self.localOutDir, self.label+".eups")
#        self.eupsSetupFiles.append(eupsSetupFile)
#        fp = open(eupsSetupFile, 'w')
#        ups = eups.Eups()
#        products = ups.getSetupProducts()
#        for product in products:
#            fp.write("%s %s\n" % (product.name, product.version))
#        fp.close()
#
#        
#        # merge in override config
#        config = self.defaultConfig
#        if overrideConfigs is not None:
#            for overrideConfig in overrideConfigs:
#                config.merge(pipette.config.Config(overrideConfig))
#                #config = pipette.config.configuration(config, overrideConfig)
#
#        #srcConf = config['measure']['source']
#        #srcConf['shape'] = "HSM_BJ"
#
#        #shapeConf = config['measure']['shape']
#        #shapeConf['HSM_BJ'] = pexPolicy.Policy()
#        #shapeConf['HSM_BJ']['enabled'] = True
#        
#        #do = config['do']
#        #do['phot'] = True
#        #do['ast']  = True
#        #do['cal']  = True
#
#        if len(self.dataTuples) == 0:
#            sys.stderr.write(
#                "WARNING: Requested data not found in registry.  Available frames:\n" +
#                str(self.availableDataTuples)
#                )
#            
#
#        for dataTuple in self.dataTuples:
#
#            # put these values in a Dict with the appropriate keys
#            dataId = self._tupleToDataId(dataTuple)
#            
#            # see if we already have the outputs
#            isWritten = self.outButler.datasetExists('src', dataId)
#
#            # set isWritten if we're a testBot ... assume we shouldn't try to run
#            #  ... can always override with 'force'
#            if re.search("^testBot", self.label):
#                isWritten = True
#            
#            thisFrame = "%s=%s" % (",".join(self.dataIdNames), str(dataTuple))
#
#            if force or (not isWritten):
#                
#                # create a log that prints to a file
#                idString = self.dataTupleToString(dataTuple)
#                if not os.path.exists(self.logDir):
#                    os.mkdir(self.logDir)
#                logFile = os.path.join(self.logDir, idString+".log")
#
#                if os.path.exists(logFile):
#                    os.remove(logFile)
#                    
#                self.logFiles.append(logFile)
#
#                log = pexLog.Log.getDefaultLog().createChildLog("testQA.TestData", pexLog.Log.INFO)
#                log.addDestination(pexLog.FileDestination(logFile, True))
#
#                # run, if necessary
#                
#                print "Running:  %s" % (thisFrame)
#                try:
#                    self.runPipette(self.rerun, dataId, config, log)
#                except Exception, e:
#                    exc_type, exc_value, exc_traceback = sys.exc_info()
#                    s = traceback.format_exception(exc_type, exc_value,
#                                                   exc_traceback)
#                    self.uncaughtExceptionDict[idString] = "Running "+idString+"\n" + "".join(s)
#
#                    log.log(log.WARN, idString + ": Unrecoverable exception. (see traceback) - "+str(e))
#
#            else:
#                print "%s exists, skipping. (use force=True to force a run)"  % (thisFrame) 

                
    #######################################################################
    #
    #######################################################################
    def getLogFiles(self):
        """Get all the logfiles in our log directory."""
        pattern = os.path.join(self.logDir, "*.log")
        return glob.glob(pattern)
    

    #######################################################################
    #
    #######################################################################
    def getEupsSetupFiles(self):
        """Get all the eups setup files in our log directory."""
        pattern = os.path.join(self.localOutDir, "*.eups")
        return glob.glob(pattern)

    
    #######################################################################
    #
    #######################################################################
    def getUncaughtExceptionDict(self):
        """
        Get any exceptions thrown from within Pipette.
        Return as a dictionary with dataId values as keys.
        """
        return self.uncaughtExceptionDict
    
    #######################################################################
    #
    #######################################################################
    def dataTupleToString(self, dataTuple):
        """Represent a dataTuple as a string."""
        
        s = []
        for i in range(len(self.dataIdNames)):
            name = self.dataIdNames[i]
            value = re.sub("[,]", "", str(dataTuple[i]))
            s.append(name + value)
        return "-".join(s)


    #######################################################################
    #
    #######################################################################
    def setAstrometryNetData(astrometryNetData):
        """Accessor for astrometryNetData package we were asked to use."""
        self.astrometryNetData = astrometryNetData


    
    #######################################################################
    #
    #######################################################################
    def getSourceSet(self, dataIdRegex):
        """Get sources for requested data as one sourceSet."""
        dataTuplesToFetch = self._regexMatchDataIds(dataIdRegex, self.dataTuples)

        # get the datasets corresponding to the request
        sourceSet = []
        for dataTuple in dataTuplesToFetch:
            dataId = self._tupleToDataId(dataTuple)

            # make sure we actually have the output file
            isWritten = self.outButler.datasetExists('src', dataId)
            if isWritten:
                persistableSourceVector = self.outButler.get('src', dataId)
                sourceSetTmp = persistableSourceVector.getSources()

                if self.outButler.datasetExists('calexp', dataId):
                    postIsrCcd = self.outButler.get('calexp', dataId)
                    calib = postIsrCcd.getCalib()
                    
                    fmag0, fmag0err = calib.getFluxMag0()
                    for s in sourceSetTmp:
                        apFlux  = s.getApFlux()
                        psfFlux = s.getPsfFlux()
                        s.setApFlux(apFlux/fmag0)
                        s.setPsfFlux(psfFlux/fmag0)
                
                sourceSet += sourceSetTmp
            else:
                print str(dataTuple) + " output file missing.  Skipping."
                
        return sourceSet
            

    


    #######################################################################
    # utility to go through a list of data Tuples and return
    #  the ones which match regexes for the corresponding data type
    # so user can say eg. raft='0,\d', visit='855.*', etc
    #######################################################################
    def _regexMatchDataIds(self, dataIdRegexDict, availableDataTuples):
        """ """

        # if it's been requested as a generic 'dataSet'
        if dataIdRegexDict.has_key('dataSet'):
            availableDataTuplesGrouped = {}
            for dataTuple in availableDataTuples:
                discrim = ""
                for i in range(len(dataTuple)):
                    if self.dataIdDiscrim[i]:
                        discrim += str(dataTuple[i])
                if availableDataTuplesGrouped.has_key(discrim):
                    availableDataTuplesGrouped[discrim].append(dataTuple)
                else:
                    availableDataTuplesGrouped[discrim] = [dataTuple]

            idxReq = dataIdRegexDict['dataSet']
            if not isinstance(idxReq, list):
                idxReq = [idxReq]
            dataTuples = []
            for idx in idxReq:
                if len(availableDataTuplesGrouped.keys()) > idx:
                    discrimRequested = sorted(availableDataTuplesGrouped.keys())[idx]
                    dataTuples += availableDataTuplesGrouped[discrimRequested]
            return dataTuples
                    
        else:
            # go through the list of what's available, and compare to what we're asked for
            # Put matches in a list of tuples, eg. [(vis1,sna1,raf1,sen1),(vis2,sna2,raf2,sen2)] 
            dataTuples = []
            for dataTuple in availableDataTuples:

                # start true and fail if any dataId keys fail ... eg. 'visit' doesn't match
                match = True
                for i in range(len(self.dataIdNames)):
                    dataIdName = self.dataIdNames[i]   # eg. 'visit', 'sensor', etc
                    regexForThisId = dataIdRegexDict.get(dataIdName, '.*') # default to '.*' or 'anything'
                    dataId = dataTuple[i]

                    # if it doesn't match, this frame isn't to be run.
                    if not re.search(regexForThisId,  str(dataId)):
                        match = False

                if match:
                    dataTuples.append(dataTuple)
                
        return dataTuples
                

    
    #######################################################################
    # utility to convert a data tuple to a dictionary using dataId keys
    #######################################################################
    def _tupleToDataId(self, dataTuple):
        """ """
        dataId = {}
        for i in range(len(self.dataIdNames)):
            dataIdName = self.dataIdNames[i]
            dataId[dataIdName] = dataTuple[i]
        return dataId


    

    
    
#######################################################################
#
#
#
#######################################################################
class ImSimTestData(TestData):
    """ """
    
    #######################################################################
    #
    ####################################################################### 
    def __init__(self, label, **kwargs):
        """ """

        import lsst.obs.lsstSim           as obsLsst
        import lsstSim
        
        mapper         = obsLsst.LsstSimMapper
        dataInfo       = [['visit',1], ['snap', 0], ['raft',0], ['sensor',0]]
        
        # find the label in the testbed path
        testbedDir, testdataDir = qaDataUtils.findDataInTestbed(label)
        
        defaultConfig   = lsstSim.getConfig()
        if not 'roots' in defaultConfig:
            defaultConfig['roots'] = pexPolicy.Policy()
        roots           = defaultConfig['roots']
        roots['data']   = testdataDir
        roots['calib']  = testdataDir
        roots['output'] = testdataDir

        rerun = None
        TestData.__init__(self, label, rerun, mapper, dataInfo, defaultConfig, kwargs)

        
    #######################################################################
    #
    #######################################################################
    def runPipette(self, rerun, dataId, config, log):
        """ """
        import lsstSim
        lsstSim.run(rerun, dataId['visit'], dataId['snap'], dataId['raft'], dataId['sensor'],
                    config, log=log)


#######################################################################
#
#
#
#######################################################################
class HscSimTestData(TestData):
    """ """

    
    #######################################################################
    #
    #######################################################################
    def __init__(self, label, **kwargs):
        """ """
        
        import lsst.obs.hscSim            as obsHsc
        import runHsc
        mapper         = obsHsc.HscSimMapper
        dataInfo       = [['visit',1], ['ccd', 0]]
        
        # find the label in the testbed path
        testbedDir, testdataDir = qaDataUtils.findDataInTestbed(label)

        #defaultConfig = runHsc.getConfigQa()
        argv = []
        argv.extend(["--instrument=hsc", "--frame=0", "--ccd=0", "--rerun=test"])
        defaultConfig, opts, args = runHsc.getConfig(argv=argv)

        if not 'roots' in defaultConfig:
            defaultConfig['roots'] = pexPolicy.Policy()
        roots           = defaultConfig['roots']
        roots['data']   = os.path.join(testdataDir, "HSC")
        roots['calib']  = os.path.join(testdataDir, "CALIB")
        roots['output'] = os.path.join(testdataDir, "HSC")

        rerun = "pipeQA"
        TestData.__init__(self, label, rerun, mapper, dataInfo, defaultConfig, kwargs)


    #######################################################################
    #
    #######################################################################
    def runPipette(self, rerun, dataId, config, log):
        """ """
        import runHsc
        #runHsc.run(rerun, dataId['visit'], dataId['ccd'], config, log=log)
        runHsc.doRun(rerun=rerun, frameId=dataId['visit'], ccdId=dataId['ccd'],
                   doMerge=True, doBreak=False,
                   instrument="hsc",
                   output =config['roots']['output'],
                   calib  =config["roots"]['calib'],
                   data   =config["roots"]['data'],
                   log    = log
        )

    

#######################################################################
#
#
#
#######################################################################
class SuprimeTestData(TestData):
    """ """

    
    #######################################################################
    #
    #######################################################################
    def __init__(self, label, **kwargs):
        """ """
        
        import lsst.obs.suprimecam        as obsSuprimecam
        import suprimecam
        mapper         = obsSuprimecam.SuprimecamMapper
        dataInfo       = [['visit',1], ['ccd', 0]]
        
        # find the label in the testbed path
        testbedDir, testdataDir = qaDataUtils.findDataInTestbed(label)

        defaultConfig   = suprimecam.getConfig()
        if not 'roots' in defaultConfig:
            defaultConfig['roots'] = pexPolicy.Policy()
        roots = defaultConfig['roots']            
        roots['data']   = os.path.join(testdataDir, "SUPA")
        roots['calib']  = os.path.join(testdataDir, "SUPA", "CALIB")
        roots['output'] = os.path.join(testdataDir, "SUPA")

        rerun = "pipeQA"
        
        TestData.__init__(self, label, rerun, mapper, dataInfo, defaultConfig, kwargs)


    #######################################################################
    #
    #######################################################################
    def runPipette(self, rerun, dataId, config, log):
        """ """
        import suprimecam
        suprimecam.run(rerun, dataId['visit'], dataId['ccd'], config, log=log)



        
#######################################################################
#
#
#
#######################################################################
class CfhtTestData(TestData):
    """ """
    
    #######################################################################
    #
    #######################################################################
    def __init__(self, label, **kwargs):
        """ """

        import lsst.obs.cfht              as obsCfht
        import megacam
        mapper         = obsCfht.CfhtMapper
        dataInfo       = [['visit',1], ['ccd', 0]]
        
        # find the label in the testbed path
        testbedDir, testdataDir = qaDataUtils.findDataInTestbed(label)

        defaultConfig   = megacam.getConfig()
        if not 'roots' in defaultConfig:
            defaultConfig['roots'] = pexPolicy.Policy()
        roots           = defaultConfig['roots']
        roots['data']   = testdataDir
        roots['calib']  = os.path.join(testdataDir, "calib")
        roots['output'] = testdataDir

        rerun = None
        
        TestData.__init__(self, label, rerun, mapper, dataInfo, defaultConfig, kwargs)


    #######################################################################
    #
    #######################################################################
    def runPipette(self, rerun, dataId, config, log):
        """ """
        import megacam
        megacam.run(rerun, dataId['visit'], dataId['ccd'], config, log=log)


        



        

#######################################################################
#
#
#
#######################################################################
def makeTestData(label, **kwargs):
        
    testbedDir, testdataDir = qaDataUtils.findDataInTestbed(label)
        
    regFile = 'registry.sqlite3'
    registry = os.path.join(testdataDir, regFile)
    cfhtCalibRegistry = os.path.join(testdataDir, "calib", "calibRegistry.sqlite3")
    hscSimCalibRegistry = os.path.join(testdataDir, "CALIB", "calibRegistry.sqlite3")

    testCfht       = os.path.exists(cfhtCalibRegistry)
    testHscSim     = os.path.exists(hscSimCalibRegistry)
    testSuprimecam = os.path.exists(os.path.join(testdataDir, "SUPA"))
    
    # define some tests to distinguish which type of data we have
    lookup = {
        "lsstSim" : [ImSimTestData,  not testCfht and not testHscSim and not testSuprimecam],
        "cfht" :    [CfhtTestData, testCfht],
        "hscSim"  : [HscSimTestData, testHscSim], 
        "suprime":  [SuprimeTestData, testSuprimecam],
        }


    ######################################################
    # Do our best to figure out what we've been handed.
    validLookups = {}
    for key, array in lookup.items():

        TestDataClass, distinguishTest = array

        if distinguishTest:
            validLookups[key] = array

    #print testbedDir, testdataDir
    
    nValid = len(validLookups.keys())
    if nValid > 1:
        raise Exception("Registries consistent with multiple mappers (" +
                        ",".join(validLookups.keys()) + ").  Can't decide which to use.")
    elif nValid == 0:
        raise Exception("Can't find registries usable with any mappers.")
    else:
        key = validLookups.keys()[0]
        TestDataClass, distinguishTest = validLookups[key]
        return TestDataClass(label, **kwargs)
