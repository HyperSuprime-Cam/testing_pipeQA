import sys, os, re

from  QaDataUtils import QaDataUtils


###################################################
# Factory for QaData
###################################################
def makeQaData(label, rerun=None, retrievalType=None, camera=None, **kwargs):
    """Factory to make a QaData object for either Butler data, or Database data.

    @param label         identifier for the data - either a directory in TESTBED_PATH/SUPRIME_DATA_DIR or a DB name
    @param rerun         data rerun to retrieve
    @param retrievalType 'butler', 'db', or None (will search first for butler, then database)
    @param camera        Specify which camera is to be used
    """
    
    if retrievalType is None:
        
        # if there's only one possibility, use that
        
        # see if there's a testbed directory called 'label'
        # if TESTBED_PATH isn't set, skip this and assume it's a db
        validButler = False
        if os.environ.has_key('TESTBED_PATH') or os.environ.has_key('SUPRIME_DATA_DIR'):
            qaDataUtils = QaDataUtils()            
            testbedDir, testdataDir = qaDataUtils.findDataInTestbed(label, raiseOnFailure=False)
            if (not testbedDir is None) and (not testdataDir is None):
                validButler = True

        # see if we can connect to a database with name 'label'
        # NOTE: must update if/when non-lsst databases get used
        validDb = True

        try:
            if not camera is None and  camera.lower() in ["suprimecam", "hsc"]:
                from HscDatabaseQuery import DbInterface, DatabaseIdentity
                identity = DatabaseIdentity(label)
                dbInterface = DbInterface(identity)              
            else:
                from DatabaseQuery import LsstSimDbInterface, DatabaseIdentity
                dbInterface = LsstSimDbInterface(DatabaseIdentity(label))

        # database not found will throw this error
        except OperationalError, e:
            validDb = False
            
        # raise other exceptions (like if we forgot to setup mysqlclient pyscopg2
        except Exception, e:
            validDb = False

            
        if validButler and not validDb:
            retrievalType = 'butler'
        if validDb and not validButler:
            retrievalType = 'db'
        if validDb and validButler:
            raise Exception("The label "+label+" is present as both a testbed directory and a database."\
                            "Please specify retrievalType='butler', or retrievalType='db'.")
        if not validDb and not validButler:
            print "Unable to find "+label+" as a testbed directory or a database."
            raise


    print "RetrievalType=", retrievalType
    print "camera=", camera
    
    if re.search("^[Bb]utler$", retrievalType):
        from ButlerQaData  import makeButlerQaData
        return makeButlerQaData(label, rerun, camera=camera, **kwargs)
    
    if re.search("^([Dd][Bb]|[Dd]atabase)$", retrievalType):


        import CameraInfo as qaCamInfo
        cameraInfos = {
    #       "cfht": qaCamInfo.CfhtCameraInfo(), # XXX CFHT camera geometry is currently broken following #1767
            "hsc"            : [qaCamInfo.HscCameraInfo, []],
            "suprimecam"     : [qaCamInfo.SuprimecamCameraInfo, []],
            "suprimecam-old" : [qaCamInfo.SuprimecamCameraInfo, [True]],
            "sdss"           : [qaCamInfo.SdssCameraInfo, []],
            "coadd"          : [qaCamInfo.CoaddCameraInfo, []],
            "lsstSim"        : [qaCamInfo.LsstSimCameraInfo, []],
            }


        cameraToUse = None
        if not camera is None:
            cam, args = cameraInfos[camera]
            cameraToUse = cam(*args)
        else:
            cam, args = cameraInfos['lsstSim']
            cameraToUse = cam(*args)
            camera = 'lsstSim'
            
        if re.search("^(hsc|suprimecam|suprimecam-old)$", camera):
            from HscDbQaData      import HscDbQaData
            return HscDbQaData(label, rerun, cameraToUse)
        else:
            from DbQaData      import DbQaData
            return DbQaData(label, rerun, cameraToUse)




