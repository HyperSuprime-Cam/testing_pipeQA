import sys, os, re

from  QaDataUtils import QaDataUtils
import CameraInfo as qaCamInfo


###################################################
# Factory for QaData
###################################################
def makeQaData(label, rerun=None, retrievalType=None, camera=None, **kwargs):
    """Factory to make a QaData object for either Butler data, or Database data.

    @param label         data identifier - either a directory in TESTBED_PATH/SUPRIME_DATA_DIR or a DB name
    @param rerun         data rerun to retrieve
    @param retrievalType 'butler', 'db', or None (will search first for butler, then database)
    @param camera        Specify which camera is to be used
    """
    
    print "RetrievalType=", retrievalType
    print "camera=", camera

    cameraInfos = {
        # "cfht": qaCamInfo.CfhtCameraInfo(), # CFHT camera geometry broken following #1767
        "hsc"            : [qaCamInfo.HscCameraInfo,        []],
        "suprimecam"     : [qaCamInfo.SuprimecamCameraInfo, []],
        "suprimecam-mit" : [qaCamInfo.SuprimecamCameraInfo, [True]],
        "sdss"           : [qaCamInfo.SdssCameraInfo,       []],
        "coadd"          : [qaCamInfo.CoaddCameraInfo,      []],
        "lsstSim"        : [qaCamInfo.LsstSimCameraInfo,    []],
        }


    cameraToUse = None
    
    # default to lsst
    if camera is None:
        cam, args = cameraInfos['lsstSim']
        cameraToUse = cam(*args)
        camera = 'lsstSim'
    else:
        cam, args = cameraInfos[camera]
        cameraToUse = cam(*args)
            
    # we should never get here as we default to LSST
    if cameraToUse is None:
        raise RuntimeError("Can't load camera:" + str(camera))

    
    #####################
    # make a butler QaData
    if retrievalType.lower() == "butler":
        
        qaDataUtils = QaDataUtils()
        testbedDir, testdataDir = qaDataUtils.findDataInTestbed(label)
        from ButlerQaData  import ButlerQaData
        print "label:       ", label
        print "rerun:       ", rerun
        print "camera:      ", cameraToUse.name
        print "testdataDir: ", testdataDir
        return ButlerQaData(label, rerun, cameraToUse, testdataDir, **kwargs)

    
    #####################
    # make a db QaData
    if retrievalType.lower() in ['db', 'database']:

        if camera in ["hsc","suprimecam","suprimecam-mit"]:
            from HscDbQaData      import HscDbQaData
            return HscDbQaData(label, rerun, cameraToUse)
        else:
            from DbQaData         import DbQaData
            return DbQaData(label, rerun, cameraToUse)




