#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

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
        
        if os.environ.has_key('TESTBED_PATH'):
            testdataDir = os.path.join(os.getenv("TESTBED_PATH"), label)
        else:
            raise Exception("Must specify environment variable TESTBED_PATH.")

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
            return HscDbQaData(label, rerun, cameraToUse, **kwargs)
        else:
            from DbQaData         import DbQaData
            return DbQaData(label, rerun, cameraToUse, **kwargs)




