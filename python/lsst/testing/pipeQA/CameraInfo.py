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

import os, re
import copy
import eups
import numpy

import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
    


####################################################################
#
# CameraInfo base class
#
####################################################################


class CameraInfo(object):
    
    def __init__(self, name, dataInfo, camera=None):
        """
        @param name         A name for this camera
        @param dataInfo     List of dataId [name,distinguisher] pairs (distinguisher now depricated)
        @param camera       LSST camera object for this camera
        """
        
        self.name        = name
        self.dataInfo    = dataInfo
        self.camera      = camera

        self.rafts       = {}
        self.sensors     = {}
        self.detectors   = {}
        self.nSensor     = 0

        self.raftCcdKeys = []
        
        if self.camera is None:
            return

        self.rawName = "raw"

        self.ccdNamesBySerial = {}

        for r in self.camera:
            raft = cameraGeom.cast_Raft(r)
            raftName = raft.getId().getName().strip()
            self.detectors[raftName] = raft
            self.rafts[raftName] = raft
            for c in raft:
                ccd       = cameraGeom.cast_Ccd(c)
                ccdId     = ccd.getId()
                ccdName   = ccdId.getName().strip()
                ccdSerial = ccdId.getSerial()
                
                self.detectors[ccdName] = ccd
                self.sensors[ccdName] = ccd
                self.nSensor += 1
                self.raftCcdKeys.append([raftName, ccdName])
                self.ccdNamesBySerial[ccdSerial] = ccdName

                
        self.dataIdTranslationMap = {
            # standard  : camera
            'visit'  : 'visit',
            'snap'   : 'snap',
            'raft'   : 'raft',
            'sensor' : 'sensor',
            }
        self.dataIdDbNames = {
            'visit' : 'visit',
            'raft'  : 'raftName',
            'sensor' : 'ccdName',
            'snap'  : 'snap',
            }

    def standardToDbName(self, name):
        return self.dataIdDbNames[self.dataIdTranslationMap[name]]
        
    def dataIdCameraToStandard(self, dataIdIn):
        """Put this camera dataId in standard visit,raft,sensor format"""

        dataId = copy.copy(dataIdIn)
        
        transMap = self.dataIdTranslationMap

        for sKey, cKey in transMap.items():

            # it may be a simple string
            if isinstance(cKey, str):
                if dataId.has_key(cKey) and cKey != sKey:
                    dataId[sKey] = dataId[cKey]
                    del dataId[cKey]

            # or it may be a list e.g. visit might be coded ('-' sep.) to map to run-frame
            elif isinstance(cKey, list):
                haveKeys = True
                values = []
                for c in cKey:
                    if not dataId.has_key(c):
                        haveKeys = False
                    else:
                        values.append(dataId[c])

                if haveKeys:
                    dataId[sKey] = "-".join(map(str, values))
                    for c in cKey:
                        del dataId[c]
                else:
                    raise KeyError, "Can't map "+",".join(cKey)+" to "+sKey+". "+str(dataId)
                

        return dataId

        
        
    def dataIdStandardToCamera(self, dataIdIn):
        """ Convert an input dataId (visit,raft,sensor) to use key,value pairs for this specific camera. """

        dataId = copy.copy(dataIdIn)
        # check all keys in the translation map
        transMap = self.dataIdTranslationMap
        for inKey, outKey in transMap.items():

            # if we have this key, and it's has to be remapped ...
            if dataId.has_key(inKey) and inKey != outKey:

                # it may be a simple string
                if isinstance(outKey, str):
                    dataId[outKey] = dataId[inKey]
                    
                # or it may be a list e.g. visit might be coded ('-' sep.) to map to run-frame
                elif isinstance(outKey, list):
                    
                    # value for inKey must be '-' joined or '.*'
                    if dataId[inKey] == '.*':
                        dataId[inKey] = '-'.join(['.*']*len(outKey))
                    if not re.search('-', dataId[inKey]):
                        raise ValueError, "Combined keys must be dash '-' separated."
                    inValues = dataId[inKey].split('-')

                    # inValues must have the same number of values as outKey
                    if len(inValues) != len(outKey):
                        raise IndexError, "Can't map %s.  %d keys != %d values" % (inKey, len(outKey), len(inValues))
                    for i in range(len(outKey)):
                        dataId[outKey[i]] = inValues[i]
                        
                del dataId[inKey]
        return dataId


    def setDataId(self, dataId, key, value):
        if dataId.has_key(key):
            dataId[key] = value
        else:
            cKey = self.dataIdTranslationMap[key]
            if isinstance(cKey, str):
                dataId[cKey] = value
            elif isinstance(cKey, list):
                # make sure value can be split
                if not re.search("-", value):
                    raise ValueError, "Cannot split split "+value+" into "+str(cKey)

                values = value.split("-")
                if len(values) != len(cKey):
                    raise IndexError, "Number of values and keys don't match: "+ str(cKey)+" "+str(values)

                for i in range(len(cKey)):
                    dataId[cKey[i]] = values[i]
                    
        
    def getDetectorName(self, raft, ccd):
        ccdId = self.detectors[ccd].getId()
        name = re.sub("\s+", "_", ccdId.getName())
        serial = "%04d" % (ccdId.getSerial())
        return name + "--" + serial

    
    def getRoots(self, data, calib, output):
        """Store data directories in a dictionary
        
        @param data    Input directory containing registry.sqlite file
        @param calib   Input calibration directory containing calibRegistry.sqlite file
        @param output  Output directory
        """
        return {'data': data, 'calib': calib, 'output': output}


    def __str__(self):
        return self.name
    
    # some butler's use a 'rerun', others don't ... base class should default to None
    def getDefaultRerun(self):
        """Get our rerun."""
        return None

    def getSensorCount(self):
        """Get the number of sensors (ie. ccds) for this camera."""
        return self.nSensor
                

    def raftKeyToName(self, raft):
        for r in self.rafts.keys():
            if re.search(raft, r):
                return r
        return None

    def ccdKeyToName(self, ccd):
        for c in self.sensors.keys():
            if re.search("^R:\d,\d S:\d,\d$", c):
                if re.search("^R:\d,\d S:"+ccd, c):
                    return c
            else:
                if re.search(ccd, c):
                    return c
        return None

    def getRaftAndSensorNames(self, dataId):
        raftName = ''
        if self.dataIdTranslationMap.has_key('raft'):
            raftName = "R:"+str(dataId[self.dataIdTranslationMap['raft']])
        ccdName = raftName
        if self.dataIdTranslationMap.has_key('sensor'):
            ccdName += " S:"+str(dataId[self.dataIdTranslationMap['sensor']])
        return raftName, ccdName


    def getSensorName(self, raft, ccd):
        ccdName = raftName + " S:"+str(rowDict[self.cameraInfo.standardToDbName('sensor')])
        pass
    
    def getBbox(self, raftName, ccdName):

        rxc, ryc = 0.0, 0.0
        if self.rafts.has_key(raftName):
            raft = self.rafts[raftName]
            # NOTE: all ccd coords are w.r.t. the *center* of the raft, not its LLC
            rc   = raft.getCenter().getPixels(raft.getPixelSize())
            rxc     = rc.getX()
            ryc     = rc.getY()

        ccd   = self.sensors[ccdName]
        cc    = ccd.getCenter().getPixels(ccd.getPixelSize())
        cxc     = cc.getX()
        cyc     = cc.getY()
        orient  = ccd.getOrientation()
        nQuart  = ccd.getOrientation().getNQuarter()
        yaw     = orient.getYaw()

        cbbox   = ccd.getAllPixels(True)
        cwidth  = cbbox.getMaxX() - cbbox.getMinX()
        cheight = cbbox.getMaxY() - cbbox.getMinY()
        if abs(yaw.asRadians() - numpy.pi/2.0) < 1.0e-3:  # nQuart == 1 or nQuart == 3:
            ctmp = cwidth
            cwidth = cheight
            cheight = ctmp

        cx0     = rxc + cxc - cwidth/2
        cy0     = ryc + cyc - cheight/2
        cx1     = cx0 + cwidth
        cy1     = cy0 + cheight

        box = [cx0, cy0, cx1, cy1]
        return box

    




####################################################################
####################################################################    




####################################################################
#
# LsstSimCameraInfo class
#
####################################################################
class LsstSimCameraInfo(CameraInfo):
    
    def __init__(self):
        """ """
        dataInfo       = [['visit',1], ['snap', 0], ['raft',0], ['sensor',0]]

        if os.environ.has_key('OBS_LSSTSIM_DIR'):
            simdir        = os.environ['OBS_LSSTSIM_DIR']
            cameraGeomPaf = os.path.join(simdir, "description", "Full_STA_geom.paf")
            cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
            camera           = cameraGeomUtils.makeCamera(cameraGeomPolicy)
        else:
            camera = None
        
        CameraInfo.__init__(self, "lsstSim", dataInfo, camera)

        self.doLabel = False

        self.dataIdTranslationMap = {
            # input  : return
            'visit'  : 'visit',
            'snap'   : 'snap',
            'raft'   : 'raft',
            'sensor'    : 'sensor',
            }

        self.dataIdDbNames = {
            'visit' : 'visit',
            'raft'  : 'raftName',
            'sensor' : 'ccdName',
            'snap'  : 'snap',
            }

        

    
####################################################################
#
# CfhtCameraInfo class
#
####################################################################
class CfhtCameraInfo(CameraInfo):
    def __init__(self):
        
        dataInfo       = [['visit',1], ['ccd', 0]]
        
        if os.environ.has_key('OBS_CFHT_DIR'):
            simdir         = os.environ['OBS_CFHT_DIR']
            cameraGeomPaf = os.path.join(simdir, "megacam", "description", "Full_Megacam_geom.paf")
            cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
            camera           = cameraGeomUtils.makeCamera(cameraGeomPolicy)
        else:
            camera = None

        CameraInfo.__init__(self, "cfht", dataInfo, camera)
        self.doLabel = True
        

        

####################################################################
#
# HscCameraInfo class
#
####################################################################
class HscCameraInfo(CameraInfo):
    
    def __init__(self):
        """ """
        dataInfo       = [['visit',1], ['ccd', 0]]

        if os.environ.has_key('OBS_SUBARU_DIR'):
            simdir         = os.environ['OBS_SUBARU_DIR']
            cameraGeomPaf = os.path.join(simdir, "hsc", "hsc_geom.paf")
            if not os.path.exists(cameraGeomPaf):
                cameraGeomPaf = os.path.join(simdir, "hsc", "description", "hsc_geom.paf")
                if not os.path.exists(cameraGeomPaf):
                    raise Exception("Unable to find cameraGeom Policy file: %s" % (cameraGeomPaf))
            cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
            camera           = cameraGeomUtils.makeCamera(cameraGeomPolicy)
        else:
            camera = None
            
        CameraInfo.__init__(self, "hsc", dataInfo, camera)

        self.dataIdTranslationMap = {
            'visit'  : 'visit',
            'sensor' : 'ccd',
            }

        self.dataIdDbNames = {
            'visit'  : 'visit',
            'ccd'    : 'ccdname',
            }

        
        self.doLabel = False
        

    def getRaftAndSensorNames(self, dataId):
        raftName = ''
        ccd = dataId[self.dataIdTranslationMap['sensor']]
        ccdName = self.ccdNamesBySerial[int(ccd)]
        return raftName, ccdName
        

####################################################################
#
# SuprimecamCameraInfo class
#
####################################################################
class SuprimecamCameraInfo(CameraInfo):

    def __init__(self, mit=False):
        
        dataInfo       = [['visit',1], ['ccd', 0]]

        if os.environ.has_key('OBS_SUBARU_DIR'):
            simdir         = os.environ['OBS_SUBARU_DIR']
            pafBase = "Full_Suprimecam_MIT_geom.paf" if mit else "Full_Suprimecam_geom.paf"
            cameraGeomPaf = os.path.join(simdir, "suprimecam", pafBase)
                                         
            if not os.path.exists(cameraGeomPaf):
                raise Exception("Unable to find cameraGeom Policy file: %s" % (cameraGeomPaf))
            
            cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
            camera           = cameraGeomUtils.makeCamera(cameraGeomPolicy)
        else:
            camera           = None

        CameraInfo.__init__(self, "suprimecam", dataInfo, camera)

        self.doLabel = True
        self.mit = mit

        self.dataIdTranslationMap = {
            # standard  : camera
            'visit'  : 'visit',
            'sensor' : 'ccd',
            }
        self.dataIdDbNames = {
            'visit' : 'visit',
            'ccd'   : 'ccdname',
            }

    def getRaftAndSensorNames(self, dataId):
        
        raftName = ''
        ccd = dataId[self.dataIdTranslationMap['sensor']]
        ccdName = self.ccdNamesBySerial[int(ccd)]
        return raftName, ccdName
        
        
        



####################################################################
#
# SdssCameraInfo class
#
####################################################################
class SdssCameraInfo(CameraInfo):
    def __init__(self):
        
        messingWithNames = True
        if messingWithNames:
            dataInfo       = [['run', 1], ['filter', 0], ['field',1], ['camcol', 0]]
        else:
            dataInfo       = [['run', 1], ['band', 0], ['frame',1], ['camcol', 0]]

        if os.environ.has_key('OBS_SDSS_DIR') and obsSdss is not None:
            simdir         = os.environ['OBS_SDSS_DIR']
            camera = obsSdss.makeCamera.makeCamera(name='SDSS')
        else:
            camera           = None

        CameraInfo.__init__(self, "sdss", dataInfo, camera)
        self.rawName = "fpC"
        
        self.doLabel = True

        if messingWithNames:
            self.dataIdTranslationMap = {
                'visit' : ['run', 'field'],
                'raft'  : 'camcol',
                'sensor'   : 'filter',
                }

            self.dataIdDbNames = {
                'run' : 'run',
                'field' : 'field',
                'camcol' : 'camcol',
                'filter'   : 'filterName',
                }
        else:
            self.dataIdTranslationMap = {
                'visit' : ['run', 'frame'],
                'raft'  : 'camcol',
                'sensor'   : 'band',
                }

            self.dataIdDbNames = {
                'run' : 'run',
                'frame' : 'frame',
                'camcol' : 'camcol',
                'band'   : 'filterName',
                }
            

    def getRaftAndSensorNames(self, dataId):
        raftName = str(dataId[self.dataIdTranslationMap['raft']])
        ccdName =  str(dataId[self.dataIdTranslationMap['sensor']]) + raftName
        return raftName, ccdName

    




####################################################################
#
# SdssCameraInfo class
#
####################################################################
class CoaddCameraInfo(CameraInfo):

    def __init__(self):
        dataInfo       = [['tract', 1], ['patch', 1], ['filterName', 1]]

        simdir        = os.environ['TESTING_PIPEQA_DIR']
        cameraGeomPaf = os.path.join(simdir, "policy", "Full_coadd_geom.paf")
        cameraGeomPolicy = cameraGeomUtils.getGeomPolicy(cameraGeomPaf)
        camera           = cameraGeomUtils.makeCamera(cameraGeomPolicy)

        CameraInfo.__init__(self, "coadd", dataInfo, camera)
        
        self.doLabel = True

        self.dataIdTranslationMap = {
            'visit' : ['tract','patch','filterName'],
            'raft'  : None,
            'sensor'   : None, #'filter',
            }

        self.dataIdDbNames = {
            'patch' : 'patch',
            'tract' : 'tract',
            'filterName' : 'filterName',
            }
            

    def setFilterless(self):
        self.dataIdTranslationMap['visit'] = ['tract', 'patch']
        del self.dataIdDbNames['filterName']
        self.dataInfo = self.dataInfo[0:2]
        
    def getRaftAndSensorNames(self, dataId):
        ccdName =  'pseudo' #str(dataId['tract']) + '-' + str(dataId['patch'])
        return None, ccdName

    

    

    
    

def getCameraInfoAvailable():
    """Get a list of available CameraInfo objects."""
    
    available = []

    def tryLoad(cameraInfo):
        haveCam = True
        ci = cameraInfo()
        if ci.camera is None:
            haveCam = False
        return haveCam

    all = [
        SdssCameraInfo,
        CoaddCameraInfo,
        LsstSimCameraInfo,
        #CfhtCameraInfo,
        HscCameraInfo,
        SuprimecamCameraInfo,
        ]

    for camInfo in all:
        if tryLoad(camInfo):
            available.append(camInfo())

    return available
    
