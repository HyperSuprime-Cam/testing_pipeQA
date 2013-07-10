#!/usr/bin/env python
# Original filename: bin/pipeQa.py
#
# Author: 
# Email: 
# Date: Mon 2011-04-11 13:11:01
# 
# Summary: 
# 

import sys
import re
import optparse
import os
import datetime

import lsst.testing.pipeQA         as pipeQA
import lsst.testing.pipeQA.figures as qaFig

#############################################################
#
# Main body of code
#
#############################################################

def main():

    ts = pipeQA.TestSet("allFpaFigures", group="selftest", useCache=True, wwwCache=True)

    # need at least one test
    t = pipeQA.Test("dummy", 0, [-1, 1], "A dummy test")
    
    ts.addMetadata("keyname", "value")

    #cameraInfos = pipeQA.getCameraInfoAvailable()
    cameraInfos = [pipeQA.HscCameraInfo(), pipeQA.SuprimecamCameraInfo()] #LsstSimCameraInfo()]
    for camInfo in cameraInfos:
        print "trying camera: ", camInfo.name
        sfig = qaFig.FpaQaFigure(camInfo)
        vfig = qaFig.VectorFpaQaFigure(camInfo)
        
        # set values to range 0 to 1
        key0 = sfig.data.keys()[0]
        n = len(sfig.data.keys()) * len(sfig.data[key0].keys())
        i = 0
        for raft in sorted(sfig.data.keys()):
            ccdDict = sfig.data[raft]
            for ccd, value in ccdDict.items():
                sfig.data[raft][ccd] = 1.0*i
                vfig.data[raft][ccd] = [2.0*3.142*1.0*i/n, 1500*i/n, 1.0*i]
                i += 1
        
        sfig.makeFigure()
        ts.addFigure(sfig, camInfo.name+"_scalar.png", "Scalar FPA for camera: "+camInfo.name)
        vfig.makeFigure()
        ts.addFigure(vfig, camInfo.name+"_vector.png", "Vector FPA for camera: "+camInfo.name)

    ts.addTest(t)


if __name__ == '__main__':
    main()
    
