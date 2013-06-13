#!/usr/bin/env python
# Original filename: bin/pipeQa.py
#
# Author: 
# Email: 
# Date: Mon 2011-04-11 13:11:01
# 
# Summary: 
# 

"""
%prog [options]
"""

import sys
import re
import argparse
import os
import datetime

import numpy

import lsst.daf.persistence           as dafPersist
import lsst.testing.displayQA         as dispQA
import lsst.testing.pipeQA.figures as qaFig
import lsst.testing.pipeQA as pqa

import hsc.pipe.base.camera as hscCamera

import lsst.obs.hscSim               as hscSim

class QuickInfo(object):
    def __init__(self, filename, displayName):
        self.filename = filename
        self.displayName = displayName
        

#############################################################
#
# Main body of code
#
#############################################################

def main(root, rerun, visits, ccds, camera="hsc", testRegex=".*",
         doCopy=False, doConvert=False):

    datadir = os.path.join(root, rerun)
    if not os.path.exists(datadir):
        raise RuntimeError("Data directory does not exist:" + datadir)
    butler = dafPersist.Butler(datadir)
    if camera == 'suprimecam':
        camInfo = pqa.SuprimecamCameraInfo(mit=False)
    elif camera in ['hsc', 'hscsim']:
        camInfo = pqa.HscCameraInfo()
    else:
        raise ValueError, "Unknown camera."

    if False:
        detList = []
        for ccdName,ccd in camInfo.sensors.items():
            ix, iy = ccd.getCenterPixel()
            detList.append([ccdName, ix, iy])
        prev = None
        for arr in sorted(detList, key=lambda d: (d[2], d[1]))[::-1]:
            if prev and prev != arr[2]:
                print ""
            print arr[0], " ",
            prev = arr[2]
    
    figures = [
        QuickInfo("ossThumb",            "oss.Thumb"),
        QuickInfo("flattenedThumb",      "flattened.Thumb"),
        QuickInfo("plotEllipseGrid",     "ellipse.Grid"),
        QuickInfo("plotEllipticityMap",  "elly.Map"),
        QuickInfo("plotSeeingMap",       "seeing.Map"),
        QuickInfo("plotEllipseMap",      "ellipse.Map"),
        QuickInfo("plotFwhmGrid",        "fwhm.Grid"),
        QuickInfo("plotSeeingRobust",    "seeing.Robust"),
        QuickInfo("plotEllipticityGrid", "elly.Grid"),
        QuickInfo("plotMagHist",         "magHist"),
        QuickInfo("plotSeeingRough",     "seeing.Rough"),
        QuickInfo("plotPsfSrcGrid",      "psfCont.Src"),
        QuickInfo("plotPsfModelGrid",    "psfCont.Model"),
        ]
    
    
    for visit in visits:
        print "Running " + str(visit)
        
        for info in figures:

            f, flabel = info.filename, info.displayName

            if not re.search(testRegex, f):
                continue
            
            ts = dispQA.TestSet(group=str(visit), label=flabel, wwwCache=True)
            print "   ..."+f

            nobjs = {}
            for ccd in ccds:

                ####
                # butler.get() fails without filter and pointing ... why?
                dataId = {'visit': visit, 'ccd': ccd} #, 'filter':'W-S-I+', 'pointing' : 41}
                #### 

                
                raftName, ccdName = camInfo.getRaftAndSensorNames(dataId)
                areaLabel = camInfo.getDetectorName(raftName, ccdName)

                print "   CCD "+str(ccd)+ " " + areaLabel

                meta = butler.get('calexp_md', dataId)
                nobj = 1 #meta.get('NOBJ_MATCHED')
                nobjs[areaLabel] = nobj

                ts.addTest("nObj", nobj, [0.0, 1.0e4], "Number of matched objects", areaLabel=areaLabel)

                camel = f+"_filename"

                figNames = butler.get(camel, dataId)
                caption = "Figure " + f
                if os.path.exists(figNames[0]):
                    ts.addFigureFile(figNames[0], caption, areaLabel=areaLabel, doCopy=doCopy, doConvert=doConvert)
                else:
                    print "Warning: File does not exist ... "+ "\n".join(figNames)


            base = "fpa"+f
            dumData, dumMap = ts.unpickle(base, [None, None])
            fpa = qaFig.FpaQaFigure(camInfo, data=dumData, map=dumMap)

            for r, cDict in fpa.data.items():
                for c, value in cDict.items():
                    areaLabel = camInfo.getDetectorName(r, c)
                    if areaLabel in nobjs:
                        fpa.data[r][c] = nobjs[areaLabel]
                        fpa.map[r][c] = "%d" % (nobjs[areaLabel])

            fpa.makeFigure(vlimits=[0, 1.0e4], title=f, showUndefined=True,
                           failLimits=[0.0, 1.0e4], cmap="gist_heat_r")
            ts.addFigure(fpa, "fpa"+f+".png", "Figure of "+f, navMap=True)

            
def parseRange(s):
    values = []
    for ss in s.split(":"):
        ss_list = map(int, ss.split("-"))
        if len(ss_list) == 1:
            values += ss_list
        elif len(ss_list) == 2:
            values += range(ss_list[0], ss_list[1]+1)
    return values

        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("root",   help="Data root directory.")
    parser.add_argument("rerun",  help="Rerun for the data to be run.")
    parser.add_argument("visits", help="Visits to run.  Use colon to separate values.  Use dash to denote range.")
    parser.add_argument("-C",  "--camera", default="hsc",
                        help="Specify camera")
    parser.add_argument("-c",  "--ccds", default="0-9",
                        help="Ccds to run.  Use colon to separate values.  Use dash to denote range.")
    parser.add_argument("-t", "--testRegex", default=".*",
                        help="regular expression to match tests to be run")
    parser.add_argument("--do-copy", dest="doCopy", action="store_true", help="if set, make copies of files onto pipeQa's web direcotry ($WWW_ROOT/$WWW_RERUN) rather than making symlinks.")
    parser.add_argument("--do-convert", dest="doConvert", action="store_true", help="if set, make thumbnail files by converting original png files. Otherwise, just use symlinks.")
    
    args = parser.parse_args()

    visits = parseRange(args.visits)
    ccds   = parseRange(args.ccds)
    
    main(args.root, args.rerun, visits, ccds, testRegex=args.testRegex, camera=args.camera,
         doCopy=args.doCopy, doConvert=args.doConvert)
