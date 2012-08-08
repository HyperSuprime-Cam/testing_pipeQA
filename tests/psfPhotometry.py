#!/usr/bin/env python

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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

"""
Test to verify quality of PSF photometry on test frames
"""
import os
import unittest
import lsst.testing.pipeQA as pipeQA
import lsst.testing.displayQA as dispQA
import numpy

import lsst.afw.detection as afwDet

import matplotlib
import matplotlib.figure as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas

#import lsst.meas.extensions.shapeHSM as shapeHSM


def testPsfPhotometry():
    """Verify PSF photometry."""

    anData = "cpl_hsc_devel"
    anData = "hsc-dc1-2010-08-04.1-starsonly"
    anData = "imsim_20100625"
    #anData = "imsim_20100716"
    
    #########################
    # run the appropriate data set for this test
    pr = pipeQA.PipeRunner()
    if True:
        td1 = pipeQA.makeTestData("imsimTestData001",
                                  dataId={'visit':'85501867', 'snap':'0', 'raft':'1,1', 'sensor':'1,1'},
                                  verifyChecksum=False, astrometryNetData=anData,
                                  outDir='local')
        pr.addTestData(td1)
        td2 = pipeQA.makeTestData("imsimTestData001",
                                  dataId={'visit':'85502008', 'snap':'0', 'raft':'1,1', 'sensor':'1,1'},
                                  verifyChecksum=False, astrometryNetData=anData,
                                  outDir='local')
        pr.addTestData(td2)

        #hsmConfig = os.path.join(os.getenv('MEAS_EXTENSIONS_SHAPEHSM_DIR'), "policy", "hsmShape.paf")
        #qaConfig = os.path.join(os.getenv('TESTING_PIPEQA_DIR'), "policy", "lsstSim.paf")
        #pr.run(force=False, overrideConfig=[hsmConfig, qaConfig])
        pr.run(force=False, overrideConfig=[])


    ##########################
    # testing
    ts = dispQA.TestSet(mainDisplayFull=True)
    ts.importLogs(pr.getLogFiles())
    ts.importEupsSetups(pr.getEupsSetupFiles())
    ts.importExceptionDict(pr.getUncaughtExceptionDict())

    ##########################
    # get the data we want from the piperunner
    ss1 = pr.getSourceSet(dataId={'visit':'85501867', 'raft':'1,1'})
    ss2 = pr.getSourceSet(dataId={'visit':'85502008', 'raft':'1,1'})
    tolerance = numpy.radians(1.0)
    matchList = afwDet.matchRaDec(ss1, ss2, tolerance)

    
    ##########################
    # demo a plot
    
    # get the matches
    mag0 = numpy.array([])
    dmagAp = numpy.array([])
    dmagPsf = numpy.array([])
    x = numpy.array([])
    y = numpy.array([])
    for match in matchList:

        s0, s1 = match.first, match.second

        print "ixx: ", s0.getIxx()
        print "e1: ", s0.getE1(), "   shear1: ", s0.getShear1()
        
        m0 = -2.5*numpy.log10(s0.getApFlux())
        m1 = -2.5*numpy.log10(s1.getApFlux())
        dmagAp = numpy.append(dmagAp, m0 - m1)
        
        m0 = -2.5*numpy.log10(s0.getPsfFlux())
        m1 = -2.5*numpy.log10(s1.getPsfFlux())
        dmagPsf = numpy.append(dmagPsf, m0 - m1)

        dRa = s0.getRa() - s1.getRa()
        dDec = s0.getDec() - s1.getDec()

        x = numpy.append(x, s1.getXAstrom())
        y = numpy.append(y, s1.getYAstrom())
        mag0 = numpy.append(mag0, m0)


        
    ##########################
    # demo adding a test
    ts.addTest("test1", 1, [0, None], "verify 1 > 0")
    ts.addTest("test2", 1, [0, None], "verify something else")
    ts.addTest("test3", 1, [2, None], "verify yet something else")
    ts.addTest("test4", 1, [0, 2], "verify ;laksdkjr ")


        
    fig = figure.Figure(figsize=(6.0,6.0))
    canvas = FigCanvas(fig)
    axes = fig.add_subplot(1,1,1)
    axes.plot(mag0, dmagPsf, 'bo', markersize=1.0)
    axes.set_ylabel("$\Delta$mag")
    axes.set_xlabel("mag")
    axes.set_xlim([16.0,23.0])
    axes.set_ylim([-1.5,1.5])

    ts.addFigure(fig, "dmag.png", "Clearly a zeropoint problem here.")
    

    
if __name__ == "__main__":
    testPsfPhotometry()

