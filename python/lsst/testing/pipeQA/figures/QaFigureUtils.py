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

import numpy
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.patches import Ellipse
import matplotlib.cm as cm


import lsst.afw.cameraGeom as cameraGeom


def cameraToRectangles(camera):
    rectangles = {}
    centers = {}
    raftBoundaries = []
    ccdBoundaries = {}

    
    for r in camera:
        raft = cameraGeom.cast_Raft(r)

        # NOTE: all ccd coords are w.r.t. the *center* of the raft, not its LLC
        rc     = raft.getCenter().getPixels(raft.getPixelSize())
        rxc, ryc = rc.getX(), rc.getY()

        xmin = +1e10
        ymin = +1e10
        xmax = -1e10
        ymax = -1e10
        for c in raft:
            ccd   = cameraGeom.cast_Ccd(c)
            label = ccd.getId().getName()

            cc       = ccd.getCenter().getPixels(ccd.getPixelSize())
            cxc, cyc = cc.getX(), cc.getY()
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

            if cx0 < xmin:
                xmin = cx0
            if cx1 > xmax:
                xmax = cx1
            if cy0 < ymin:
                ymin = cy0
            if cy1 > ymax:
                ymax = cy1

            crectangle = Rectangle((cx0, cy0), cwidth, cheight, fill = False, label = label)
            rectangles[label] = crectangle
            centers[label] = (rxc+cxc, ryc+cyc)

            ccdBoundaries[label] = ((cx0, cx1), (cy0, cy1))
            

        raftBoundaries.append(((xmin, xmin), (ymin, ymax)))
        raftBoundaries.append(((xmin, xmax), (ymin, ymin)))
        raftBoundaries.append(((xmin, xmax), (ymax, ymax)))
        raftBoundaries.append(((xmax, xmax), (ymin, ymax)))

    return centers, rectangles, raftBoundaries, ccdBoundaries


