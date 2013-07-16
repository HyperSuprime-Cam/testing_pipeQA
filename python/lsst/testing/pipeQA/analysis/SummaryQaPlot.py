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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import sys
import numpy as num
import matplotlib

import matplotlib.figure as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
from matplotlib import colors
import matplotlib.font_manager as fm
from  matplotlib.ticker import MaxNLocator
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle


import QaPlotUtils as qaPlotUtil

def plot(data):

    x = data['x']
    y = data['y']
    
    summary = data['summary']
    
    if len(x) == 0:
        x = num.array([0.0])
        y = num.array([0.0])

    figsize = (4.0, 4.0)
    fig = figure.Figure(figsize=figsize)
    canvas = FigCanvas(fig)
        
    if not summary:
        ax = fig.add_subplot(111)
        ax.plot(x, y)
    else:
        ax = fig.add_subplot(111)
        ax.plot(x, y)

    return fig



if __name__ == '__main__':
    filename, = sys.argv[1:2]
    data, isSummary = qaPlotUtil.unshelveGlob(filename, flag='r')
    if isSummary:
        data['summary'] = True
    fig = plot(data)
    fig.savefig(filename)
