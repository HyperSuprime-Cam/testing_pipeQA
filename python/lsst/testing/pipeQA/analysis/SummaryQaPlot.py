#!/usr/bin/env python

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
