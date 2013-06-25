#!/usr/bin/env python

import sys
import numpy
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
    xmat = data['xmat']
    ymat = data['ymat']
    limits = data['limits']
    summary = data['summary']
    nx, ny  = data['nxn']
    
    xlo, xhi, ylo, yhi = limits
    xwid, ywid = xhi - xlo, yhi - ylo

    if not summary:
        xbLo, xbHi, ybLo, ybHi = data['bbox']
        x -= xbLo
        y -= ybLo
        xmat -= xbLo
        ymat -= ybLo
        
        
    
    # handle no-data possibility
    if len(x) == 0:
        x = numpy.array([0.0])
        y = numpy.array([0.0])

    summaryLabel = "matched"
    if len(xmat) == 0:
        if len(x) == 0 or not summary:
            xmat = numpy.array([0.0])
            ymat = numpy.array([0.0])
        else:
            summaryLabel = "detected"
            xmat = x
            ymat = y

    figsize = (4.0, 4.0)

    ####################
    # create the plot
    fig = figure.Figure(figsize=figsize)
    canvas = FigCanvas(fig)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(left=0.19) #, bottom=0.15)

    ncol = None
    if summary:
        ms = 0.1
        if len(xmat) < 10000:
            ms = 0.5
        if len(xmat) < 1000:
            ms = 1.0
        ax.plot(xmat, ymat, "k.", ms=ms, label=summaryLabel)
        ncol = 1
    else:
        ax.plot(x, y, "k.", ms=2.0, label="detected")
        ax.plot(xmat, ymat, "ro", ms=4.0, label="matched",
                mfc='None', markeredgecolor='r')
        ncol = 2


    ax.set_xlim([xlo, xhi])
    ax.set_ylim([ylo, yhi])
    ax.set_xlabel("x [pixel]", size='x-small')
    ax.set_ylabel("y [pixel]", size='x-small')
    ax.legend(prop=fm.FontProperties(size ="xx-small"), ncol=ncol, loc="upper center")
    for tic in ax.get_xticklabels() + ax.get_yticklabels():
        tic.set_size("x-small")

    # don't bother with this stuff for the final summary plot
    if not summary:
        # show the regions
        for i in range(nx):
            xline = (i+1)*xwid/nx
            ax.axvline(xline, color="k")
        for i in range(ny):
            yline = (i+1)*ywid/ny
            ax.axhline(yline, color="k")

        # add map areas to allow mouseover tooltip showing pixel coords
        dx, dy = 20, 20  # on a 4kx4k ccd, < +/-20 pixels is tough to hit with a mouse
        #for i in range(len(x)):
        #    area = x[i]-dx, y[i]-dy, x[i]+dx, y[i]+dy
        #    fig.addMapArea("no_label_info", area, "nolink:%.1f_%.1f"%(x[i],y[i]))

        ax.set_title("Matched Detections by CCD Sector", size='small')
    else:
        ax.set_title("Matched Detections", size='small')


    return fig


if __name__ == '__main__':
    filename, = sys.argv[1:2]
    data, isSummary = qaPlotUtil.unshelveGlob(filename, flag='r')
    if isSummary:
        data['summary'] = True
        data['limits'] = data['alllimits']
    fig = plot(data)
    fig.savefig(filename)
