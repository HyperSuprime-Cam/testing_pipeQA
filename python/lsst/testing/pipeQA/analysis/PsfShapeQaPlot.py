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
import matplotlib.cm as cm


import QaPlotUtils as qaPlotUtil

def plot(data):

    t = data['t']
    x = data['x']
    y = data['y']
    dx = data['dx']
    dy = data['dy']
    color = data['color']
    limits = data['limits']
    vLen = data['vLen']
    summary = data['summary']
    vlim = data['vlim']
    fwhm = data['fwhm']


    if not summary:
        xbLo, xbHi, ybLo, ybHi = data['bbox']
        x -= xbLo
        y -= ybLo

    norm = colors.Normalize(vmin=vlim[0], vmax=vlim[1])
    sm = cm.ScalarMappable(norm, cmap=cm.jet)
    
    figsize = (5.0, 4.0)
    xlo, xhi, ylo, yhi = limits

    if len(x) == 0:
        x     = numpy.array([0.0])
        y     = numpy.array([0.0])
        dx    = numpy.array([0.0])
        dy    = numpy.array([0.0])
        color = numpy.array((0.0, 0.0, 0.0))
        fwhm  = numpy.array([0.0])

    xlim = [xlo, xhi]
    ylim = [ylo, yhi]

    fig = figure.Figure(figsize=figsize)
    canvas = FigCanvas(fig)
    
    fig.subplots_adjust(left=0.15, bottom=0.15)
    ax = fig.add_subplot(111)



    if summary:
        vmin, vmax = sm.get_clim()
        color = fwhm
        q = ax.quiver(x, y, vLen*dx, vLen*dy, color=sm.to_rgba(color),
                      scale=4.0*vLen, angles='xy', pivot='middle',
                      headlength=1.0, headwidth=1.0, width=0.002) 
        ax.quiverkey(q, 0.9, -0.12, 0.1*vLen, "e=0.1", coordinates='axes',
                     fontproperties={'size':"small"}, labelpos='E', color='k')
        q.set_array(color)
        cb = fig.colorbar(q)
        cb.ax.set_xlabel("FWHM$_{\mathrm{xc,yc}}$", size="small")
        cb.ax.xaxis.set_label_position('top')
        for tick in cb.ax.get_yticklabels():
            tick.set_size("x-small")
        ax.set_title("PSF Shape")
    else:
        q = ax.quiver(x, y, vLen*dx, vLen*dy, color=color, scale=4.0*vLen, angles='xy', pivot='middle',
                      headlength=1.0, headwidth=1.0, width=0.002)
        ax.quiverkey(q, 0.9, -0.12, 0.1*vLen, "e=0.1", coordinates='axes',
                     fontproperties={'size':"small"}, labelpos='E', color='k')
        ax.set_title("PSF Shape (FWHM$_{\mathrm{xc,yc}}$=%.2f)"%(fwhm[0]))

    ax.set_xlabel("x [pixels]")

    ax.set_ylabel("y [pixels]")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    for tic in ax.get_xticklabels() + ax.get_yticklabels():
        tic.set_size("x-small")
    for tic in ax.get_xticklabels():
        tic.set_rotation(22)
    for tic in ax.get_yticklabels():
        tic.set_rotation(45)

    return fig


if __name__ == '__main__':
    filename, = sys.argv[1:2]
    data, isSummary = qaPlotUtil.unshelveGlob(filename, flag='r')
    if isSummary:
        data['summary'] = True
        data['limits'] = data['alllimits']
        data['vLen'] = 5.0*data['vLen']
    fig = plot(data)
    fig.savefig(filename)
