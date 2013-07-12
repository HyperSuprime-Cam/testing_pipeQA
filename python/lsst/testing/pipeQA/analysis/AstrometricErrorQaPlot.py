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
    dx = data['dx']
    dy = data['dy']
    limits = data['limits']
    gridVectors = data['gridVectors']


    xlo, xhi, ylo, yhi = limits
    xwid, ywid = xhi - xlo, yhi - ylo


    if not gridVectors:
        xbLo, xbHi, ybLo, ybHi = data['bbox']
        x -= xbLo
        y -= ybLo


    figsize = (6.5, 3.25)
    conv = colors.ColorConverter()
    black = conv.to_rgb('k')
    red = conv.to_rgb('r')
    green = conv.to_rgb('g')

    # if there's no data, dump a single point at 0,0
    if not len(x) > 0:
        x = numpy.array([0.0])
        y = numpy.array([0.0])
        dx = numpy.array([0.0])
        dy = numpy.array([0.0])


    if len(x) > 1:
        xlim = xlo, xhi
        ylim = ylo, yhi
    else:
        xlim = [0, 1.0]
        ylim = [0, 1.0]

    xmin, xmax = xlim
    ymin, ymax = ylim
        
    r = numpy.sqrt(dx**2 + dy**2)
    rmax = 0.6 # r.max()

    fig = figure.Figure(figsize=figsize)
    canvas = FigCanvas(fig)
    fig.subplots_adjust(left=0.1)

    ################
    # ccd view
    ax = fig.add_subplot(121)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + 0.1*box.height, box.width, 0.9*box.height])

    # for the 'all' plot, just show avg vectors in grid cells
    if gridVectors:
        nx, ny = 8, 8
        xstep, ystep = (xlim[1]-xlim[0])/nx, (ylim[1]-ylim[0])/ny
        if xstep == 0:
            xstep = 1.0;
        if ystep == 0:
            ystep = 1.0;

        xgrid = [[0.0]*nx for i in range(ny)]
        ygrid = [[0.0]*nx for i in range(ny)]
        ngrid = [[0]*nx for i in range(ny)]
        for i in range(len(x)):
            ix, iy = int((x[i]-xlim[0])/xstep), int((y[i]-ylim[0])/ystep)
            if ix >= 0 and ix < nx and iy >= 0 and iy < ny:
                xgrid[ix][iy] += dx[i]
                ygrid[ix][iy] += dy[i]
                ngrid[ix][iy] += 1

        xt, yt, dxt, dyt, nt = [],[],[],[],[]
        for ix in range(nx):
            for iy in range(ny):
                xt.append(xstep*(ix+0.5))
                yt.append(ystep*(iy+0.5))
                nt.append("%d" % (ngrid[ix][iy]))
                if ngrid[ix][iy] > 0:
                    xval = xgrid[ix][iy]/ngrid[ix][iy]
                    yval = ygrid[ix][iy]/ngrid[ix][iy]
                    dxt.append(xval)
                    dyt.append(yval)
                else:
                    dxt.append(0.0)
                    dyt.append(0.0)

        q = ax.quiver(numpy.array(xt), numpy.array(yt), numpy.array(dxt), numpy.array(dyt),
                      color='k', scale=1.0, angles='xy', pivot='middle', width=0.004)
        ax.quiverkey(q, 0.9, -0.2, 0.1, "100 mas", coordinates='axes',
                     fontproperties={'size':'small'})
        for i in range(len(xt)):
            ax.text(xt[i]+0.1*xstep, yt[i]+0.1*ystep, nt[i], size=5)
        xlim = [0.0, xlim[1] - xlim[0]]
        ylim = [0.0, ylim[1] - ylim[0]]
    else:
        ax.scatter(x, y, 0.5, color='r')
        q = ax.quiver(x, y, dx, dy, color='k', scale=10.0, angles='xy')
        ax.quiverkey(q, 0.9, -0.2, 1.0, "1 arcsec", coordinates='axes',
                     fontproperties={'size':"small"})

    ax.xaxis.set_major_locator(MaxNLocator(8))
    ax.yaxis.set_major_locator(MaxNLocator(8))

    ax.set_xlabel("x [pixels]")
    ax.set_ylabel("y [pixels]")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    for tic in ax.get_xticklabels() + ax.get_yticklabels():
        tic.set_size("x-small")


    ################
    # rose view
    xmargin = 0.07
    ymargin = 0.12
    spacer = 0.03
    left, bottom, width, height = 0.5+spacer, 0.35+spacer, 0.5-2*xmargin, 0.65-ymargin-spacer
    ax = fig.add_axes([left, bottom, width, height])

    #get the figure width/heigh in inches to correct aspect ratio
    f_w, f_h = fig.get_size_inches()
    xlimRose = [-width*f_w/(height*f_h)*rmax, f_w*width/(f_h*height)*rmax]
    limRose = [-rmax, rmax]
    # this is much faster than calling plot() in a loop, and quiver() scale length buggy
    if gridVectors:
        ybin = 50
        xbin = 50
        qaPlotUtil.make_densityplot(ax, dx, dy, xlims=limRose, ylims=limRose, bins=(xbin,ybin),
                                   log=True)
        qaPlotUtil.make_densityContour(ax, dx, dy, xlims=limRose, ylims=limRose, bins=(xbin,ybin),
                                      log=True, percentiles=True, normed=False, levels=[0.5])
        c0 = Circle((0.0, 0.0), radius=0.0, facecolor='none', edgecolor=green, zorder=3, label="50%")
        ax.vlines(0.0, limRose[0], limRose[1], linestyle='dashed')
        ax.hlines(0.0, xlimRose[0], xlimRose[1], linestyle='dashed')
        ax.add_patch(c0)
    else:
        z = numpy.zeros(len(dx))
        xy2 = zip(dx, dy)
        xy1 = zip(z, z)
        lines = zip(xy1, xy2)
        p = LineCollection(lines, colors=red*len(lines), zorder=1, label="_nolegend_")
        ax.add_collection(p)
        ax.scatter(dx, dy, s=0.05, color='k', zorder=2,label="_nolegend_")

        isort = r.argsort()
        i50 = isort[len(r)/2]
        r50 = r[i50]
        c50 = Circle((0.0, 0.0), radius=r50, facecolor='none', edgecolor=green, zorder=3, label="50%")
        ax.add_patch(c50)

        ax.vlines(0.0, limRose[0], limRose[1], linestyle='dashed')
        ax.hlines(0.0, xlimRose[0], xlimRose[1], linestyle='dashed')


    fp = fm.FontProperties(size="xx-small")
    ax.legend(prop=fp)

    ax.xaxis.set_label_position('top')
    ax.yaxis.set_label_position('right')
    ax.xaxis.set_ticks_position('top')
    ax.yaxis.set_ticks_position('right')

    ax.set_xlabel("dRa [arcsec]", size='x-small')
    ax.set_ylabel("dDec [arcsec]", size='x-small')
    ax.set_xlim(xlimRose)
    ax.set_ylim(limRose)
    for tic in ax.get_xticklabels() + ax.get_yticklabels():
        tic.set_size("x-small")


    ################
    # hist view
    left, bottom, width, height = 0.5+spacer, 0.0+ymargin, 0.5-2*xmargin, 0.35-ymargin
    ax0 = fig.add_axes([left, bottom, width, height])

    rNmax = 1.0
    if len(r) > 1:
        rsort = numpy.sort(r)
        nclip = int(0.9*len(r))
        rstd = numpy.std(rsort[0:nclip])
        binWidth = 0.5*rstd*(20.0/len(r))**0.2
        nBin = (rmax - r.min())/binWidth
        rN, rBin, xx = ax0.hist(r, bins=nBin)

        # add a median arrow
        rmed = numpy.median(r)
        w = numpy.where(rBin > rmed)[0]
        if len(w) > 0 and len(rN) > w[0]-1:
            histMed = 1.1*rN[w[0]-1]
        else:
            histMed = 0.0
        rNmax = rN.max()
        ax0.arrow(rmed, 1.2*rNmax, 0.0, histMed-1.2*rNmax, facecolor='r', edgecolor='r', lw=0.5)
        ax0.text(1.2*rmed, 0.5*(histMed+1.2*rNmax), "median",
                 verticalalignment="center", color='k', size='x-small')


    ax0.yaxis.set_ticks_position('right')
    ax0.yaxis.set_label_position('right')

    ax0.set_xlabel("r [arcsec]", size='x-small')
    ax0.set_ylabel("N", size='x-small')
    ax0.set_xlim([0, 0.9*rmax])
    ax0.set_ylim([0, 1.4*rNmax])

    for tic in ax0.get_xticklabels() + ax0.get_yticklabels(): # + ax.get_xticklabels():
        tic.set_size("xx-small")

    return fig



if __name__ == '__main__':
    filename, = sys.argv[1:2]
    data, isSummary = qaPlotUtil.unshelveGlob(filename, flag='r')
    if isSummary:
        data['gridVectors'] = True
    fig = plot(data)
    fig.savefig(filename)
