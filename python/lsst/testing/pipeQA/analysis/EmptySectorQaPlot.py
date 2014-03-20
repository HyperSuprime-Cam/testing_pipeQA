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

    ra = data['ra']
    dec = data['dec']
    ramat = data['ramat']
    decmat = data['decmat']
    w = numpy.where(numpy.isfinite(ra) & numpy.isfinite(dec))
    wm = numpy.where(numpy.isfinite(ramat) & numpy.isfinite(decmat))
    ra_lo = min(ra[w].min(), ramat[wm].min())
    ra_hi = max(ra[w].max(), ramat[wm].max())
    dec_lo = min(dec[w].min(), decmat[wm].min())
    dec_hi = max(dec[w].max(), decmat[wm].max())

    buff = 0.1
    if summary:
        ra_lo -= buff
        ra_hi += buff
        dec_lo -= buff
        dec_hi += buff
        
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

    figsize = (7.0, 3.5)

    ####################
    # create the plot
    fig = figure.Figure(figsize=figsize)
    canvas = FigCanvas(fig)
    ax = fig.add_subplot(121)
    rax = fig.add_subplot(122)
    fig.subplots_adjust(bottom=0.15) #left=0.19) #, bottom=0.15)
    
    ncol = None
    if summary:
        ms = 0.5
        if len(xmat) < 10000:
            ms = 1.0
        if len(xmat) < 1000:
            ms = 2.0
        ax.plot(xmat, ymat, "k.", ms=ms, label=summaryLabel)
        rax.plot(ramat, decmat, "k.", ms=ms, label=summaryLabel)
        ncol = 1
    else:
        ax.plot(x, y, "k.", ms=2.0, label="detected")
        ax.plot(xmat, ymat, "ro", ms=4.0, label="matched",
                mfc='None', markeredgecolor='r')
        rax.plot(ra, dec, "k.", ms=2.0, label="detected")
        rax.plot(ramat, decmat, "ro", ms=4.0, label="matched",
                mfc='None', markeredgecolor='r')
        ncol = 2


    ax.set_xlim([xlo, xhi])
    ax.set_ylim([ylo, yhi])
    ax.set_xlabel("x [pixel]", size='x-small')
    ax.set_ylabel("y [pixel]", size='x-small')
    #ax.legend(prop=fm.FontProperties(size ="xx-small"), ncol=ncol, loc="upper center")
    for tic in ax.get_xticklabels() + ax.get_yticklabels():
        tic.set_size("x-small")

    rax_twin = rax.twinx() 
    rax.get_xaxis().get_major_formatter().set_useOffset(False)
    rax_twin.get_yaxis().get_major_formatter().set_useOffset(False)
    rax.set_xlim([ra_lo, ra_hi])
    rax.set_ylim([dec_lo, dec_hi])
    rax_twin.set_ylim([dec_lo, dec_hi])
    rax.get_yaxis().set_visible(False)
    rax.set_xlabel("R.A. [deg]", size='x-small')
    rax_twin.set_ylabel("Decl. [deg]", size='x-small')
    #rax.legend(prop=fm.FontProperties(size ="xx-small"), ncol=ncol, loc="upper center")
    for tic in rax.get_xticklabels() + rax_twin.get_yticklabels():
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
        if False:
            dx, dy = 20, 20  # on a 4kx4k ccd, < +/-20 pixels is tough to hit with a mouse
            for i in range(len(x)):
                area = x[i]-dx, y[i]-dy, x[i]+dx, y[i]+dy
                fig.addMapArea("no_label_info", area, "nolink:%.1f_%.1f"%(x[i],y[i]))

                
        fig.suptitle("Matched Detections by CCD Sector", size='small')
    else:
        fig.suptitle("Matched Detections", size='small')


    return fig


if __name__ == '__main__':
    filename, = sys.argv[1:2]
    data, isSummary = qaPlotUtil.unshelveGlob(filename, flag='r')
    if isSummary:
        data['summary'] = True
        data['limits'] = data['alllimits']
    fig = plot(data)
    fig.savefig(filename)
