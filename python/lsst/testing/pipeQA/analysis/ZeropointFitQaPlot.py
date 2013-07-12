#!/usr/bin/env python

import sys
import numpy  as num
import matplotlib

import matplotlib.figure as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
from matplotlib import colors
import matplotlib.font_manager as fm
from  matplotlib.ticker import MaxNLocator
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle

from matplotlib.font_manager import FontProperties
from matplotlib.patches import Ellipse

import QaPlotUtils as qaPlotUtil

def plot(data):


    mrefGmag  = data['mrefGmag']
    mimgGmag  = data['mimgGmag']
    mimgGmerr = data['mimgGmerr']
    mrefSmag  = data['mrefSmag']
    mimgSmag  = data['mimgSmag']
    mimgSmerr = data['mimgSmerr']
    urefmag   = data['urefmag']
    uimgmag   = data['uimgmag']
    zeropt    = data['zeropt']
    title     = data['title']
    figsize   = data['figsize']
    fluxType  = data['fluxType']
    
    # Just to get the histogram results

    fig = figure.Figure(figsize=figsize)
    canvas = FigCanvas(fig)

    legLines  = []
    legLabels = []

    axis = fig.add_axes([0.225, 0.225, 0.675, 0.550])


    ################    --------------------------------  ##################
    # can't fig.savefig() with this Ellipse stuff ??? ... i love matplotlib
    if False:
        for i in range(len(mrefGmag)):
            a = Ellipse(xy=num.array([mimgGmag[i], mrefGmag[i]]),
                        width=mimgGmerr[i]/2., height=mimgGmerr[i]/2.,
                        alpha=0.5, fill=True, ec='g', fc='g', zorder=10)
            axis.add_artist(a)
    else:
        pass
        #axis.plot(mimgGmag, mrefGmag, 'g.', zorder=10, alpha=0.5)
    #########################################################################


    mimgGplot = axis.plot(mimgGmag, mrefGmag, '.', color='g', mfc='g', mec='g',
                          zorder=10, label = 'Matched Galaxies', ms=2.5)

    legLines.append(mimgGplot[0])
    legLabels.append("Matched Galaxies")


    ################    --------------------------------  ##################
    # can't fig.savefig() with this Ellipse stuff ??? ... i love matplotlib
    if False:
        for i in range(len(mrefSmag)):
            a = Ellipse(xy=num.array([mimgSmag[i], mrefSmag[i]]),
                        width=mimgSmerr[i]/2., height=mimgSmerr[i]/2.,
                        alpha=0.5, fill=True, ec='b', fc='b', zorder=12)
            axis.add_artist(a)
    else:
        pass
        #axis.plot(mimgSmag, mrefSmag, "b.", zorder=12, alpha=0.5)
    ########################################################################

    mimgSplot = axis.plot(mimgSmag, mrefSmag, '.', color='b', mfc='b', mec='b',
                          zorder=12, label = 'Matched Stars', ms=2.5)

    legLines.append(mimgSplot[0])
    legLabels.append("Matched Stars")

    if len(mrefGmag) == 0 and len(mrefSmag) == 0:
        xmin, xmax, ymin, ymax = -15, -8, 16, 28
    else:
        xmin, xmax, ymin, ymax = axis.axis()

    # Plot zpt
    xzpt = num.array((xmin, xmax))
    pzpt = axis.plot(xzpt, xzpt - zeropt, 'k--', label = 'Zeropoint')
    legLines.append(pzpt)
    legLabels.append("Zeropoint")

    maxN2 = 999

    nu, bu, pu = None, None, None

    # Unmatched & matched reference objects
    ax2  = fig.add_axes([0.1,   0.225, 0.125, 0.550], sharey=axis)
    if len(urefmag) > 0:
        nu, bu, pu = ax2.hist(urefmag, bins=num.arange(ymin, ymax, 0.25),
                              orientation='horizontal', facecolor = 'r',
                              log = True, alpha = 0.5, zorder = 1)
        if num.max(nu) > maxN2: maxN2 = num.max(nu)
        legLines.append(pu[0])
        legLabels.append("Unmatched Sources")

    if len(mrefGmag) > 0 and len(mrefSmag) > 0:
        nu, bu, pu = ax2.hist(num.concatenate((mrefGmag,mrefSmag)), bins=num.arange(ymin, ymax, 0.25),
                 orientation='horizontal', log = True, facecolor = 'b', alpha = 0.5, zorder = 2)
    elif len(mrefGmag):
        nu, bu, pu = ax2.hist(mrefGmag, bins=num.arange(ymin, ymax, 0.25),
                 orientation='horizontal', log = True, facecolor = 'b', alpha = 0.5, zorder = 2)
    elif len(mrefSmag) > 0:
        nu, bu, pu =  ax2.hist(mrefSmag, bins=num.arange(ymin, ymax, 0.25),
                 orientation='horizontal', log = True, facecolor = 'b', alpha = 0.5, zorder = 2)
    ax2.set_xlabel('N', fontsize = 10)
    ax2.set_ylabel('Reference catalog mag', fontsize = 10)


    if nu is not None and num.max(nu) > maxN2: maxN2 = num.max(nu)

    maxN3 = 999

    nm, bm, pm = None, None, None

    # Unmatched & matched stellar objects
    ax3  = fig.add_axes([0.225, 0.1,   0.675, 0.125], sharex=axis)
    ax3.get_yaxis().set_ticks_position('right')
    ax3.get_yaxis().set_label_position('right')
    if len(mimgGmag) > 0 and len(mimgSmag) > 0:
        nm, bm, pm = ax3.hist(num.concatenate((mimgGmag,mimgSmag)), bins=num.arange(xmin, xmax, 0.25),
                              log = True, facecolor = 'b', alpha = 0.5, zorder = 2)
        legLines.append(pm[0])
        legLabels.append("Matched Sources")
    elif len(mimgGmag) > 0:
        nm, bm, pm = ax3.hist(mimgGmag, bins=num.arange(xmin, xmax, 0.25),
                              log = True, facecolor = 'b', alpha = 0.5, zorder = 2)
        legLines.append(pm[0])
        legLabels.append("Matched Sources")
    elif len(mimgSmag) > 0:
        nm, bm, pm = ax3.hist(mimgSmag, bins=num.arange(xmin, xmax, 0.25),
                              log = True, facecolor = 'b', alpha = 0.5, zorder = 2)
        legLines.append(pm[0])
        legLabels.append("Matched Sources")

    if nm is not None and num.max(nm) > maxN3: maxN3 = num.max(nm)

    if len(uimgmag) > 0:
        try:
            ax3.hist(uimgmag, bins=num.arange(xmin, xmax, 0.25),
                     log = True, facecolor = 'r', alpha = 0.5, zorder = 1)
        except:
            pass
    if abs(zeropt) < 1.0e-5:
        ax3.set_xlabel('Image calibrated %s mag' % (fluxType), fontsize = 10)
    else:
        ax3.set_xlabel('Image instrumental %s mag' % (fluxType), fontsize = 10)
    ax3.set_ylabel('N', rotation = 180, fontsize = 10)

    # Mag - Zpt
    ax4  = fig.add_axes([0.225, 0.775, 0.675, 0.125], sharex=axis)
    if len(mimgSmag):
        mimgSeb = ax4.errorbar(mimgSmag, mimgSmag - zeropt - mrefSmag,
                               yerr = mimgSmerr,
                               fmt = 'bo', ms = 2, alpha = 0.25, capsize = 0, elinewidth = 0.5)
        mimgSeb[2][0].set_alpha(0.25) # alpha for error bars

    if len(mimgGmag):
        mimgGeb = ax4.errorbar(mimgGmag, mimgGmag - zeropt - mrefGmag,
                               yerr = mimgGmerr,
                               fmt = 'go', ms = 2, alpha = 0.25, capsize = 0, elinewidth = 0.5)
        mimgGeb[2][0].set_alpha(0.25) # alpha for error bars

    ax4.get_yaxis().set_ticks_position('right')
    ax4.get_yaxis().set_label_position('right')
    ax4.set_ylabel('Cal-Ref', fontsize = 10, rotation = 270)
    ax4.axhline(y = 0, c='k', linestyle='--', alpha = 0.25)

    # Cleaning up figure
    qaPlotUtil.qaSetp(axis.get_xticklabels()+axis.get_yticklabels(), visible=False)
    qaPlotUtil.qaSetp(ax2.get_xticklabels()+ax2.get_yticklabels(), fontsize = 8)
    qaPlotUtil.qaSetp(ax2.get_xticklabels(), rotation=90)
    qaPlotUtil.qaSetp(ax3.get_xticklabels()+ax3.get_yticklabels(), fontsize = 8)
    qaPlotUtil.qaSetp(ax4.get_xticklabels(), visible=False)
    qaPlotUtil.qaSetp(ax4.get_yticklabels(), fontsize = 8)

    fig.legend(legLines, legLabels,
                   numpoints=1, prop=FontProperties(size='x-small'), loc = 'center right')

    fig.suptitle('%s' % (title), fontsize = 12)

    numerator   = (mimgSmag - zeropt - mrefSmag)
    if len(numerator) == 0:
        numerator = num.array([0.0])
    denominator = mimgSmerr
    med         = num.median(numerator) 
    ax4.axhline(y = med, c='k', linestyle=':', alpha = 0.5)

    # Final axis limits
    ax2.set_xlim(0.75, 3.0*maxN2)
    ax3.set_ylim(0.75, 3.0*maxN3)
    ax4.set_ylim(-0.24, 0.24)

    #axis.axis((xmax, xmin, ymax, ymin))
    axis.set_ylim(26, 14)
    axis.set_xlim(26 + zeropt, 14 + zeropt)

    return fig



if __name__ == '__main__':
    filename, = sys.argv[1:2]
    data, isSummary = qaPlotUtil.unshelveGlob(filename, flag='r')
    if isSummary:
        data['summary'] = True
    fig = plot(data)
    fig.savefig(filename)
