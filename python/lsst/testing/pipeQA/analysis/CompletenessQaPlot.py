#!/usr/bin/env python


import sys
import shelve

import numpy as num
import matplotlib

import matplotlib.figure as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
from matplotlib import colors
import matplotlib.font_manager as fm
from  matplotlib.ticker import MaxNLocator
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle

import matplotlib.ticker as ticker
from matplotlib.font_manager import FontProperties
import matplotlib.patches as patches

import QaPlotUtils as qaPlotUtil


def plot(data):

    title            = data['title']
    orphan           = data['orphan']
    depth            = data['depth']
    matchedStar      = data['matchedStar']
    blendedStar      = data['blendedStar']
    undetectedStar   = data['undetectedStar']
    matchedGalaxy    = data['matchedGalaxy']
    blendedGalaxy    = data['blendedGalaxy']
    undetectedGalaxy = data['undetectedGalaxy']
    bins             = data['bins']

    figsize = (4.0, 4.0)
    fig = figure.Figure(figsize=figsize)
    canvas = FigCanvas(fig)
    
    sp1 = fig.add_subplot(211)
    fig.subplots_adjust(left=0.13)
    sp2 = fig.add_subplot(212, sharex = sp1)

    # Stacked histogram
    orphanHist         = num.histogram(orphan, bins=bins)
    matchedStarHist    = num.histogram(matchedStar, bins=bins)
    blendedStarHist    = num.histogram(blendedStar, bins=bins)
    undetectedStarHist = num.histogram(undetectedStar, bins=bins)
    # For bar, you send the coordinate of the left corner of the bar
    barbins    = orphanHist[1][:-1]
    width      = 1.0 * (orphanHist[1][1] - orphanHist[1][0])
    orphanBar  = sp1.bar(barbins, orphanHist[0], width=width, color='r', alpha = 0.5, label = 'Orphan', capsize = 1)
    bottom     = orphanHist[0]
    matchedBar = sp1.bar(barbins, matchedStarHist[0], width=width, color='g', alpha=0.5, label='Matched',
                         bottom=bottom, capsize=1)
    bottom    += matchedStarHist[0]
    blendedBar = sp1.bar(barbins, blendedStarHist[0], width=width, color='cyan', alpha=0.5, label='Blended',
                         bottom=bottom, capsize=1)
    bottom    += blendedStarHist[0]
    unmatBar   = sp1.bar(barbins, undetectedStarHist[0], width=width, color='b', alpha=0.5, label='Unmatched',
                         bottom=bottom, capsize=1)

    ymax = num.max(orphanHist[0] + matchedStarHist[0] + blendedStarHist[0] + undetectedStarHist[0])
    sp1.set_ylim([0, 1.4*ymax])

    sp1x2           = sp1.twinx()
    allStars        = num.concatenate((matchedStar, blendedStar, undetectedStar))
    foundStars      = num.concatenate((matchedStar, blendedStar))
    histAll         = num.histogram(allStars, bins=bins)
    histFound       = num.histogram(foundStars, bins=bins)

    magbins = 0.5 * (histAll[1][1:] + histAll[1][:-1])
    w       = num.where(histAll[0] != 0)
    x       = magbins[w]
    n       = 1.0 * histFound[0][w]
    d       = 1.0 * histAll[0][w]
    y       = n / d  
    sp1x2.plot(x, y)
    sp1x2.set_ylim([0.0, 1.4])
    sp1x2.set_ylabel('(Match+Blend)/Tot', fontsize=8)
    sp1x2.axhline(y = 0.5, c='k', linestyle='-.', alpha = 0.75)
    sp1x2.axhline(y = 1.0, c='k', linestyle='-.', alpha = 0.75)
    sp1x2.axvline(x = depth, c='k', linestyle='-', alpha = 0.75)
    sp1x2.text(depth-1.0, 1.2, "%.2f" % (depth), fontsize=8,
               horizontalalignment='right', verticalalignment='center')
    fa = patches.FancyArrow(depth-0.8, 1.2, 0.8, 0.0, length_includes_head=True,
                            overhang=0.2, head_length=0.2, head_width=0.04)
    sp1x2.add_patch(fa)

    qaPlotUtil.qaSetp(sp1x2.get_xticklabels(), visible=False)
    qaPlotUtil.qaSetp(sp1x2.get_yticklabels(), fontsize = 6)

    sp1.set_ylabel('N Stars', fontsize=10)
    qaPlotUtil.qaSetp(sp1.get_xticklabels()+sp1.get_yticklabels(), fontsize = 6)

    ##############

    orphanHist        = num.histogram(orphan, bins=bins)
    matchedGalHist    = num.histogram(matchedGalaxy, bins=bins)
    blendedGalHist    = num.histogram(blendedGalaxy, bins=bins)
    undetectedGalHist = num.histogram(undetectedGalaxy, bins=bins)
    orphanBar  = sp2.bar(barbins, orphanHist[0], width=width, color='r', alpha = 0.5, label = 'Orphan', capsize = 1, log=False)
    bottom     = orphanHist[0]
    matchedBar = sp2.bar(barbins, matchedGalHist[0], width=width, color='g', alpha=0.5, label='Matched',
                         bottom=bottom, capsize=1, log=False)
    bottom    += matchedGalHist[0]
    blendedBar = sp2.bar(barbins, blendedGalHist[0], width=width, color='cyan', alpha=0.5, label='Blended',
                         bottom=bottom, capsize=1, log=False)
    bottom    += blendedGalHist[0]
    unmatBar   = sp2.bar(barbins, undetectedGalHist[0], width=width, color='b', alpha=0.5, label='Unmatched',
                         bottom=bottom, capsize=1, log=False)

    sp2.set_xlabel('Mag', fontsize=10)
    sp2.set_ylabel('N Gals', fontsize=10)
    qaPlotUtil.qaSetp(sp2.get_xticklabels()+sp2.get_yticklabels(), fontsize = 6)
    qaPlotUtil.qaSetp(sp2.get_yticklabels(), rotation = 45.0)
    sp2.legend(numpoints = 1, prop=FontProperties(size='x-small'), loc = 'upper left')

    sp1.set_xlim(14, 26)

    fig.suptitle('%s Stacked histogram' % (title), fontsize = 11)

    return fig



if __name__ == '__main__':
    filename, = sys.argv[1:2]
    data, isSummary = qaPlotUtil.unshelveGlob(filename, flag='r')
    if isSummary:
        data['title'] = 'All Sensors'
    fig = plot(data)
    fig.savefig(filename)

