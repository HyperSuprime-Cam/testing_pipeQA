#!/usr/bin/env python

import sys, copy
import numpy
import matplotlib

import matplotlib.figure             as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
from matplotlib                      import colors
import matplotlib.font_manager       as fm
from  matplotlib.ticker              import MaxNLocator
from matplotlib.collections          import LineCollection
from matplotlib.patches              import Circle


import QaPlotUtils                   as qaPlotUtil
import QaAnalysisUtils               as qaAnaUtil


def plot(data):

    mag0      = data['mag0']
    diff0     = data['diff0']
    star0     = data['star0']
    derr0     = data['derr0']
    areaLabel = data['areaLabel']
    raft      = data['raft']
    ccd       = data['ccd']
    figsize   = data['figsize']
    xlim      = data['xlim']
    ylim      = data['ylim']
    xlim2     = data['xlim2']
    ylim2     = data['ylim2']
    ylimStep  = data['ylimStep']
    tag1      = data['tag1']
    tag       = data['tag']
    mode      = data['mode']

    x         = data['x']
    y         = data['y']
    trend     = data['trend']
    magCut    = data['magCut']

    figType       = data['figType']

    if figType == 'standard':
        return standardFigure(mag0, diff0, star0, derr0, areaLabel, raft, ccd,
                              figsize, xlim, ylim, xlim2, ylim2, ylimStep,
                              tag1, tag, mode, x, y, trend, magCut)
    if figType == 'derr':
        return derrFigure(mag0, diff0, star0, derr0, areaLabel,
                          raft, ccd, figsize, xlim, ylim, xlim2, ylim2,
                          ylimStep, tag1, tag, mode, x, y, trend, magCut)

    if figType == 'summary':
        return summaryFigure(mag0, diff0, star0, derr0, areaLabel,
                             raft, ccd, figsize, xlim, ylim, xlim2, ylim2,
                             ylimStep, tag1, tag, mode, x, y, trend, magCut)

    
def standardFigure(*args):
    
    mag0, diff0, star0, derr0, areaLabel, raft, ccd, figsize, xlim, ylim, xlim2, ylim2, ylimStep, \
        tag1, tag, mode, x, y, trend, magCut = args
    
    eps = 1.0e-5

    conv = colors.ColorConverter()
    size = 2.0

    xlimDefault = [14.0, 25.0]

    red = conv.to_rgba('r')
    black = conv.to_rgba('k')

    x0    = x.copy()
    y0    = y.copy()
    modeIdx = {'all': 0, 'galaxies': 1, 'stars': 2 }
    lineFit = copy.copy(trend)[modeIdx[mode]]

    trendCoeffs = lineFit[0], lineFit[2]
    trendCoeffsLo = lineFit[0]+lineFit[1], lineFit[2]-lineFit[3]
    trendCoeffsHi = lineFit[0]-lineFit[1], lineFit[2]+lineFit[3]
    #print trendCoeffs        
    if len(mag0) == 0:
        mag0 = numpy.array([xlim[1]])
        diff0 = numpy.array([eps])
        x0    = numpy.array([eps])
        y0    = numpy.array([eps])
        star0 = numpy.array([0])

    #################
    # data for one ccd
    if mode == 'fourPanel':
        figsize = (6.5, 5.0)


    fig = figure.Figure(figsize=figsize)
    canvas = FigCanvas(fig)
    fig.subplots_adjust(left=0.09, right=0.93, bottom=0.125)

    if mode == 'fourPanel':
        ax_1s = fig.add_subplot(221)
        ax_2s = fig.add_subplot(222)
        ax_1g = fig.add_subplot(223)
        ax_2g = fig.add_subplot(224)
        axSets        = [ [ax_1s, ax_2s], [ax_1g, ax_2g] ]
        starGalLabels = ["stars", "galaxies"]
        whereStarGals = [numpy.where(star0 == 0)[0], numpy.where(star0 > 0)[0] ]
    else:
        ax_1 = fig.add_subplot(131)
        ax_2 = fig.add_subplot(132)
        ax_3 = fig.add_subplot(133)
        axSets        = [ [ax_1, ax_2, ax_3] ]
        starGalLabels = [mode]
        if mode == 'stars':
            whereStarGals = [numpy.where(star0 > 0)[0] ]
        if mode == 'galaxies':
            whereStarGals = [numpy.where(star0 == 0)[0] ]
        if mode == 'all':
            starGalLabels = ["all data"]
            whereStarGals = [numpy.where(star0 > -1)[0] ]


    for iSet in range(len(axSets)):
        ax_1, ax_2, ax_3   = axSets[iSet]
        starGalLabel = starGalLabels[iSet]
        whereStarGal = whereStarGals[iSet]

        mag  = mag0[whereStarGal]
        diff = diff0[whereStarGal]
        x    = x0[whereStarGal]
        y    = y0[whereStarGal]
        star = star0[whereStarGal]

        if len(x) == 0:
            mag = numpy.array([eps])
            diff = numpy.array([eps])
            x = numpy.array([eps])
            y = numpy.array([eps])
            star = numpy.array([0])


        whereCut = numpy.where((mag < magCut))[0]
        whereOther = numpy.where((mag > magCut))[0]

        xTrend = numpy.array(xlim)
        ax_1.plot(xTrend, numpy.array([eps, eps]), "-k", lw=1.0)

        ax_1.text(1.02*xlim[0], 0.87*ylim[1], starGalLabel, size='x-small', horizontalalignment='left')
        for ax in [ax_1, ax_2]:
            ax.plot(mag[whereOther], diff[whereOther], "k.", ms=size, label=ccd)
            ax.plot(mag[whereCut], diff[whereCut], "r.", ms=size, label=ccd)
            ax.set_xlabel(tag1, size='small')

        norm = colors.Normalize(vmin=-0.05, vmax=0.05, clip=True)
        cdict = {'red': ((0.0, 0.0, 0.0),
                         (0.5, 0.0, 0.0),
                         (1.0, 1.0, 1.0)),
                 'green': ((0.0, 0.0, 0.0),
                           (0.5, 0.0, 0.0),
                           (1.0, 0.0, 0.0)),
                 'blue': ((0.0, 1.0, 1.0),
                          (0.5, 0.0, 0.0),
                          (1.0, 0.0, 0.0))}
        my_cmap = colors.LinearSegmentedColormap('my_colormap',cdict,256)

        midSize = 2.0
        maxSize = 6*midSize
        minSize = midSize/2
        magLo = xlim[0] if xlim[0] > mag0.min() else mag0.min()
        magHi = magCut + 2.0 #xlim[1] if xlim[1] < mag0.max() else mag0.max()
        sizes = minSize + (maxSize - minSize)*(magHi - mag)/(magHi - magLo)
        sizes = numpy.clip(sizes, minSize, maxSize)
        xyplot = ax_3.scatter(x, y, s=sizes, c=diff, marker='o',
                              cmap=my_cmap, norm=norm, edgecolors='none')

        if len(x0) > 1:
            ax_3.set_xlim([x0.min(), x0.max()])
            ax_3.set_ylim([y0.min(), y0.max()])
        else:
            ax_3.set_xlim(xlimDefault)
            ax_3.set_ylim(xlimDefault)

        ax_3.set_xlabel("x", size='x-small')
        #ax_3.set_ylabel("y", labelpad=20, y=1.0, size='x-small', rotation=0.0)
        cb = fig.colorbar(xyplot)
        cb.ax.set_xlabel("$\delta$ mag", size="x-small")
        cb.ax.xaxis.set_label_position('top')
        for tick in cb.ax.get_yticklabels():
            tick.set_size("x-small")
        for t in ax_3.get_xticklabels() + ax_3.get_yticklabels():
            t.set_rotation(45.0)
            t.set_size('xx-small')
        #for t in ax_3.get_yticklabels():
        #    t.set_size('xx-small')

        lineVals = numpy.lib.polyval(trendCoeffs, xTrend)
        lineValsLo = numpy.lib.polyval(trendCoeffsLo, xTrend)
        lineValsHi = numpy.lib.polyval(trendCoeffsHi, xTrend)


        if abs(trendCoeffs[0] - 99.0) > 1.0e-6:
            lmin, lmax = lineVals.min(), lineVals.max()
            if False:
                if lmin < ylim[0]:
                    ylim[0] = -(int(abs(lmin)/ylimStep) + 1)*ylimStep
                if lmax > ylim[1]:
                    ylim[1] = (int(lmax/ylimStep) + 1)*ylimStep

            ax_1.plot(xTrend, lineVals, "-r", lw=1.0)
            ax_1.plot(xTrend, lineValsLo, "--r", lw=1.0)
            ax_1.plot(xTrend, lineValsHi, "--r", lw=1.0)
            ax_2.plot(xTrend, lineVals, "-r", lw=1.0)
            ax_2.plot(xTrend, lineValsLo, "--r", lw=1.0)
            ax_2.plot(xTrend, lineValsHi, "--r", lw=1.0)

        ax_2.plot([xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]],
                  [ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]], '-k')
        ax_1.set_ylabel(tag, size="small")
        ax_1.set_xlim(xlim if xlim[0] != xlim[1] else xlimDefault)
        ax_2.set_xlim(xlim2 if xlim2[0] != xlim2[1] else xlimDefault)
        ax_1.set_ylim(ylim if ylim[0] != ylim[1] else [-0.1, 0.1])
        ax_2.set_ylim(ylim2 if ylim2[0] != ylim2[1] else [-0.1, 0.1])

        # move the y axis on right panel
        #ax_3.yaxis.set_label_position('right')
        #ax_3.yaxis.set_ticks_position('right')

        dmag = 0.1
        ddiff1 = 0.02
        ddiff2 = ddiff1*(ylim2[1]-ylim2[0])/(ylim[1]-ylim[0]) # rescale for larger y range

        if False:
            for j in xrange(len(mag)):
                info = "nolink:x:%.2f_y:%.2f" % (x[j], y[j])
                area = (mag[j]-dmag, diff[j]-ddiff1, mag[j]+dmag, diff[j]+ddiff1)
                fig.addMapArea(areaLabel, area, info, axes=ax_1)
                area = (mag[j]-dmag, diff[j]-ddiff2, mag[j]+dmag, diff[j]+ddiff2)
                fig.addMapArea(areaLabel, area, info, axes=ax_2)

        for ax in [ax_1, ax_2]:
            for tick in ax.get_xticklabels() + ax.get_yticklabels():
                tick.set_size("x-small")

    return fig




def derrFigure(*args):
    mag0, diff0, star0, derr0, areaLabel, raft, ccd, figsize, xlim, ylim, xlim2, ylim2, ylimStep, \
          tag1, tag, mode, x, y, trend, magCut = args

    eps = 1.0e-5

    conv  = colors.ColorConverter()
    size  = 2.0
    red   = conv.to_rgba('r')
    black = conv.to_rgba('k')
    mode = "stars"  # this better be the case!        
    if len(mag0) == 0:
        mag0 = numpy.array([xlim[1]])
        diff0 = numpy.array([eps])
        derr0 = numpy.array([eps])
        star0 = numpy.array([0])

    fig = figure.Figure(figsize=figsize)
    canvas = FigCanvas(fig)
    fig.subplots_adjust(left=0.09, right=0.93, bottom=0.125)

    if len(mag0) == 0:
        mag0 = numpy.array([xlim[1]])
        diff0 = numpy.array([eps])
        derr0 = numpy.array([eps])
        star0 = numpy.array([0])

    whereStarGal = numpy.where(star0 > 0)[0]
    mag  = mag0[whereStarGal]
    diff = diff0[whereStarGal]
    derr = derr0[whereStarGal]
    star = star0[whereStarGal]

    if len(mag) == 0:
        mag = numpy.array([eps])
        diff = numpy.array([eps])
        derr = numpy.array([eps])
        star = numpy.array([0])

    whereCut = numpy.where((mag < magCut))[0]
    whereOther = numpy.where((mag > magCut))[0]

    xlim = [14.0, 25.0]
    ylim3 = [0.001, 0.99]

    #####

    sp1 = fig.add_subplot(221)
    sp1.plot(mag[whereCut], diff[whereCut], "r.", ms=size, label=ccd)
    sp1.plot(mag[whereOther], diff[whereOther], "k.", ms=size, label=ccd)
    sp1.set_ylabel(tag, fontsize = 10)

    #####
    sp2 = fig.add_subplot(222, sharex = sp1)

    def noNeg(xIn):
        x = numpy.array(xIn)
        if len(x) > 1:
            xMax = x.max()
        else:
            xMax = 1.0e-4
        return x.clip(1.0e-5, xMax)


    sp2.plot(mag[whereCut],   noNeg(derr[whereCut]), "r.", ms=size, label=ccd)
    sp2.plot(mag[whereOther], noNeg(derr[whereOther]), "k.", ms=size, label=ccd)
    sp2.set_ylabel('Error Bars', fontsize = 10)

    #####
    sp3 = fig.add_subplot(223, sharex = sp1)

    binmag  = []
    binstd  = []
    binmerr = []
    xmin = xlim[0]
    xmax = xlim[1]
    bins1 = numpy.arange(xmin, xmax, 0.5)
    for i in range(1, len(bins1)):
        idx = numpy.where((mag>bins1[i-1])&(mag<=bins1[i]))[0]
        if len(idx) == 0:
            continue
        #avgMag  = afwMath.makeStatistics(mag[idx], afwMath.MEAN).getValue(afwMath.MEAN)
        #stdDmag = 0.741 * afwMath.makeStatistics(diff[idx], afwMath.IQRANGE).getValue(afwMath.IQRANGE)
        #avgEbar = afwMath.makeStatistics(derr[idx], afwMath.MEAN).getValue(afwMath.MEAN)
        avgMag = mag[idx].mean()
        stdDmag = diff[idx].std()
        avgEbar = derr[idx].mean()
        binmag.append(avgMag)
        binstd.append(stdDmag)
        binmerr.append(avgEbar)
    # Shows the 2 curves   
    sp3.plot(binmag, noNeg(binstd), 'r-', label="Phot RMS")
    sp3.plot(binmag, noNeg(binmerr), 'b--', label="Avg Error Bar")
    sp3.set_xlabel(tag1, fontsize = 10)

    #####
    sp4 = fig.add_subplot(224, sharex = sp1)

    binmag = numpy.array(binmag)
    binstd = numpy.array(binstd)
    binmerr = numpy.array(binmerr)


    idx         = numpy.where( (binstd > binmerr) )[0]
    errbarmag   = binmag[idx]
    errbarresid = numpy.sqrt(binstd[idx]**2 - binmerr[idx]**2)

    whereCut    = numpy.where((errbarmag < magCut))[0]
    whereOther  = numpy.where((errbarmag > magCut))[0]

    sp4.plot(errbarmag[whereCut], noNeg(errbarresid[whereCut]), 'ro', ms = 3, label="Err Underestimate")
    sp4.plot(errbarmag[whereOther], noNeg(errbarresid[whereOther]), 'ko', ms = 3)
    sp4.set_xlabel(tag1, fontsize = 10)

    #### CONFIG
    sp2.semilogy()
    sp3.semilogy()
    sp4.semilogy()

    sp2.yaxis.set_label_position('right')
    sp2.yaxis.set_ticks_position('right')
    sp4.yaxis.set_label_position('right')
    sp4.yaxis.set_ticks_position('right')

    sp3.legend(prop=fm.FontProperties(size="xx-small"), loc="upper left")
    sp4.legend(prop=fm.FontProperties(size="xx-small"), loc="upper left")

    qaPlotUtil.qaSetp(sp1.get_xticklabels()+sp1.get_yticklabels(), fontsize=8)
    qaPlotUtil.qaSetp(sp2.get_xticklabels()+sp2.get_yticklabels(), fontsize=8)
    qaPlotUtil.qaSetp(sp3.get_xticklabels()+sp3.get_yticklabels(), fontsize=8)
    qaPlotUtil.qaSetp(sp4.get_xticklabels()+sp4.get_yticklabels(), fontsize=8)


    sp1.set_xlim(xlim)
    sp1.set_ylim(ylim)

    sp2.set_ylim(ylim3)
    sp3.set_ylim(ylim3)
    sp4.set_ylim(ylim3)

    return fig



def summaryFigure(*summaryArgs):

    # many of these are not needed, but for simplicity all the same args get passed
    # to the plot functions standardFigure(), derrFigure(), and summaryFigure()
    mag0, diff0, star0, derr0, areaLabel, raft, ccd, figsize, xlim, ylim, xlim2, ylim2, ylimStep, \
        tag1, tag, mode, x, y, trend, magCut = summaryArgs

    allMags = mag0
    allDiffs = diff0
    allStars = star0
    allDerr = derr0
    
    size = 2.0

    # dmag vs mag
    fig0 = figure.Figure(figsize=figsize)
    canvan = FigCanvas(fig0)
    
    fig0.subplots_adjust(left=0.125, bottom=0.125)
    rect0_1  = [0.125, 0.35, 0.478-0.125, 0.9-0.35]
    rect0_2  = [0.548, 0.35, 0.9-0.548, 0.9-0.35]
    rect0_1b = [0.125, 0.125, 0.478-0.125, 0.2]
    rect0_2b = [0.548, 0.125, 0.9-0.548, 0.2]
    ax0_1 = fig0.add_axes(rect0_1)
    ax0_2 = fig0.add_axes(rect0_2)
    ax0_1b = fig0.add_axes(rect0_1b, sharex = ax0_1)
    ax0_2b = fig0.add_axes(rect0_2b, sharex = ax0_2)

    w = numpy.where( (allMags < magCut) & (allStars > 0))[0]
    wStar = numpy.where(allStars > 0)[0]
    wGxy = numpy.where(allStars == 0)[0]

    if len(w) > 0:
        lineFit = qaAnaUtil.robustPolyFit(allMags[w], allDiffs[w], 1)
        trendCoeffs = lineFit[0], lineFit[2]
        trendCoeffsLo = lineFit[0]+lineFit[1], lineFit[2]-lineFit[3]
        trendCoeffsHi = lineFit[0]-lineFit[1], lineFit[2]+lineFit[3]
    else:
        trendCoeffs = [0.0, 0.0]
        trendCoeffsLo = [0.0, 0.0]
        trendCoeffsHi = [0.0, 0.0]

    ####################
    # data for all ccds
    haveData = True
    if len(allMags) == 0:
        haveData = False
        allMags = numpy.array([xlim[1]])
        allDiffs = numpy.array([0.0])
        #allColor = [black]
        #allLabels = ["no_valid_data"]
        trendCoeffs = [0.0, 0.0]
        trendCoeffsLo = [0.0, 0.0]
        trendCoeffsHi = [0.0, 0.0]
    else:
        #
        # Lower plots

        mStar = allMags[wStar]
        dStar = allDiffs[wStar]
        eStar = allDerr[wStar]

        binmag  = []
        binstd  = []
        binmerr = []
        xmin = xlim[0]
        xmax = xlim[1]
        bins1 = numpy.arange(xmin, xmax, 0.5)
        for i in range(1, len(bins1)):
            idx = numpy.where((mStar>bins1[i-1])&(mStar<bins1[i]))[0]
            if len(idx) == 0:
                continue
            #avgMag  = afwMath.makeStatistics(mStar[idx], afwMath.MEAN).getValue(afwMath.MEAN)
            avgMag = mStar[idx].mean()
            #stdDmag = 0.741 * afwMath.makeStatistics(dStar[idx], afwMath.IQRANGE).getValue(afwMath.IQRANGE)
            stdDmag = dStar[idx].std()
            #avgEbar = afwMath.makeStatistics(eStar[idx], afwMath.MEAN).getValue(afwMath.MEAN)
            avgEbar = eStar[idx].mean()
            binmag.append(avgMag)
            binstd.append(stdDmag)
            binmerr.append(avgEbar)
        # Shows the 2 curves
        ax0_1b.plot(binmag, binstd, 'r-', label="Phot RMS")
        ax0_1b.plot(binmag, binmerr, 'b--', label="Avg Error Bar")

        binmag = numpy.array(binmag)
        binstd = numpy.array(binstd)
        binmerr = numpy.array(binmerr)

        medresid = 0.0
        idx         = numpy.where( (binstd > binmerr) & (binmag < magCut) )[0]
        if len(idx) == 0:
            medresid = 0.0
        else:
            brightmag   = binmag[idx]
            brightresid = numpy.sqrt(binstd[idx]**2 - binmerr[idx]**2)
            medresid    = numpy.median(brightresid)
            if numpy.isnan(medresid) or medresid == None:
               medresid = 0.0

        #label       = "Quad error"
        #comment     = "Median value to add in quadrature to star phot error bars (mag < %.2f)" % (magCut)
        #testSet.addTest( testCode.Test(label, medresid, derrLimits, comment, areaLabel="all"))

        # SECOND SUBPANEL

        binmag  = []
        binstd  = []
        binmerr = []
        xmin2 = xlim2[0]
        xmax2 = xlim2[1]
        bins2 = numpy.arange(xmin2, xmax2, 0.5)
        for i in range(1, len(bins2)):
            idx = numpy.where((mStar>bins2[i-1])&(mStar<bins2[i]))[0]
            if len(idx) == 0:
                continue
            #avgMag  = afwMath.makeStatistics(mStar[idx], afwMath.MEAN).getValue(afwMath.MEAN)
            #stdDmag = 0.741 * afwMath.makeStatistics(dStar[idx], afwMath.IQRANGE).getValue(afwMath.IQRANGE)
            #avgEbar = afwMath.makeStatistics(eStar[idx], afwMath.MEAN).getValue(afwMath.MEAN)
            avgMag = mStar[idx].mean()
            stdDmag = dStar[idx].std()
            avgEbar = eStar[idx].mean()
            binmag.append(avgMag)
            binstd.append(stdDmag)
            binmerr.append(avgEbar)

        binmag = numpy.array(binmag)
        binstd = numpy.array(binstd)
        binmerr = numpy.array(binmerr)

        idx         = numpy.where( (binstd > binmerr) )[0]
        errbarmag   = binmag[idx]
        errbarresid = numpy.sqrt(binstd[idx]**2 - binmerr[idx]**2)
        ax0_2b.plot(errbarmag, errbarresid, 'ko', ms = 3, label="Err Underestimate")


        # Lower plots
        #

    #allColor = numpy.array(allColor)
    for ax in [ax0_1, ax0_2]:
        ax.plot(xlim2, [0.0, 0.0], "-k", lw=1.0)  # show an x-axis at y=0
        ax.plot(allMags[wGxy], allDiffs[wGxy], "bo", ms=size, label="galaxies", alpha=0.5)
        ax.plot(allMags[wStar], allDiffs[wStar], "ro", ms=size, label="stars", alpha=0.5)
        # 99 is the 'no-data' values
        if abs(trendCoeffs[0] - 99.0) > 1.0e-6:
            ax.plot(xlim2, numpy.polyval(trendCoeffs, xlim2), "-k", lw=1.0)
            ax.plot(xlim2, numpy.polyval(trendCoeffsLo, xlim2), "--k", lw=1.0)
            ax.plot(xlim2, numpy.polyval(trendCoeffsHi, xlim2), "--k", lw=1.0)

    # show outliers railed at the ylims
    wStarOutlier = numpy.where((numpy.abs(allDiffs) > ylim2[1]) & (allStars > 0))[0]
    wGxyOutlier = numpy.where((numpy.abs(allDiffs) > ylim2[1]) & (allStars == 0))[0]
    clip = 0.99
    ax0_2.plot(allMags[wStarOutlier], numpy.clip(allDiffs[wStarOutlier], clip*ylim2[0], clip*ylim2[1]),
               'r.', ms=size)
    ax0_2.plot(allMags[wGxyOutlier], numpy.clip(allDiffs[wGxyOutlier], clip*ylim2[0], clip*ylim2[1]),
               'g.', ms=size)
    ax0_2.legend(prop=fm.FontProperties(size="xx-small"), loc="upper left")

    dmag = 0.1
    ddiff1 = 0.01
    ddiff2 = ddiff1*(ylim2[1]-ylim2[0])/(ylim[1]-ylim[0]) # rescale for larger y range
    if False:
        for j in xrange(len(allMags)):
            area = (allMags[j]-dmag, allDiffs[j]-ddiff1, allMags[j]+dmag, allDiffs[j]+ddiff1)
            fig0.addMapArea(allLabels[j], area, "%.3f_%.3f"% (allMags[j], allDiffs[j]), axes=ax0_1)
            area = (allMags[j]-dmag, allDiffs[j]-ddiff2, allMags[j]+dmag, allDiffs[j]+ddiff2)
            fig0.addMapArea(allLabels[j], area, "%.3f_%.3f"% (allMags[j], allDiffs[j]), axes=ax0_2)


    del allMags
    del allDiffs
    #del allColor
    #del allLabels

    # move the yaxis ticks/labels to the other side
    ax0_2.yaxis.set_label_position('right')
    ax0_2.yaxis.set_ticks_position('right')
    ax0_2b.yaxis.set_label_position('right')
    ax0_2b.yaxis.set_ticks_position('right')

    ax0_1b.set_xlabel(tag1, fontsize=12)
    ax0_2b.set_xlabel(tag1, fontsize=12)

    ax0_2.plot([xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]],
               [ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]], '-k')

    ax0_1b.semilogy()
    ax0_2b.semilogy()

    if haveData:
        ax0_1b.legend(prop=fm.FontProperties(size="xx-small"), loc="upper left")
        ax0_2b.legend(prop=fm.FontProperties(size="xx-small"), loc="upper left")


    qaPlotUtil.qaSetp(ax0_1.get_xticklabels()+ax0_2.get_xticklabels(), visible=False)
    qaPlotUtil.qaSetp(ax0_1.get_yticklabels()+ax0_2.get_yticklabels(), fontsize=11)
    qaPlotUtil.qaSetp(ax0_1b.get_yticklabels()+ax0_2b.get_yticklabels(), fontsize=10)

    for ax in [ax0_1, ax0_2]:
        ax.set_ylabel(tag)

    ax0_1.set_xlim(xlim)
    ax0_2.set_xlim(xlim2)
    ax0_1.set_ylim(ylim)
    ax0_2.set_ylim(ylim2)

    ax0_1b.set_ylim(0.001, 0.99)
    ax0_2b.set_ylim(0.001, 0.99)

    return fig0





if __name__ == '__main__':
    filename, mode, figType = sys.argv[1:4]
    data, isSummary = qaPlotUtil.unshelveGlob(filename, flag='r')
    data['mode'] = mode
    data['figType'] = figType
    fig = plot(data)
    fig.savefig(filename)
