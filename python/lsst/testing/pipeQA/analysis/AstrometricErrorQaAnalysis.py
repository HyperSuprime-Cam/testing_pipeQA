import sys, os, re
import lsst.meas.algorithms        as measAlg
import lsst.testing.pipeQA.figures as qaFig
import numpy

import time

import lsst.afw.math                as afwMath
import lsst.testing.pipeQA.TestCode as testCode

import QaAnalysis as qaAna
import RaftCcdData as raftCcdData
import QaAnalysisUtils as qaAnaUtil
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtil

import lsst.testing.pipeQA.source as pqaSource

import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.font_manager as fm
from  matplotlib.ticker import MaxNLocator
from matplotlib.collections import LineCollection
from matplotlib.patches import Circle

class AstrometricErrorQaAnalysis(qaAna.QaAnalysis):

    def __init__(self, maxErr, **kwargs):
        qaAna.QaAnalysis.__init__(self, **kwargs)
        self.limits = [0.0, maxErr]

        self.description = """
         For each CCD, the left figure shows the distance offset between the
         measured centroid of matched objects and the catalog position of these
         objects, represented as an arrow.  The top right panel provides the
         view of these vectors stacked at the position of the reference object,
         with the green circle representing the radius that contains 50% of the
         matches.  The bottom panel provides a histogram of the astrometric
         offsets, with the median indicated.  The summary FPA figure provides
         the median vector (offset and angle) for each chip.
        """
        
    def free(self):
        del self.x
        del self.y
        del self.dDec
        del self.dRa
        del self.filter
        
        del self.detector
        del self.matchListDictSrc

        del self.medErrArcsec
        del self.medThetaRad

        del self.sCatDummy
        del self.srefCatDummy
        
    def test(self, data, dataId):

        # get data
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')
        self.detector         = data.getDetectorBySensor(dataId)
        self.filter           = data.getFilterBySensor(dataId)
        
        #self.clusters = data.getSourceClusters(dataId)

        # compute the mean ra, dec for each source cluster

        self.dRa  = raftCcdData.RaftCcdVector(self.detector)
        self.dDec = raftCcdData.RaftCcdVector(self.detector)
        self.x    = raftCcdData.RaftCcdVector(self.detector)
        self.y    = raftCcdData.RaftCcdVector(self.detector)

        self.sCatDummy = pqaSource.Catalog()
        sCatDummy = self.sCatDummy.catalog
        sCatSchema = sCatDummy.getSchema()
        self.srefCatDummy  = pqaSource.RefCatalog()
        srefCatDummy = self.srefCatDummy.catalog
        srefCatSchema = srefCatDummy.getSchema()
        
        xKey      = sCatSchema.find('XAstrom').key
        yKey      = sCatSchema.find('YAstrom').key
        raKey     = sCatSchema.find('Ra').key
        decKey    = sCatSchema.find('Dec').key
        refRaKey  = srefCatSchema.find('Ra').key
        refDecKey = srefCatSchema.find('Dec').key

        filter = None
        for key in self.matchListDictSrc.keys():
            raft = self.detector[key].getParent().getId().getName()
            ccd  = self.detector[key].getId().getName()
            filter = self.filter[key].getName()

            matchList = self.matchListDictSrc[key]['matched']
            for m in matchList:
                sref, s, dist = m
                ra, dec, raRef, decRef = \
                    [numpy.radians(x) for x in [s.getD(raKey), s.getD(decKey),
                                                sref.getD(refRaKey), sref.getD(refDecKey)]]
                
                dDec = decRef - dec
                dRa  = (raRef - ra)*abs(numpy.cos(decRef))
                
                if not (s.getD(sCatSchema.find('FlagPixInterpCen').key)):
                    self.dRa.append(raft, ccd, dRa)
                    self.dDec.append(raft, ccd, dDec)
                    self.x.append(raft, ccd, s.getD(xKey))
                    self.y.append(raft, ccd, s.getD(yKey))
                    
                    
        testSet = self.getTestSet(data, dataId)
        testSet.addMetadata({"Description": self.description})

        self.medErrArcsec = raftCcdData.RaftCcdData(self.detector)
        self.medThetaRad  = raftCcdData.RaftCcdData(self.detector)

        for raft,  ccd in self.dRa.raftCcdKeys():
            dRa  = self.dRa.get(raft, ccd)
            dDec = self.dDec.get(raft, ccd)
            
            errArcsec = 206265.0*numpy.sqrt(dRa**2 + dDec**2)
            #errArcsec = 3600.0*numpy.sqrt(dRa**2 + dDec**2)
            thetaRad  = numpy.arctan2(dDec, dRa)

            if len(errArcsec) > 0:
                stat  = afwMath.makeStatistics(errArcsec, afwMath.NPOINT | afwMath.MEDIAN)
                medErrArcsec = stat.getValue(afwMath.MEDIAN)
                stat  = afwMath.makeStatistics(thetaRad, afwMath.NPOINT | afwMath.MEDIAN)
                medThetaRad = stat.getValue(afwMath.MEDIAN)
                n = stat.getValue(afwMath.NPOINT)
            else:
                medErrArcsec = -1.0
                medThetaRad = 0.0
                n = 0

            self.medErrArcsec.set(raft, ccd, medErrArcsec)
            self.medThetaRad.set(raft, ccd, medThetaRad)
            
            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            label = "median astrometry error "
            comment = "median sqrt(dRa^2+dDec^2) (arcsec, nstar=%d)" % (n)
            test = testCode.Test(label, medErrArcsec, self.limits, comment, areaLabel=areaLabel)
            testSet.addTest(test)

        
    def plot(self, data, dataId, showUndefined=False):

        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)


        isFinalDataId = False
        if len(data.brokenDataIdList) == 0 or data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True

        medAstBase = "medAstError"
        medAstData, medAstMap = testSet.unpickle(medAstBase, default=[None, None])

        # fpa figure
        astFig = qaFig.VectorFpaQaFigure(data.cameraInfo, data=medAstData, map=medAstMap)

        vLen = 5000 # length in pixels for 1 arcsec error vector
        for raft, ccdDict in astFig.data.items():
            for ccd, value in ccdDict.items():
                if not self.medErrArcsec.get(raft, ccd) is None:
                    astErrArcsec = self.medErrArcsec.get(raft, ccd)
                    thetaRad = self.medThetaRad.get(raft, ccd)
                    astFig.data[raft][ccd] = [thetaRad, vLen*astErrArcsec, astErrArcsec]
                    astFig.map[raft][ccd] = "\"/theta=%.2f/%.0f" % (astErrArcsec, numpy.degrees(thetaRad))
                
        testSet.pickle(medAstBase, [astFig.data, astFig.map])
        
        if not self.delaySummary or isFinalDataId:
            print "plotting FPAs"
            astFig.makeFigure(showUndefined=showUndefined, cmap="Reds", vlimits=[0.0, 2.0*self.limits[1]],
                              title="Median astrometric error", cmapOver='#ff0000', failLimits=self.limits,
                              cmapUnder="#ff0000")
            testSet.addFigure(astFig, medAstBase+".png", "Median astrometric error",  navMap=True)
            
        del astFig


        cacheLabel = "astromError"
        shelfData = {}
        
        for raft, ccd in self.dRa.raftCcdKeys():
            ra = self.dRa.get(raft, ccd)
            dec = self.dDec.get(raft, ccd)
            eLen = 206265.0*numpy.sqrt(ra**2 + dec**2)
            #eLen = 3600.0*numpy.sqrt(ra**2 + dec**2)
            t = numpy.arctan2(dec, ra)
            
            dx = eLen*numpy.cos(t)
            w = numpy.where(numpy.abs(dx) < 10.0)
            
            dx = dx[w]
            dy = (eLen*numpy.sin(t))[w]
            x = (self.x.get(raft, ccd))[w]
            y = (self.y.get(raft, ccd))[w]

            print "plotting ", ccd

            fig = self.standardFigure(x, y, dx, dy)
            label = data.cameraInfo.getDetectorName(raft, ccd)
            testSet.addFigure(fig, "astromError.png", "Astrometric error"+label, areaLabel=label)
            del fig

            shelfData[ccd] = [x, y, dx, dy]

        if self.useCache:
            testSet.shelve(cacheLabel, shelfData)

        if not self.delaySummary or isFinalDataId:
            print "plotting Summary figure"

            # unstash the values
            if self.useCache:
                shelfData = testSet.unshelve(cacheLabel)

            xAll  = numpy.array([])
            yAll  = numpy.array([])
            dxAll = numpy.array([])
            dyAll = numpy.array([])
            for k,v in shelfData.items():
                #allCcds.append(k)
                x, y, dx, dy = v
                xAll   = numpy.append(xAll  , x )
                yAll   = numpy.append(yAll  , y )
                dxAll  = numpy.append(dxAll , dx)
                dyAll  = numpy.append(dyAll , dy)

            allFig = self.standardFigure(xAll, yAll, dxAll, dyAll, gridVectors=True)
            del xAll, yAll, dxAll, dyAll
            
            label = "all"
            testSet.addFigure(allFig, "astromError.png", "Astrometric error"+label, areaLabel=label)
            del allFig
            

        
    def standardFigure(self, x, y, dx, dy, gridVectors=False):


        #
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

        # round up to nearest 1024 for limits
        xmax, ymax = x.max(), y.max()
        if len(x) > 1:
            xlim = [0, 1024*int(xmax/1024.0 + 0.5)]
            ylim = [0, 1024*int(ymax/1024.0 + 0.5)]
        else:
            xlim = [0, 1.0]
            ylim = [0, 1.0]

        r = numpy.sqrt(dx**2 + dy**2)
        rmax = 1.2 # r.max()

        fig = qaFig.QaFigure(size=figsize)
        fig.fig.subplots_adjust(left=0.1)

        ################
        # ccd view
        ax = fig.fig.add_subplot(121)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + 0.1*box.height, box.width, 0.9*box.height])

        # for the 'all' plot, just show avg vectors in grid cells
        if gridVectors:
            nx, ny = 8, 8
            xstep, ystep = xlim[1]/nx, ylim[1]/ny
            xgrid = [[0.0]*nx for i in range(ny)]
            ygrid = [[0.0]*nx for i in range(ny)]
            ngrid = [[0]*nx for i in range(ny)]
            for i in range(len(x)):
                ix, iy = int(x[i]/xstep), int(y[i]/ystep)
                if ix >= 0 and ix < nx and iy >= 0 and iy < ny:
                    xgrid[ix][iy] += dx[i]
                    ygrid[ix][iy] += dy[i]
                    ngrid[ix][iy] += 1

            xt, yt, dxt, dyt = [],[],[],[]
            for ix in range(nx):
                for iy in range(ny):
                    xt.append(xstep*(ix+0.5))
                    yt.append(ystep*(iy+0.5))
                    if ngrid[ix][iy] > 0:
                        xval = xgrid[ix][iy]/ngrid[ix][iy]
                        yval = ygrid[ix][iy]/ngrid[ix][iy]
                        dxt.append(xval)
                        dyt.append(yval)
                    else:
                        dxt.append(0.0)
                        dyt.append(0.0)
                        
            #for i in range(len(xt)):
            #    print xt[i], yt[i], dxt[i], dyt[i]
            q = ax.quiver(numpy.array(xt), numpy.array(yt), numpy.array(dxt), numpy.array(dyt),
                          color='k', scale=1.0, angles='xy', pivot='middle', width=0.004)
            ax.quiverkey(q, 0.9, -0.2, 0.1, "100 mas", coordinates='axes',
                         fontproperties={'size':'small'})
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
        ax = fig.fig.add_axes([left, bottom, width, height])

        #get the figure width/heigh in inches to correct
        #aspect ratio
        f_w, f_h = fig.fig.get_size_inches()
        xlimRose = [-width*f_w/(height*f_h)*rmax, f_w*width/(f_h*height)*rmax]
        limRose = [-rmax, rmax]
        # this is much faster than calling plot() in a loop, and quiver() scale length buggy
        if gridVectors:
            ybin = 50
            xbin = 50 #int(1.0*ybin*f_w*width/(f_h*height))
            qaFigUtil.make_densityplot(ax, dx, dy, xlims=limRose, ylims=limRose, bins=(xbin,ybin),
                                       log=True)
            qaFigUtil.make_densityContour(ax, dx, dy, xlims=limRose, ylims=limRose, bins=(xbin,ybin),
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

            #r = numpy.sqrt(dx**2 + dy**2)
            #rmax = r.max()
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
        ax0 = fig.fig.add_axes([left, bottom, width, height])

        rNmax = 1.0
        if len(r) > 1:
            binWidth = 0.5*numpy.std(r)*(20.0/len(r))**0.2
            nBin = (r.max() - r.min())/binWidth
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

