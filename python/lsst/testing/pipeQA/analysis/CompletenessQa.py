import re
import numpy as num

import time
import lsst.testing.pipeQA.TestCode as testCode
import QaAnalysis as qaAna
import lsst.testing.pipeQA.figures as qaFig
import lsst.meas.algorithms as measAlg
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtils
import RaftCcdData as raftCcdData

import lsst.testing.pipeQA.source as pqaSource

import matplotlib.ticker as ticker
from matplotlib.font_manager import FontProperties
import matplotlib.patches as patches

# Until we can make it more robust
hasMinuit = False
try:
    import minuit2
except:
    hasMinuit = False
    

class CompletenessQa(qaAna.QaAnalysis):
    def __init__(self, completenessMagMin, completenessMagMax, **kwargs):
        qaAna.QaAnalysis.__init__(self, **kwargs)
        self.limits = [completenessMagMin, completenessMagMax]
        self.bins   = num.arange(14, 27, 0.5)

        self.description = """
         For each CCD, the top figure shows the counts of 4 classes of
         objects, plotted in histograms as a function of magnitude, for stars
         only.  "Matched" objects are detections that match 1-to-1 with with
         reference catalog (through Source Association); "Blended" objects
         match N-to-1 with the reference catalog; "Orphan" objects do not match
         any entry in the reference catalog, and may include false positives
         and asteroids; while "Unmatched" stars are those present in the
         reference catalog but not detected in the science image.  The measure
         of completeness is represented by (matched + blended) / (matched +
         blended + unmatched), which is the blue line.  The magnitude below
         which the completeness is < 50% is indicated with the crosshairs.  The
         bottom panel shows the same for galaxies only.  Orphans are included
         in both plots.  The summary FPA figure provides a visual
         representation of this photometric depth.
        """
         
    def free(self):
        del self.detector
        del self.filter
        del self.matchListDictSrc

        del self.orphan
        del self.matchedStar
        del self.blendedStar
        del self.undetectedStar
        del self.matchedGalaxy
        del self.blendedGalaxy
        del self.undetectedGalaxy
        del self.depth
        if hasMinuit:
            del self.fit
        del self.faintest

    def limitingMag(self, raftId, ccdId):
        if hasMinuit:
            try:
                return self.limitingMagMinuit(raftId, ccdId)
            except:
                pass

        matchedStar     = num.array(self.matchedStar.get(raftId, ccdId))
        blendedStar     = num.array(self.blendedStar.get(raftId, ccdId))
        undetectedStar  = num.array(self.undetectedStar.get(raftId, ccdId))

        allStars        = num.concatenate((matchedStar, blendedStar, undetectedStar))
        foundStars      = num.concatenate((matchedStar, blendedStar))
        histAll         = num.histogram(allStars, bins=self.bins)
        histFound       = num.histogram(foundStars, bins=self.bins)

        magbins = 0.5 * (histAll[1][1:] + histAll[1][:-1])
        w       = num.where(histAll[0] != 0)
        x       = magbins[w]
        n       = 1.0 * histFound[0][w]
        d       = 1.0 * histAll[0][w]
        y       = n / d

        binsize = self.bins[1] - self.bins[0]
        x = num.append(x, x[-1] + binsize)
        y = num.append(y, 0.0)

        for i in num.arange(len(y) - 1, 1, -1):
            if y[i] <= 0.5 and y[i-1] > 0.5:
                return (0.5 - y[i-1]) / (y[i] - y[i-1]) * (x[i] - x[i-1]) + x[i-1]
        return 0.0
        
    def limitingMagMinuit(self, raftId, ccdId):
        # Model is of the form:
        # 0.5 + -1.0 / num.pi * num.arctan(A * x + B)
        
        import minuit2
        matchStarList   = self.matchStarSrc.get(raftId, ccdId)
        unmatchStarList = self.unmatchCatStar.get(raftId, ccdId)
        
        histStarSrc     = num.histogram(matchStarList, bins = self.bins)
        maxSrcIdx       = num.argsort(histStarSrc[0])[-1]
        histUnmatchStar = num.histogram(unmatchStarList, bins = self.bins)

        magbins   = 0.5 * (histStarSrc[1][1:] + histStarSrc[1][:-1])
        histRatio = num.zeros(len(histStarSrc[0]))

        w     = num.where((histStarSrc[0] + histUnmatchStar[0]) != 0)
        x     = magbins[w]

        # approximate
        d     = 1.0 * histStarSrc[0][w]
        u     = 1.0 * histUnmatchStar[0][w]
        dd    = num.sqrt(d)
        n     = d + u
        y     = d / n
        dy    = dd / n

        idx = num.where(dy != 0)
        x   = x[idx]
        y   = y[idx]
        dy  = dy[idx]
              
        def fcn(A, B):
            model  = 0.5 + -1.0 / num.pi * num.arctan(A * x + B)
            chi    = (model - y) / dy
            return num.sum(chi**2)
        
        m = minuit2.Minuit2(fcn)
        m.values['A'] = 1
        m.values['B'] = -10
        m.migrad()

        mx = num.arange(min(x), max(x), 0.1)
        my = 0.5 + -1.0 / num.pi * num.arctan(m.values['A'] * mx + m.values['B'])
        mindx = num.argsort((num.abs(my-0.5)))[0]

        self.fit.set(raftId, ccdId, [m.values['A'], m.values['B']])
        return mx[mindx]

    def test(self, data, dataId, fluxType = "psf"):

        testSet = self.getTestSet(data, dataId)
        testSet.addMetadata({"Description": self.description})
        
        self.fluxType = fluxType
        self.detector         = data.getDetectorBySensor(dataId)
        self.filter           = data.getFilterBySensor(dataId)
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')

        self.orphan           = raftCcdData.RaftCcdVector(self.detector)
        self.matchedStar      = raftCcdData.RaftCcdVector(self.detector)
        self.blendedStar      = raftCcdData.RaftCcdVector(self.detector)
        self.undetectedStar   = raftCcdData.RaftCcdVector(self.detector)
        self.matchedGalaxy    = raftCcdData.RaftCcdVector(self.detector)
        self.blendedGalaxy    = raftCcdData.RaftCcdVector(self.detector)
        self.undetectedGalaxy = raftCcdData.RaftCcdVector(self.detector)
        self.depth            = raftCcdData.RaftCcdData(self.detector)

        if hasMinuit:
            self.fit = raftCcdData.RaftCcdData(self.detector, initValue=[0.0, 0.0]) 



        sCatDummy = pqaSource.Catalog().catalog
        sCatSchema = sCatDummy.getSchema()
        srefCatDummy  = pqaSource.RefCatalog().catalog
        srefCatSchema = srefCatDummy.getSchema()
        
        psfKey = sCatSchema.find('PsfFlux').key
        psfErrKey = sCatSchema.find('PsfFluxErr').key
        apKey = sCatSchema.find('ApFlux').key
        apErrKey = sCatSchema.find('ApFluxErr').key
        extKey = sCatSchema.find('Extendedness').key
        
        refPsfKey = srefCatSchema.find('PsfFlux').key

        self.faintest = 0.0
        for key in self.detector.keys():
            raftId     = self.detector[key].getParent().getId().getName()
            ccdId      = self.detector[key].getId().getName()
            filterName = self.filter[key].getName()

            if self.matchListDictSrc.has_key(key):
                # Detections
                matchSet = [
                    [self.matchListDictSrc[key]['matched'], self.matchedStar, self.matchedGalaxy],
                    [self.matchListDictSrc[key]['blended'], self.blendedStar, self.blendedGalaxy]
                    ]
                for mset in matchSet:
                    mdict, starvec, galvec = mset
                    
                    stars    = []
                    galaxies = []
                    for m in mdict:
                        sref, s, dist = m
                        if fluxType == "psf":
                            fref  = sref.getD(refPsfKey)
                            f     = s.getD(psfKey)
                            ferr  = s.getD(psfErrKey)
                        else:
                            fref  = sref.getD(refPsfKey)
                            f     = s.getD(apKey)
                            ferr  = s.getD(apErrKey)


                        if (fref > 0.0 and f > 0.0):
                            # Use known catalog mag
                            mrefmag  = -2.5*num.log10(fref)
                            if num.isfinite(mrefmag):
                                if mrefmag > self.faintest:
                                    self.faintest == mrefmag
                                    
                                if s.getD(extKey) > 0.0:
                                    galaxies.append(mrefmag)
                                else:
                                    stars.append(mrefmag)

                    starvec.set(raftId, ccdId, num.array(stars))
                    galvec.set(raftId, ccdId, num.array(galaxies))
    
                # Non-detections
                undetectedStars = []
                undetectedGalaxies = []

                for nondet in self.matchListDictSrc[key]['undetected']:
                    mag = nondet.getMag(filterName)
                    if mag > self.faintest:
                        self.faintest = mag
                    if nondet.getIsStar():
                        undetectedStars.append(mag)
                    else:
                        undetectedGalaxies.append(mag)
                self.undetectedStar.set(raftId, ccdId, num.array(undetectedStars))
                self.undetectedGalaxy.set(raftId, ccdId, num.array(undetectedGalaxies))
                    
                # Orphans
                orphans = []
                for orphan in self.matchListDictSrc[key]['orphan']:
                    if self.fluxType == "psf":
                        f = orphan.getD(psfKey)
                    else:
                        f = orphan.getD(apKey)
                    if f > 0.0:
                        orphmag = -2.5*num.log10(f)
                        orphans.append(orphmag)
                        if orphmag > self.faintest:
                            self.faintest = orphmag
                self.orphan.set(raftId, ccdId, num.array(orphans))

                ############ Calculate limiting mag
                
                maxDepth = self.limitingMag(raftId, ccdId)
                self.depth.set(raftId, ccdId, maxDepth)
                
                areaLabel = data.cameraInfo.getDetectorName(raftId, ccdId)
                label = "photometric depth "
                comment = "magnitude where star completeness drops below 0.5"
                test = testCode.Test(label, maxDepth, self.limits, comment, areaLabel=areaLabel)
                testSet.addTest(test)
                
    def plot(self, data, dataId, showUndefined = False):
        
        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)
        isFinalDataId = False
        if len(data.brokenDataIdList) == 0 or data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True

        # override self.limits to use the faintest mag as an upper limit, if upper limit is 99.0
        # put the new limits in a different variable so we do this on a visit-by-visit basis
        # ie. ... don't set it permanently in self.limits or every visit will use that.
        limitsToUse = [self.limits[0], self.limits[1]]
        if abs(self.limits[1] - 99.0) < 1.0e-3:
            limitsToUse[1] = self.faintest
            
        # fpa figure
        filebase = "completenessDepth"
        depthData, depthMap = testSet.unpickle(filebase, default=[None, None])
        depths = []
        depthFig = qaFig.FpaQaFigure(data.cameraInfo, data=depthData, map=depthMap)
        for raft, ccdDict in depthFig.data.items():
            for ccd, value in ccdDict.items():
                if not self.depth.get(raft, ccd) is None:
                    depth = self.depth.get(raft, ccd)
                    depths.append(depth)
                    depthFig.data[raft][ccd] = depth
                    if num.isfinite(depth):
                        depthFig.map[raft][ccd] = 'mag=%.2f'%(depth)
                    else:
                        depthFig.map[raft][ccd] = 'mag=nan'
                        
        testSet.pickle(filebase, [depthFig.data, depthFig.map])

        blue = '#0000ff'
        red  = '#ff0000'
        
        if False:
            if len(depths) >= 2:
                vmin = max(num.min(depths), self.limits[0])
                vmax = min(num.max(depths), self.limits[1])
            else:
                vmin = self.limits[0]
                vmax = self.limits[1]

            if vmax <= vmin:
                vmin = self.limits[0]
                vmax = self.limits[1]
        else:
            vmin, vmax = 1.0*limitsToUse[0], 1.0*limitsToUse[1]

        if not self.delaySummary or isFinalDataId:
            print "plotting FPAs"
            depthFig.makeFigure(showUndefined=showUndefined, cmap="RdBu_r", vlimits=[vmin, vmax],
                                title="Photometric Depth", cmapOver=red, cmapUnder=blue,
                                failLimits=limitsToUse)
            testSet.addFigure(depthFig, filebase+".png", "Estimate of photometric depth",  navMap=True)

        del depthFig


        cacheLabel = "completeness"
        shelfData = {}
        
        # Each CCD
        for raft, ccd in self.depth.raftCcdKeys():
            orphan           = self.orphan.get(raft, ccd)
            matchedStar      = self.matchedStar.get(raft, ccd)
            blendedStar      = self.blendedStar.get(raft, ccd)
            undetectedStar   = self.undetectedStar.get(raft, ccd)
            matchedGalaxy    = self.matchedGalaxy.get(raft, ccd)
            blendedGalaxy    = self.blendedGalaxy.get(raft, ccd)
            undetectedGalaxy = self.undetectedGalaxy.get(raft, ccd)
            depth            = self.depth.get(raft, ccd)

            print "Plotting ", ccd
            label = data.cameraInfo.getDetectorName(raft, ccd)
            fig = self.standardFigure(label, orphan, depth,
                                      matchedStar, blendedStar, undetectedStar,
                                      matchedGalaxy, blendedGalaxy, undetectedGalaxy)

            testSet.addFigure(fig, "completeness.png", "Photometric detections "+label, areaLabel=label)
            del fig


            shelfData[ccd] = [orphan, depth,
                              matchedStar, blendedStar, undetectedStar,
                              matchedGalaxy, blendedGalaxy, undetectedGalaxy]
            
        if self.useCache:
            testSet.shelve(cacheLabel, shelfData)

        if not self.delaySummary or isFinalDataId:
            print "plotting Summary figure"

            # unstash the values
            if self.useCache:
                shelfData = testSet.unshelve(cacheLabel)

            orphanAll           = num.array([])
            depthAll            = num.array([])
            matchedStarAll      = num.array([])
            blendedStarAll      = num.array([])
            undetectedStarAll   = num.array([])
            matchedGalaxyAll    = num.array([])
            blendedGalaxyAll    = num.array([])
            undetectedGalaxyAll = num.array([])
            for k,v in shelfData.items():
                orphan, depth = v[0:2]
                matchedStar, blendedStar, undetectedStar = v[2:5]
                matchedGalaxy, blendedGalaxy, undetectedGalaxy = v[5:8]
                orphanAll           = num.append(orphanAll           , orphan)
                depthAll            = num.append(depthAll            , depth)
                matchedStarAll      = num.append(matchedStarAll      , matchedStar)
                blendedStarAll      = num.append(blendedStarAll      , blendedStar)
                undetectedStarAll   = num.append(undetectedStarAll   , undetectedStar)
                matchedGalaxyAll    = num.append(matchedGalaxyAll    , matchedGalaxy)
                blendedGalaxyAll    = num.append(blendedGalaxyAll    , blendedGalaxy)
                undetectedGalaxyAll = num.append(undetectedGalaxyAll , undetectedGalaxy)

            allFig = self.standardFigure("All Sensors", orphanAll, depthAll.mean(),
                                         matchedStarAll, blendedStarAll, undetectedStarAll,
                                         matchedGalaxyAll, blendedGalaxyAll, undetectedGalaxyAll)

            del orphanAll, depthAll, \
                matchedStarAll, blendedStarAll, undetectedStarAll, \
                matchedGalaxyAll, blendedGalaxyAll, undetectedGalaxyAll
            
            label = "all"
            testSet.addFigure(allFig, "completeness.png", "Photometric detections "+label, areaLabel=label)
            del allFig
            

    def standardFigure(self, title, orphan, depth,
                       matchedStar, blendedStar, undetectedStar,
                       matchedGalaxy, blendedGalaxy, undetectedGalaxy):


        fig = qaFig.QaFigure()
        sp1 = fig.fig.add_subplot(211)
        fig.fig.subplots_adjust(left=0.13)
        sp2 = fig.fig.add_subplot(212, sharex = sp1)

        # Stacked histogram
        orphanHist         = num.histogram(orphan, bins=self.bins)
        matchedStarHist    = num.histogram(matchedStar, bins=self.bins)
        blendedStarHist    = num.histogram(blendedStar, bins=self.bins)
        undetectedStarHist = num.histogram(undetectedStar, bins=self.bins)
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
        histAll         = num.histogram(allStars, bins=self.bins)
        histFound       = num.histogram(foundStars, bins=self.bins)

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
        
        qaFigUtils.qaSetp(sp1x2.get_xticklabels(), visible=False)
        qaFigUtils.qaSetp(sp1x2.get_yticklabels(), fontsize = 6)

        sp1.set_ylabel('N Stars', fontsize=10)
        if False: #1.4*ymax > 1000:
            sp1.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2g"))
            for t in sp1.get_yticklabels():
                print t.get_text()
                t.set_text(re.sub("\+0", "", t.get_text()))
        qaFigUtils.qaSetp(sp1.get_xticklabels()+sp1.get_yticklabels(), fontsize = 6)

        ##############

        orphanHist        = num.histogram(orphan, bins=self.bins)
        matchedGalHist    = num.histogram(matchedGalaxy, bins=self.bins)
        blendedGalHist    = num.histogram(blendedGalaxy, bins=self.bins)
        undetectedGalHist = num.histogram(undetectedGalaxy, bins=self.bins)
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
        qaFigUtils.qaSetp(sp2.get_xticklabels()+sp2.get_yticklabels(), fontsize = 6)
        qaFigUtils.qaSetp(sp2.get_yticklabels(), rotation = 45.0)
        sp2.legend(numpoints = 1, prop=FontProperties(size='x-small'), loc = 'upper left')
        #sp2.set_ylim(0.75, 999)
        #sp2.semilogy()

        sp1.set_xlim(14, 26)

        fig.fig.suptitle('%s Stacked histogram' % (title), fontsize = 11)

        return fig


