import re
import numpy                                     as num

import lsst.meas.algorithms                      as measAlg
import lsst.pex.config                           as pexConfig
import lsst.pipe.base                            as pipeBase

from   .QaAnalysisTask                           import QaAnalysisTask
import lsst.testing.pipeQA.TestCode              as testCode
import lsst.testing.pipeQA.figures               as qaFig
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtils
import RaftCcdData                               as raftCcdData

import QaPlotUtils                               as qaPlotUtil


# Until we can make it more robust
hasMinuit = False
try:
    import minuit2
except:
    hasMinuit = False
    

class CompletenessQaConfig(pexConfig.Config):
    cameras        = pexConfig.ListField(dtype=str,
                                         doc="Cameras to run CompletenessQaTask",
                                         default=("lsstSim", "cfht", "sdss", "coadd"))
    completeMinMag = pexConfig.Field(dtype=float, doc="Minimum photometric depth", default = 20.0)
    completeMaxMag = pexConfig.Field(dtype=float, doc="Maximum reasonable photometric depth", default = 25.0)

    
class CompletenessQaTask(QaAnalysisTask):
    ConfigClass  = CompletenessQaConfig
    _DefaultName = "completenessQa"

    def __init__(self, **kwargs):
        QaAnalysisTask.__init__(self, **kwargs)
        self.limits = [self.config.completeMinMag, self.config.completeMaxMag]
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
                            fref  = sref.getD(data.k_rPsf)
                            f     = s.getD(data.k_Psf)
                            ferr  = s.get(data.k_PsfE)
                        else:
                            fref  = sref.getD(data.k_rPsf)
                            f     = s.getD(data.k_Ap)
                            ferr  = s.getD(data.k_ApE)


                        if (fref > 0.0 and f > 0.0):
                            # Use known catalog mag
                            mrefmag  = -2.5*num.log10(fref)
                            if num.isfinite(mrefmag):
                                if mrefmag > self.faintest:
                                    self.faintest == mrefmag
                                    
                                if s.get(data.k_ext) > 0.0:
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
                        f = orphan.getD(data.k_Psf)
                    else:
                        f = orphan.getD(data.k_Ap)
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
        depths = []
        depthFig = qaFig.FpaQaFigure(data.cameraInfo, data=depthData, map=depthMap)

        if self.summaryProcessing != self.summOpt['summOnly']:
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
                        label = data.cameraInfo.getDetectorName(raft, ccd)
                        testSet.pickle(filebase+label, [depthFig.data, depthFig.map])

        blue = '#0000ff'
        red  = '#ff0000'

        
        if (self.summaryProcessing in [self.summOpt['summOnly'], self.summOpt['delay']]) and isFinalDataId:

            for raft, ccdDict in depthFig.data.items():
                for ccd, value in ccdDict.items():
                    if not self.depth.get(raft, ccd) is None:
                        label = data.cameraInfo.getDetectorName(raft, ccd)
                        depthDataTmp, depthMapTmp = testSet.unpickle(filebase+label, default=[None, None])
                        depthFig.mergeValues(depthDataTmp, depthMapTmp)
                        
            vmin, vmax = 1.0*limitsToUse[0], 1.0*limitsToUse[1]
            if vmin > vmax:
                vmax = vmin + (self.limits[1] - self.limits[0])
            self.log.log(self.log.INFO, "plotting FPAs")
            depthFig.makeFigure(showUndefined=showUndefined, cmap="RdBu_r", vlimits=[vmin, vmax],
                                title="Photometric Depth", cmapOver=red, cmapUnder=blue,
                                failLimits=limitsToUse)
            testSet.addFigure(depthFig, filebase+".png", "Estimate of photometric depth",  navMap=True)

        del depthFig

        cacheLabel = "completeness"


        if self.summaryProcessing != self.summOpt['summOnly']:
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

                self.log.log(self.log.INFO, "Plotting %s" % (ccd))
                label = data.cameraInfo.getDetectorName(raft, ccd)

                dataDict = {
                    'title'                     :        label           ,
                    'orphan'                    :        orphan          ,
                    'depth'                     :        depth           ,
                    'matchedStar'               :        matchedStar     ,
                    'blendedStar'               :        blendedStar     ,
                    'undetectedStar'            :        undetectedStar  ,
                    'matchedGalaxy'             :        matchedGalaxy   ,
                    'blendedGalaxy'             :        blendedGalaxy   ,
                    'undetectedGalaxy'          :        undetectedGalaxy,
                    'bins'                      :        self.bins.tolist(),
                    }

                import CompletenessQaPlot as plotModule
                caption = "Photometric detections " + label
                pngFile = cacheLabel + ".png"


                if self.lazyPlot.lower() in ['sensor', 'all']:
                    testSet.addLazyFigure(dataDict, pngFile, caption,
                                          plotModule, areaLabel=label, plotargs="")
                else:
                    testSet.cacheLazyData(dataDict, pngFile, areaLabel=label)
                    fig = plotModule.plot(dataDict)
                    testSet.addFigure(fig, pngFile, caption, areaLabel=label)
                    del fig


                
        if (self.summaryProcessing in [self.summOpt['summOnly'], self.summOpt['delay']]) and isFinalDataId:

            self.log.log(self.log.INFO, "plotting Summary figure")

            label = 'all'
            import CompletenessQaPlot as plotModule
            caption = "Photometric detections " + label
            pngFile = "completeness.png"
            
            if self.lazyPlot in ['all']:
                testSet.addLazyFigure(dataDict, cacheLabel+".png", caption,
                                      plotModule, areaLabel=label, plotargs="")
            else:
                dataDict, isSummary = qaPlotUtil.unshelveGlob(cacheLabel+"-all.png", testSet=testSet)
                dataDict['title'] = 'All Sensors'
                fig = plotModule.plot(dataDict)                
                testSet.addFigure(fig, pngFile, caption, areaLabel=label)
                del fig


            self.combineOutputs(data, dataId)
                
            

