import re
import numpy as num
import time

import lsst.meas.algorithms         as measAlg
import lsst.testing.pipeQA.TestCode as testCode
import QaAnalysis as qaAna
import lsst.testing.pipeQA.figures as qaFig
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtils
import QaAnalysisUtils as qaAnaUtil
import RaftCcdData as raftCcdData

import lsst.testing.pipeQA.source as pqaSource

from matplotlib.font_manager import FontProperties
from matplotlib.patches import Ellipse

class ZeropointFitQa(qaAna.QaAnalysis):
    def __init__(self, medOffsetMin, medOffsetMax, figsize=(5.0,5.0), **kwargs):
        qaAna.QaAnalysis.__init__(self, **kwargs)
        self.figsize = figsize
        self.limits = [medOffsetMin, medOffsetMax]

        self.sCatDummy = pqaSource.Catalog()
        self.srefCatDummy = pqaSource.RefCatalog()
        
        self.description = """
         For each CCD, the central panel shows the instrumental magnitude of
         matched stars and galaxies, plotted as a function of the catalog
         magnitude of the reference objects they were matched to.  The fitted
         zeropoint is shown as the dashed line.  The bottom panel shows a
         histogram of the "matched" and "orphan" source as a function of
         instrumental magnitude, while the left panel shows a histogram of the
         "matched" and "unmatched" entries in the reference catalog.  The top
         panel shows the scatter of the matched stars and galaxies around the
         zeropoint fit (similar to pipeQa.PhotCompareQaAnalysis.psf-cat) with
         the median offset of star from the zeropoint indicated with the dotted
         line.  The summary FPA figures show the median offset of stars from
         the zeropoint across the focal plane, as well as the fitted zeropoint.
        """

    def free(self):
        del self.detector
        del self.filter
        del self.calib
        del self.matchListDictSrc

        del self.orphan
        del self.matchedStar
        del self.undetectedStar
        del self.matchedGalaxy
        del self.undetectedGalaxy
        
        del self.zeroPoint
        del self.medOffset
        
    def test(self, data, dataId, fluxType = "psf"):

        testSet = self.getTestSet(data, dataId)
        testSet.addMetadata({"Description": self.description})

        self.fluxType         = fluxType
        self.detector         = data.getDetectorBySensor(dataId)
        self.filter           = data.getFilterBySensor(dataId)
        self.calib            = data.getCalibBySensor(dataId)
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')

        # Ignore blends in this analysis
        self.orphan           = raftCcdData.RaftCcdVector(self.detector)
        self.matchedStar      = raftCcdData.RaftCcdData(self.detector)
        self.undetectedStar   = raftCcdData.RaftCcdVector(self.detector)
        self.matchedGalaxy    = raftCcdData.RaftCcdData(self.detector)
        self.undetectedGalaxy = raftCcdData.RaftCcdVector(self.detector)

        # Results
        self.zeroPoint        = raftCcdData.RaftCcdData(self.detector)
        self.medOffset        = raftCcdData.RaftCcdData(self.detector)

        badFlags = pqaSource.INTERP_CENTER | pqaSource.SATUR_CENTER | pqaSource.EDGE

        for key in self.detector.keys():
            raftId     = self.detector[key].getParent().getId().getName()
            ccdId      = self.detector[key].getId().getName()
            filterName = self.filter[key].getName()


            self.matchedStar.set(raftId, ccdId, {"Refmag": num.array([]),
                                                 "Imgmag": num.array([]),
                                                 "Imgerr": num.array([]),})
            self.matchedGalaxy.set(raftId, ccdId, {"Refmag": num.array([]),
                                                   "Imgmag": num.array([]),
                                                   "Imgerr": num.array([])})
            self.undetectedStar.set(raftId, ccdId, num.array([]))
            self.undetectedGalaxy.set(raftId, ccdId, num.array([]))
            self.orphan.set(raftId, ccdId, num.array([]))


            fmag0 = self.calib[key].getFluxMag0()[0]
            if fmag0 <= 0.0:
                self.zeroPoint.set(raftId, ccdId, NaN)
                continue
            zpt = -2.5*num.log10(fmag0)
            self.zeroPoint.set(raftId, ccdId, zpt)

            if self.matchListDictSrc.has_key(key):
                # Matched
                mdict    = self.matchListDictSrc[key]['matched']
                stars    = []
                galaxies = []

                for m in mdict:
                    sref, s, dist = m
                    if fluxType == "psf":
                        fref  = sref.getD(self.srefCatDummy.PsfFluxKey)
                        f     = s.getD(self.sCatDummy.PsfFluxKey)
                        ferr  = s.getD(self.sCatDummy.PsfFluxErrKey)
                    else:
                        fref  = sref.getD(self.srefCatDummy.PsfFluxKey)
                        f     = s.getD(self.sCatDummy.ApFluxKey)
                        ferr  = s.getD(self.sCatDummy.ApFluxErrKey)

                    # un-calibrate the magnitudes
                    f *= fmag0
                    
                    intcen = s.getD(self.sCatDummy.FlagPixInterpCenKey)
                    satcen = s.getD(self.sCatDummy.FlagPixSaturCenKey)
                    edge   = s.getD(self.sCatDummy.FlagPixEdgeKey)

                    if (fref > 0.0 and f > 0.0 and not (intcen or satcen or edge)):
                        mrefmag  = -2.5*num.log10(fref)
                        mimgmag  = -2.5*num.log10(f)
                        mimgmerr =  2.5 / num.log(10.0) * ferr / f
    
                        star = 0 if s.getD(self.sCatDummy.ExtendednessKey) else 1
                        
                        if num.isfinite(mrefmag) and num.isfinite(mimgmag):
                            if star:
                                stars.append((mrefmag, mimgmag, mimgmerr))
                            else:
                                galaxies.append((mrefmag, mimgmag, mimgmerr))
                self.matchedStar.set(raftId, ccdId, {"Refmag": num.array([x[0] for x in stars]),
                                                     "Imgmag": num.array([x[1] for x in stars]),
                                                     "Imgerr": num.array([x[2] for x in stars])})
                self.matchedGalaxy.set(raftId, ccdId, {"Refmag": num.array([x[0] for x in galaxies]),
                                                       "Imgmag": num.array([x[1] for x in galaxies]),
                                                       "Imgerr": num.array([x[2] for x in galaxies])})
            
                # Non-detections
                undetectedStars = []
                undetectedGalaxies = []
                for nondet in self.matchListDictSrc[key]['undetected']:
                    mag = nondet.getMag(filterName)
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
                        f = orphan.getD(self.sCatDummy.PsfFluxKey)
                    else:
                        f = orphan.getD(self.sCatDummy.ApFluxKey)
                    if f > 0.0:
                        # un-calibrate the magnitudes
                        f *= fmag0
                        
                        orphans.append(-2.5 * num.log10(f))
                        
                self.orphan.set(raftId, ccdId, num.array(orphans))

                # Metrics
                offset      = num.array(self.matchedStar.get(raftId, ccdId)["Imgmag"]) # make a copy
                offset     -= self.zeroPoint.get(raftId, ccdId)
                offset     -= self.matchedStar.get(raftId, ccdId)["Refmag"]
                med         = num.median(offset) 
                self.medOffset.set(raftId, ccdId, med)
                
                areaLabel = data.cameraInfo.getDetectorName(raftId, ccdId)
                label = "median offset from zeropoint"
                comment = "Median offset of calibrated stellar magnitude to zeropoint fit"
                test = testCode.Test(label, med, self.limits, comment, areaLabel=areaLabel)
                testSet.addTest(test)

                label = "median zeropoint"
                comment = "Median zeropoint measured for sensor"
                test = testCode.Test(label, zpt, [None, 0], comment, areaLabel=areaLabel)
                testSet.addTest(test)
                
            
    def plot(self, data, dataId, showUndefined=False):

        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)
        isFinalDataId = False
        if len(data.brokenDataIdList) == 0 or data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True

        # fpa figure
        zpts = []
        zptBase = "zeropoint"
        zptData, zptMap = testSet.unpickle(zptBase, default=[None, None])
        zptFig = qaFig.FpaQaFigure(data.cameraInfo, data=zptData, map=zptMap)

        offsetBase = "medZeropointOffset"
        offsetData, offsetMap = testSet.unpickle(offsetBase, default=[None, None])
        offsetFig = qaFig.FpaQaFigure(data.cameraInfo, data=offsetData, map=offsetMap)

        for raft, ccdDict in zptFig.data.items():
            for ccd, value in ccdDict.items():
                if not self.zeroPoint.get(raft, ccd) is None:
                    zpt = self.zeroPoint.get(raft, ccd)
                    zpts.append(zpt)
                    zptFig.data[raft][ccd] = zpt
                    zptFig.map[raft][ccd] = 'zpt=%.2f' % (zpt)

                    offset = self.medOffset.get(raft, ccd)
                    offsetFig.data[raft][ccd] = offset
                    offsetFig.map[raft][ccd] = 'offset=%.2f' % (offset)
                else:
                    if not zptFig.data[raft][ccd] is None:
                        zpts.append(zptFig.data[raft][ccd])
                    
                    
        testSet.pickle(zptBase, [zptFig.data, zptFig.map])
        testSet.pickle(offsetBase, [offsetFig.data, offsetFig.map])
        
        if not self.delaySummary or isFinalDataId:
            print "plotting FPAs"
            
            blue = '#0000ff'
            red = '#ff0000'
            zptFig.makeFigure(showUndefined=showUndefined, cmap="jet",
                              vlimits=[num.min(zpts)-0.05, num.max(zpts)+0.05],
                              title="Zeropoint", cmapOver=red, cmapUnder=blue)
            testSet.addFigure(zptFig, zptBase+".png", "Photometric zeropoint", navMap=True)
            del zptFig
        
            offsetFig.makeFigure(showUndefined=showUndefined, cmap="jet", vlimits=self.limits,
                                 title="Med offset from Zpt Fit", cmapOver=red, failLimits=self.limits,
                                 cmapUnder=blue)
            testSet.addFigure(offsetFig, offsetBase + ".png", "Median offset from photometric zeropoint", 
                              navMap=True)
            del offsetFig
        else:
            del zptFig, offsetFig


        cacheLabel = "zeropointFit"
        shelfData = {}
        
        # Each CCD
        for raft, ccd in self.zeroPoint.raftCcdKeys():
            zeropt = self.zeroPoint.get(raft, ccd)
            if zeropt == 0.0:
                continue
            
            print "Plotting", ccd


            # Plot all matched galaxies
            mrefGmag  = self.matchedGalaxy.get(raft, ccd)["Refmag"]
            mimgGmag  = self.matchedGalaxy.get(raft, ccd)["Imgmag"]
            mimgGmerr = self.matchedGalaxy.get(raft, ccd)["Imgerr"]

            # Plot all matched stars
            mrefSmag  = self.matchedStar.get(raft, ccd)["Refmag"]
            mimgSmag  = self.matchedStar.get(raft, ccd)["Imgmag"]
            mimgSmerr = self.matchedStar.get(raft, ccd)["Imgerr"]


            urefmag     = num.concatenate((self.undetectedStar.get(raft, ccd),
                                           self.undetectedGalaxy.get(raft, ccd)))
            uimgmag     = self.orphan.get(raft, ccd)

            label = data.cameraInfo.getDetectorName(raft, ccd)
            fig = self.standardFigure(mrefGmag, mimgGmag, mimgGmerr,
                                      mrefSmag, mimgSmag, mimgSmerr,
                                      urefmag, uimgmag,
                                      zeropt, label)

            testSet.addFigure(fig, "zeropointFit.png", "Zeropoint fit "+label, areaLabel=label)


            shelfData[ccd] = [mrefGmag, mimgGmag, mimgGmerr,
                              mrefSmag, mimgSmag, mimgSmerr,
                              urefmag, uimgmag, zeropt]
            

        if self.useCache:
            testSet.shelve(cacheLabel, shelfData)

        if not self.delaySummary or isFinalDataId:
            print "plotting Summary figure"

            # unstash the values
            if self.useCache:
                shelfData = testSet.unshelve(cacheLabel)


            mrefGmagAll   = num.array([])
            mimgGmagAll   = num.array([])
            mimgGmerrAll  = num.array([])
            mrefSmagAll   = num.array([])
            mimgSmagAll   = num.array([])
            mimgSmerrAll  = num.array([])
            urefmagAll    = num.array([])
            uimgmagAll    = num.array([])
            #zeroG         = num.array([])
            #zeroS         = num.array([])

            #n = 0
            for k,v in shelfData.items():
                mrefGmag, mimgGmag, mimgGmerr, mrefSmag, mimgSmag, mimgSmerr, urefmag, uimgmag, zeropt = v
                mrefGmagAll   = num.append(mrefGmagAll   ,mrefGmag   )
                mimgGmagAll   = num.append(mimgGmagAll   ,mimgGmag - zeropt  )
                mimgGmerrAll  = num.append(mimgGmerrAll  ,mimgGmerr  )
                mrefSmagAll   = num.append(mrefSmagAll   ,mrefSmag   )
                mimgSmagAll   = num.append(mimgSmagAll   ,mimgSmag - zeropt  )
                mimgSmerrAll  = num.append(mimgSmerrAll  ,mimgSmerr  )
                urefmagAll    = num.append(urefmagAll    ,urefmag    )
                uimgmagAll    = num.append(uimgmagAll    ,uimgmag - zeropt   )
                #zeroG         = num.append(zeroG,  num.array([zeropt]*len(mrefGmag)))
                #zeroS         = num.append(zeroS,  num.array([zeropt]*len(mrefSmag)))
                #zeroMean += zeropt
                #n+=1
            #zeroMean /= n

            allFig = self.standardFigure(mrefGmagAll, mimgGmagAll, mimgGmerrAll,
                                         mrefSmagAll, mimgSmagAll, mimgSmerrAll,
                                         urefmagAll, uimgmagAll,
                                         0.0, "All Sensors")
            
            label = "all"
            testSet.addFigure(allFig, "zeropointFit.png", "zeropoint fit "+label, areaLabel=label)
            del allFig
            

    def standardFigure(self,
                       mrefGmag, mimgGmag, mimgGmerr,
                       mrefSmag, mimgSmag, mimgSmerr,
                       urefmag, uimgmag,
                       zeropt, title):


        # Just to get the histogram results
        fig = qaFig.QaFigure(size=self.figsize)

        legLines  = []
        legLabels = []

        axis = fig.fig.add_axes([0.225, 0.225, 0.675, 0.550])


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
        ax2  = fig.fig.add_axes([0.1,   0.225, 0.125, 0.550], sharey=axis)
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


        if not nu is None and num.max(nu) > maxN2: maxN2 = num.max(nu)

        maxN3 = 999

        nm, bm, pm = None, None, None

        # Unmatched & matched stellar objects
        ax3  = fig.fig.add_axes([0.225, 0.1,   0.675, 0.125], sharex=axis)
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

        if not nm is None and num.max(nm) > maxN3: maxN3 = num.max(nm)

        if len(uimgmag) > 0:
            try:
                ax3.hist(uimgmag, bins=num.arange(xmin, xmax, 0.25),
                         log = True, facecolor = 'r', alpha = 0.5, zorder = 1)
            except:
                pass
        if abs(zeropt) < 1.0e-5:
            ax3.set_xlabel('Image calibrated %s mag' % (self.fluxType), fontsize = 10)
        else:
            ax3.set_xlabel('Image instrumental %s mag' % (self.fluxType), fontsize = 10)
        ax3.set_ylabel('N', rotation = 180, fontsize = 10)

        # Mag - Zpt
        ax4  = fig.fig.add_axes([0.225, 0.775, 0.675, 0.125], sharex=axis)
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
        qaFigUtils.qaSetp(axis.get_xticklabels()+axis.get_yticklabels(), visible=False)
        qaFigUtils.qaSetp(ax2.get_xticklabels()+ax2.get_yticklabels(), fontsize = 8)
        qaFigUtils.qaSetp(ax2.get_xticklabels(), rotation=90)
        qaFigUtils.qaSetp(ax3.get_xticklabels()+ax3.get_yticklabels(), fontsize = 8)
        qaFigUtils.qaSetp(ax4.get_xticklabels(), visible=False)
        qaFigUtils.qaSetp(ax4.get_yticklabels(), fontsize = 8)

        fig.fig.legend(legLines, legLabels,
                       numpoints=1, prop=FontProperties(size='x-small'), loc = 'center right')

        fig.fig.suptitle('%s' % (title), fontsize = 12)

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

