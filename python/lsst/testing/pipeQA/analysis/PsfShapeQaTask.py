import sys, os, re
import numpy

import lsst.afw.math                as afwMath
import lsst.afw.geom                as afwGeom
import lsst.meas.algorithms         as measAlg
import lsst.pex.config              as pexConfig
import lsst.pipe.base               as pipeBase

from   .QaAnalysisTask              import QaAnalysisTask
import lsst.testing.pipeQA.figures  as qaFig
import lsst.testing.pipeQA.TestCode as testCode
import RaftCcdData                  as raftCcdData
import QaAnalysisUtils              as qaAnaUtil
import QaPlotUtils                  as qaPlotUtil



class PsfShapeQaConfig(pexConfig.Config): 
    cameras  = pexConfig.ListField(dtype = str, doc = "Cameras to run PsfShapeQaTask",
                                   default = ("lsstSim", "hscSim", "suprimecam", "cfht", "sdss", "coadd"))
    ellipMax = pexConfig.Field(dtype = float, doc = "Maximum median ellipticity", default = 0.2)
    fwhmMax  = pexConfig.Field(dtype = float, doc = "Maximum Psf Fwhm (arcsec)", default = 2.0)
    
    
class PsfShapeQaTask(QaAnalysisTask):
    ConfigClass = PsfShapeQaConfig
    _DefaultName = "psfShapeQa"

    def __init__(self, **kwargs):
        QaAnalysisTask.__init__(self, **kwargs)
        
        self.limitsEllip = [0.0, self.config.ellipMax]
        self.limitsFwhm  = [0.0, self.config.fwhmMax]
        
        self.description = """
         For each CCD, the ellipticity of stars used in the Psf model are
         plotted as a function of position in the focal plane.  The summary FPA
         figures show the median vector (offset and angle) of this ellipticity
         for each chip, as well as the effective FWHM in arcsec for the final
         Psf model.
        """

    def free(self):
        del self.theta
        del self.ellip
        del self.y
        del self.x
        del self.filter
        del self.detector
        del self.ssDict

        del self.ellipMedians
        del self.thetaMedians
        del self.fwhm

        del self.calexpDict
        
    def test(self, data, dataId):

        # get data
        self.ssDict        = data.getSourceSetBySensor(dataId)
        self.detector      = data.getDetectorBySensor(dataId)
        self.filter        = data.getFilterBySensor(dataId)
        self.calexpDict    = data.getCalexpBySensor(dataId)
        self.wcs           = data.getWcsBySensor(dataId)

        # create containers for data in the focal plane
        self.x     = raftCcdData.RaftCcdVector(self.detector)
        self.y     = raftCcdData.RaftCcdVector(self.detector)
        self.ra    = raftCcdData.RaftCcdVector(self.detector)
        self.dec   = raftCcdData.RaftCcdVector(self.detector)
        self.ellip = raftCcdData.RaftCcdVector(self.detector)
        self.theta = raftCcdData.RaftCcdVector(self.detector)

        # compute values of interest
        filter = None
        sigmaToFwhm = 2.0*numpy.sqrt(2.0*numpy.log(2.0))

        fwhmByKey = {}
        for key, ss in self.ssDict.items():
            
            if self.detector.has_key(key):
                raft = self.detector[key].getParent().getId().getName()
                ccd  = self.detector[key].getId().getName()
            else:
                continue

            fwhmByKey[key] = 0.0

            mags = []
            for s in ss:
                flux = s.getD(data.k_Psf)
                if flux > 0 and numpy.isfinite(flux) and not s.getD(data.k_ext):
                    m = -2.5*numpy.log10(flux)
                    mags.append(m)
            mag_med = numpy.median(mags)

            
            fwhmTmp = []
            for s in ss:
                ixx = s.getD(data.k_ixx)
                iyy = s.getD(data.k_iyy)
                ixy = s.getD(data.k_ixy)

                tmp = 0.25*(ixx-iyy)**2 + ixy**2
                if tmp < 0:
                    continue

                a2 = 0.5*(ixx+iyy) + numpy.sqrt(tmp)
                b2 = 0.5*(ixx+iyy) - numpy.sqrt(tmp)

                if a2 == 0 or b2/a2 < 0:
                    continue
                
                ellip = 1.0 - numpy.sqrt(b2/a2)
                theta = 0.5*numpy.arctan2(2.0*ixy, ixx-iyy)

                # vectors have no direction, so default to pointing in +ve 'y'
                # - failing to do this caused a stats bug when alignment is near pi/2
                #   both +/- pi/2 arise but are essentially the same, ... and the mean is near zero
                if theta < 0.0:
                    theta += numpy.pi
                    
                isStar = 0 if s.getD(data.k_ext) else 1

                flux = s.getD(data.k_Psf)
                mag = 99.0
                if flux > 0:
                    mag = -2.5*numpy.log10(s.getD(data.k_Psf))
                if (numpy.isfinite(ellip) and numpy.isfinite(theta) and
                    isStar and mag < mag_med and not data.isFlagged(s)):

                    self.ellip.append(raft, ccd, ellip)
                    self.theta.append(raft, ccd, theta)
                    self.x.append(raft, ccd,   s.getD(data.k_x))
                    self.y.append(raft, ccd,   s.getD(data.k_y))
                    self.ra.append(raft, ccd,  s.getD(data.k_Ra))
                    self.dec.append(raft, ccd, s.getD(data.k_Dec))
                    fwhmTmp.append( sigmaToFwhm*numpy.sqrt(0.5*(a2 + b2)))

            if len(fwhmTmp):
                fwhmByKey[key] = numpy.mean(fwhmTmp)
            else:
                fwhmByKey[key] = 0.0
                
        # create a testset and add values
        testSet = self.getTestSet(data, dataId)
        testSet.addMetadata({"Description": self.description})

        # gets the stats for each sensor and put the values in the raftccd container
        self.ellipMedians = raftCcdData.RaftCcdData(self.detector)
        self.thetaMedians = raftCcdData.RaftCcdData(self.detector)
        for raft, ccd in self.ellip.raftCcdKeys():
            ellip = self.ellip.get(raft, ccd)
            theta = self.theta.get(raft, ccd)

            if len(ellip) > 0:
                stat = afwMath.makeStatistics(ellip, afwMath.NPOINT | afwMath.MEDIAN)
                ellipMed = stat.getValue(afwMath.MEDIAN)
                stat = afwMath.makeStatistics(theta, afwMath.NPOINT | afwMath.MEDIAN)
                thetaMed = stat.getValue(afwMath.MEDIAN)
                n      = stat.getValue(afwMath.NPOINT)
            else:
                ellipMed = -1.0
                thetaMed = 0.0
                n = 0

            # add a test for acceptible psf ellipticity
            self.ellipMedians.set(raft, ccd, ellipMed)
            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            label = "median psf ellipticity "
            comment = "median psf ellipticity (nstar=%d)" % (n)
            testSet.addTest( testCode.Test(label, ellipMed, self.limitsEllip, comment, areaLabel=areaLabel) )

            # stash the angles.  We'll use them to make figures in plot()
            self.thetaMedians.set(raft, ccd, thetaMed)
            

        # And the Fwhm
        self.fwhm  = raftCcdData.RaftCcdData(self.detector)
        for key, item in self.calexpDict.items():
            if (self.detector.has_key(key) and hasattr(self.detector[key], 'getParent') and
                hasattr(self.detector[key], 'getId')):
                raft = self.detector[key].getParent().getId().getName()
                ccd  = self.detector[key].getId().getName()
            else:
                continue

            wcs = self.wcs[key]
            fwhmTmp = float(fwhmByKey[key]*wcs.pixelScale().asArcseconds())
            self.fwhm.set(raft, ccd, fwhmTmp)
            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            label = "psf fwhm (arcsec) "
            comment = "psf fwhm (arcsec)"
            testSet.addTest( testCode.Test(label, fwhmTmp, self.limitsFwhm, comment, areaLabel=areaLabel) )


    def plot(self, data, dataId, showUndefined=False):

        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)
        isFinalDataId = False
        if len(data.brokenDataIdList) == 0 or data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True

        vLen = 3000.0  # for e=1.0

        # fpa figures
        ellipBase = "medPsfEllip"
        ellipFig = qaFig.VectorFpaQaFigure(data.cameraInfo, data=None, map=None)

        fwhmBase = "psfFwhm"
        fwhmFig = qaFig.FpaQaFigure(data.cameraInfo, data=None, map=None)

        
        if self.summaryProcessing != self.summOpt['summOnly']:
            for raft, ccdDict in ellipFig.data.items():
                for ccd, value in ccdDict.items():
                    label = data.cameraInfo.getDetectorName(raft, ccd)
                    if self.ellipMedians.get(raft, ccd) is not None:
                        ellipFig.data[raft][ccd] = [self.thetaMedians.get(raft, ccd),
                                                    10*vLen*self.ellipMedians.get(raft, ccd),
                                                    self.ellipMedians.get(raft, ccd)]
                        ellipFig.map[raft][ccd] = "ell/theta=%.3f/%.0f" % (self.ellipMedians.get(raft, ccd),
                                                                           numpy.degrees(self.thetaMedians.get(raft, ccd)))
                        testSet.pickle(ellipBase+label, [ellipFig.data, ellipFig.map])
                    if self.fwhm.get(raft, ccd) is not None:
                        fwhm = self.fwhm.get(raft, ccd)
                        fwhmFig.data[raft][ccd] = fwhm
                        fwhmFig.map[raft][ccd] = "fwhm=%.2f asec" % (fwhm)
                        testSet.pickle(fwhmBase+label, [fwhmFig.data, fwhmFig.map])


        vlimMin, vlimMax = self.limitsFwhm[0], self.limitsFwhm[1]


        if (self.summaryProcessing in [self.summOpt['summOnly'], self.summOpt['delay']]) and isFinalDataId:
            
            fwhmMin =  1e10
            fwhmMax = -1e10
            
            for raft, ccdDict in ellipFig.data.items():
                for ccd, value in ccdDict.items():
                    label = data.cameraInfo.getDetectorName(raft, ccd)
                    ellipDataTmp, ellipMapTmp = testSet.unpickle(ellipBase+label, default=[None, None])
                    ellipFig.mergeValues(ellipDataTmp, ellipMapTmp)
                    fwhmDataTmp, fwhmMapTmp = testSet.unpickle(fwhmBase+label, default=[None, None])
                    fwhmFig.mergeValues(fwhmDataTmp, fwhmMapTmp)

                    fwhm = None
                    if fwhmFig.data[raft][ccd] is not None:
                        fwhm = fwhmFig.data[raft][ccd]

                    if fwhm is not None:
                        if fwhm > fwhmMax:
                            fwhmMax = fwhm
                        if fwhm < fwhmMin:
                            fwhmMin = fwhm

            
            if fwhmMin < 1e10:
                vlimMin = numpy.max([self.limitsFwhm[0], fwhmMin])
            else:
                vlimMin = self.limitsFwhm[0]
            if fwhmMax > -1e10:
                vlimMax = numpy.min([self.limitsFwhm[1], fwhmMax])
            else:
                vlimMax = self.limitsFwhm[1]

            if vlimMax < vlimMin:
                vlimMax = vlimMin + (self.limitsFwhm[1] - self.limitsFwhm[0])

            self.log.log(self.log.INFO, "plotting FPAs")
            ellipFig.makeFigure(showUndefined=showUndefined, cmap="Reds", vlimits=self.limitsEllip,
                                title="Median PSF Ellipticity", failLimits=self.limitsEllip)
            testSet.addFigure(ellipFig, ellipBase+".png", "Median PSF Ellipticity", navMap=True)


            blue = '#0000ff'
            red = '#ff0000'
            
            fwhmFig.makeFigure(showUndefined=showUndefined, cmap="jet", vlimits=[vlimMin, vlimMax],
                               title="PSF FWHM (arcsec)", cmapOver=red, failLimits=self.limitsFwhm,
                               cmapUnder=blue)
            testSet.addFigure(fwhmFig, fwhmBase + ".png", "FWHM of Psf (arcsec)", navMap=True)

        del ellipFig, fwhmFig


            

        cacheLabel = "psfEllip"


        xlo, xhi, ylo, yhi = 1.e10, -1.e10, 1.e10, -1.e10
        
        for raft,ccd in data.cameraInfo.raftCcdKeys:
            xxlo, yylo, xxhi, yyhi = data.cameraInfo.getBbox(raft, ccd)
            if xxlo < xlo: xlo = xxlo
            if xxhi > xhi: xhi = xxhi
            if yylo < ylo: ylo = yylo
            if yyhi > yhi: yhi = yyhi


        if self.summaryProcessing != self.summOpt['summOnly']:
            i = 0
            xmin, xmax = 1.0e99, -1.0e99
            for raft, ccd in self.ellip.raftCcdKeys():
                eLen = self.ellip.get(raft, ccd)

                t = self.theta.get(raft, ccd)
                dx = eLen*numpy.cos(t)
                dy = eLen*numpy.sin(t)
                x = self.x.get(raft, ccd)
                y = self.y.get(raft, ccd)

                fwhm = self.fwhm.get(raft, ccd)

                self.log.log(self.log.INFO, "plotting %s" % (ccd))


                if data.cameraInfo.name == 'coadd':
                    xmin, ymin, xmax, ymax = x.min(), y.min(), x.max(), y.max()
                    x -= xmin
                    y -= ymin
                    xxlo, yylo, xxhi, yyhi = xmin, ymin, xmax, ymax
                    xlo, ylo, xhi, yhi = xmin, ymin, xmax, ymax
                else:
                    xxlo, yylo, xxhi, yyhi = data.cameraInfo.getBbox(raft, ccd)
                limits = [xxlo, xxhi, yylo, yyhi]

                dataDict = {
                    't' : t, 'x' : x+xxlo, 'y' : y+yylo, 'dx' : dx, 'dy' : dy,
                    'color' : 'k', 'limits' : [0, xxhi-xxlo, 0, yyhi-yylo],
                    'alllimits' : [xlo, xhi, ylo, yhi],
                    'bbox' : [xxlo, xxhi, yylo, yyhi],
                    'vLen' : vLen, 'fwhm' : numpy.array([fwhm]*len(t)), 'vlim' : [vlimMin, vlimMax],
                    'summary' : False,
                    }
                label = data.cameraInfo.getDetectorName(raft, ccd)
                import PsfShapeQaPlot as plotModule
                caption = "PSF ellipticity (e=1 shown with length %.0f pix))"%(vLen)
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
            import PsfShapeQaPlot as plotModule
            caption = "PSF ellipticity " + label
            pngFile = cacheLabel + ".png"
                

            if self.lazyPlot in ['all']:
                testSet.addLazyFigure({}, cacheLabel+".png", caption,
                                      plotModule, areaLabel=label, plotargs="")
            else:
                dataDict, isSummary = qaPlotUtil.unshelveGlob(cacheLabel+"-all.png", testSet=testSet)
                if 'x' in dataDict:
                    dataDict['summary'] = True
                    dataDict['vLen'] = 5.0*vLen
                    dataDict['limits'] = [xlo, xhi, ylo, yhi]
                    fig = plotModule.plot(dataDict)                
                    testSet.addFigure(fig, pngFile, caption, areaLabel=label)
                    del fig


            self.combineOutputs(data, dataId)
