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

import sys, os, re
import numpy

import time

import lsst.afw.geom                as afwGeom
import lsst.afw.math                as afwMath
import lsst.meas.algorithms         as measAlg
import lsst.pex.config              as pexConfig
import lsst.pipe.base               as pipeBase

from   .QaAnalysisTask              import QaAnalysisTask
import lsst.testing.pipeQA.figures  as qaFig
import lsst.testing.pipeQA.TestCode as testCode
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtil
import RaftCcdData                  as raftCcdData
import QaAnalysisUtils              as qaAnaUtil

import QaPlotUtils                  as qaPlotUtil



class AstrometricErrorQaConfig(pexConfig.Config):
    cameras = pexConfig.ListField(dtype = str,
                                  doc = "Cameras to run AstrometricErrorQaTask",
                                  default = ("lsstSim", "hscSim", "suprimecam", "cfht", "sdss", "coadd"))
    maxErr  = pexConfig.Field(dtype = float,
                              doc = "Maximum astrometric error (in arcseconds)",
                              default = 0.5)

class AstrometricErrorQaTask(QaAnalysisTask):
    ConfigClass = AstrometricErrorQaConfig
    _DefaultName = "astrometricErrorQa"

    def __init__(self, **kwargs):
        QaAnalysisTask.__init__(self, **kwargs)
        self.limits = [0.0, self.config.maxErr]

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

    def test(self, data, dataId):

        # get data
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')
        self.detector         = data.getDetectorBySensor(dataId)
        self.filter           = data.getFilterBySensor(dataId)
        

        # compute the mean ra, dec for each source cluster
        self.dRa  = raftCcdData.RaftCcdVector(self.detector)
        self.dDec = raftCcdData.RaftCcdVector(self.detector)
        self.x    = raftCcdData.RaftCcdVector(self.detector)
        self.y    = raftCcdData.RaftCcdVector(self.detector)


        filter = None
        key = None
        for key in self.matchListDictSrc.keys():
            raft = self.detector[key].getParent().getId().getName()
            ccd  = self.detector[key].getId().getName()
            filter = self.filter[key].getName()

            matchList = self.matchListDictSrc[key]['matched']
            for m in matchList:
                sref, s, dist = m
                ra, dec, raRef, decRef = \
                    [numpy.radians(x) for x in [s.getD(data.k_Ra), s.getD(data.k_Dec),
                                                sref.getD(data.k_rRa), sref.getD(data.k_rDec)]]
                

                
                dDec = decRef - dec
                dRa  = (raRef - ra)*abs(numpy.cos(decRef))
                
                if not data.isFlagged(s):
                    self.dRa.append(raft, ccd, dRa)
                    self.dDec.append(raft, ccd, dDec)
                    self.x.append(raft, ccd, s.getD(data.k_x))
                    self.y.append(raft, ccd, s.getD(data.k_y))
                    
                    
        testSet = self.getTestSet(data, dataId)
        testSet.addMetadata({"Description": self.description})
        
        self.medErrArcsec = raftCcdData.RaftCcdData(self.detector)
        self.medThetaRad  = raftCcdData.RaftCcdData(self.detector)

        for raft,  ccd in self.dRa.raftCcdKeys():
            dRa  = self.dRa.get(raft, ccd).copy()
            dDec = self.dDec.get(raft, ccd).copy()

            if len(dRa) > 0:
                dRaMed = numpy.median(dRa)
                dDecMed = numpy.median(dDec)
            else:
                dRaMed = 0.0
                dDecMed = 0.0

            sysErr = numpy.sqrt(dRaMed**2 + dDecMed**2)*afwGeom.radians
            sysErrArcsec = sysErr.asArcseconds()
            sysThetaRad  = numpy.arctan2(dDecMed, dRaMed)
            
            dRa  -= dRaMed
            dDec -= dDecMed

            rmsErr = numpy.sqrt(dRa**2 + dDec**2)
            rmsThetaRad  = numpy.arctan2(dDec, dRa)

            if len(rmsErr) > 0:
                stat  = afwMath.makeStatistics(rmsErr, afwMath.NPOINT | afwMath.MEDIAN)
                medRmsErr = stat.getValue(afwMath.MEDIAN)
                stat  = afwMath.makeStatistics(rmsThetaRad, afwMath.NPOINT | afwMath.MEDIAN)
                medRmsThetaRad = stat.getValue(afwMath.MEDIAN)
                n = stat.getValue(afwMath.NPOINT)
            else:
                medRmsErr = -1.0
                medRmsThetaRad = 0.0
                n = 0
                
            medRmsErr = medRmsErr*afwGeom.radians
            
            self.medErrArcsec.set(raft, ccd, sysErrArcsec)
            self.medThetaRad.set(raft, ccd, sysThetaRad)
            
            areaLabel = data.cameraInfo.getDetectorName(raft, ccd)
            label = "median systematic astrometry error "
            comment = "median sqrt(dRa^2+dDec^2) (arcsec, nstar=%d)" % (n)
            test = testCode.Test(label, sysErrArcsec, self.limits, comment, areaLabel=areaLabel)
            testSet.addTest(test)

            label = "median random astrometry error "
            comment = "median sqrt((dRa-dRaMed)^2+(dDec-dDecMed)^2) (arcsec, nstar=%d)" % (n)
            test = testCode.Test(label, medRmsErr.asArcseconds(), self.limits, comment, areaLabel=areaLabel)
            testSet.addTest(test)
            
        
    def plot(self, data, dataId, showUndefined=False):

        testSet = self.getTestSet(data, dataId)
        testSet.setUseCache(self.useCache)


        isFinalDataId = False
        if len(data.brokenDataIdList) == 0 or data.brokenDataIdList[-1] == dataId:
            isFinalDataId = True


        medAstBase = "medAstError"
        astFig = qaFig.VectorFpaQaFigure(data.cameraInfo, data=None, map=None)

        vLen = 5000 # length in pixels for 1 arcsec error vector
        if self.summaryProcessing != self.summOpt['summOnly']:
            for raft, ccdDict in astFig.data.items():
                for ccd, value in ccdDict.items():
                    if self.medErrArcsec.get(raft, ccd) is not None:
                        astErrArcsec = self.medErrArcsec.get(raft, ccd)
                        thetaRad = self.medThetaRad.get(raft, ccd)
                        astFig.data[raft][ccd] = [thetaRad, vLen*astErrArcsec, astErrArcsec]
                        astFig.map[raft][ccd] = "\"/theta=%.2f/%.0f" % (astErrArcsec, numpy.degrees(thetaRad))

                        label = data.cameraInfo.getDetectorName(raft, ccd)
                        testSet.pickle(medAstBase + label, [astFig.data, astFig.map])

        # only make the FPA figure if this is a 'summary' run, or the final CCD
        if (self.summaryProcessing in [self.summOpt['summOnly'], self.summOpt['delay']]) and isFinalDataId:

            for raft, ccdDict in astFig.data.items():
                for ccd, value in ccdDict.items():
                    label = data.cameraInfo.getDetectorName(raft, ccd)
                    medAstDataTmp, medAstMapTmp = testSet.unpickle(medAstBase+label, default=[None, None])
                    astFig.mergeValues(medAstDataTmp, medAstMapTmp)

            self.log.log(self.log.INFO, "plotting FPAs")
            astFig.makeFigure(showUndefined=showUndefined, cmap="Reds", vlimits=[0.0, 2.0*self.limits[1]],
                              title="Median astrometric error", cmapOver='#ff0000', failLimits=self.limits,
                              cmapUnder="#ff0000")
            testSet.addFigure(astFig, medAstBase+".png", "Median astrometric error",  navMap=True)


        del astFig

        prePlot = False

        cacheLabel = "astromError"


        if self.summaryProcessing != self.summOpt['summOnly']:


            xlo, xhi, ylo, yhi = 1.e10, -1.e10, 1.e10, -1.e10
            for raft,ccd in data.cameraInfo.raftCcdKeys:
                if data.cameraInfo.name == 'coadd':
                    xtmp, ytmp = self.x.get(raft, ccd), self.y.get(raft, ccd)
                    xxlo, yylo, xxhi, yyhi = xtmp.min(), ytmp.min(), xtmp.max(), ytmp.max()
                else:
                    xxlo, yylo, xxhi, yyhi = data.cameraInfo.getBbox(raft, ccd)
                if xxlo < xlo: xlo = xxlo
                if xxhi > xhi: xhi = xxhi
                if yylo < ylo: ylo = yylo
                if yyhi > yhi: yhi = yyhi


            for raft, ccd in self.dRa.raftCcdKeys():
                ra = self.dRa.get(raft, ccd)
                dec = self.dDec.get(raft, ccd)
                dAngle = numpy.sqrt(ra**2 + dec**2)
                eLen = 3600.0*numpy.degrees(dAngle)
                t = numpy.arctan2(dec, ra)

                dx = eLen*numpy.cos(t)
                w = numpy.where(numpy.abs(dx) < 10.0)

                dx = dx[w]
                dy = (eLen*numpy.sin(t))[w]
                x = (self.x.get(raft, ccd))[w]
                y = (self.y.get(raft, ccd))[w]                

                self.log.log(self.log.INFO, "plotting %s" % (ccd))


                if data.cameraInfo.name == 'coadd':
                    xmin, ymin, xmax, ymax = x.min(), y.min(), x.max(), y.max()
                    x -= xmin
                    y -= ymin
                    xxlo, yylo, xxhi, yyhi = xmin, ymin, xmax, ymax
                else:
                    xxlo, yylo, xxhi, yyhi = data.cameraInfo.getBbox(raft, ccd)


                import AstrometricErrorQaPlot as plotModule
                label = data.cameraInfo.getDetectorName(raft, ccd)
                dataDict = {'x': x, 'y':y, 'dx':dx, 'dy':dy,
                            'limits' : [0, xxhi-xxlo, 0, yyhi-yylo],
                            'bbox' : [0, xxhi-xxlo, 0, yyhi-yylo],
                            'gridVectors':False }
                caption = "Astrometric error" + label
                pngFile = cacheLabel+".png"


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
            
            import AstrometricErrorQaPlot as plotModule
            label = 'all'
            caption = "Astrometric error " + label
            pngFile = "astromError.png"
            
            if self.lazyPlot in ['all']:
                testSet.addLazyFigure({}, cacheLabel+".png", caption,
                                      plotModule, areaLabel=label, plotargs="")
            else:
                dataDict, isSummary = qaPlotUtil.unshelveGlob(cacheLabel+"-all.png", testSet=testSet)
                dataDict['gridVectors'] = True
                if 'x' in dataDict:
                    fig = plotModule.plot(dataDict)                
                    testSet.addFigure(fig, pngFile, caption, areaLabel=label)
                    del fig

            self.combineOutputs(data, dataId)
                    
            

