import sys, os, re
import lsst.meas.algorithms         as measAlg
import lsst.testing.pipeQA.figures  as qaFig
import numpy                        as num

import lsst.afw.math                as afwMath
import lsst.afw.coord               as afwCoord
import lsst.afw.geom                as afwGeom
import lsst.testing.pipeQA.TestCode as testCode
import lsst.testing.pipeQA.figures.QaFigureUtils as qaFigUtils

import QaAnalysis as qaAna
import RaftCcdData as raftCcdData
import QaAnalysisUtils as qaAnaUtil

import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.font_manager import FontProperties
from .AstrometricErrorQaAnalysis import AstrometricErrorQaAnalysis

class VisitToVisitAstromQaAnalysis(AstrometricErrorQaAnalysis):

    def __init__(self, maxErr, database, visits, **kwargs):
        AstrometricErrorQaAnalysis.__init__(self, maxErr, **kwargs)
        self.database      = database
        self.visits        = visits
        
        self.description = """

         For each CCD in the reference visit, its sky boundary in the
         comparison visit is queried for sources.  Differences in the
         astrometrically calibrated positions of objects are compared.

        """
    def free(self):
        AstrometricErrorQaAnalysis.free(self)
        del self.visitMatches   


    def test(self, data, dataId):
        self.matchListDictSrc = data.getMatchListBySensor(dataId, useRef='src')
        self.detector         = data.getDetectorBySensor(dataId)
        self.filter           = data.getFilterBySensor(dataId)

        if len(self.visits) == 0:
            visitList = [dataId['visit'],]
        else:
            visitList = self.visits

        # Implement for 1 visit comparison
        visit = visitList[0]
        self.visitMatches = data.getVisitMatchesBySensor(self.database, visit, dataId)

        # create containers for data we're interested in
        self.dRa   = raftCcdData.RaftCcdVector(self.detector)
        self.dDec  = raftCcdData.RaftCcdVector(self.detector)
        self.x     = raftCcdData.RaftCcdVector(self.detector)
        self.y     = raftCcdData.RaftCcdVector(self.detector)
    
        badFlags = measAlg.Flags.INTERP_CENTER | measAlg.Flags.SATUR_CENTER | measAlg.Flags.EDGE

        for key in self.matchListDictSrc.keys():
            raft = self.detector[key].getParent().getId().getName()
            ccd  = self.detector[key].getId().getName()
    
            # All detections matched to ref objects
            srcMatchList = self.matchListDictSrc[key]['matched']
            # All visit dets in footprint of original image
            visMatchList = self.visitMatches[key]

            # List of reference object ids
            srcObjIds = num.array([x[0].getId() for x in srcMatchList])
            visObjIds = num.array([x[0].getId() for x in visMatchList])
            common    = list(set(srcObjIds) & set(visObjIds))

            print key, ":", 
            if len(common) == 0:
                print "No overlap, Using pre-to-post PT1.2 objectID mapping...", 
                    
                # Try objectID hack from PT1.2 to post-PT1.2:
                visObjIds *= 2
                isStar     = (num.array([x[0].getFlagForDetection() for x in visMatchList]) & measAlg.Flags.STAR) > 0
                visObjIds += isStar
                common     = list(set(srcObjIds) & set(visObjIds))
            print "Found %d matches" % (len(common))
                    
            # Iterate over all object ids
            for i in range(len(common)):
                srcObjId  = common[i]
                idxS = num.where(srcObjIds == srcObjId)[0]
                idxV = num.where(visObjIds == srcObjId)[0] 

                # only take 1-to-1 matches
                if len(idxS) != 1 or len(idxV) != 1:
                    continue
    
                sref1 = srcMatchList[idxS[0]][0]
                srcv1 = srcMatchList[idxS[0]][1]
    
                sref2 = visMatchList[idxV[0]][0]
                srcv2 = visMatchList[idxV[0]][1]

                #coord1 = afwCoord.Coord(afwGeom.Point2D(srcv1.getRa(), srcv1.getDec()))  # coords originally in degrees
                #coord2 = afwCoord.Coord(afwGeom.Point2D(srcv2.getRa(), srcv2.getDec()))  # coords originally in degrees
                #angle  = coord1.angularSeparation(coord2)

                # Measurment flags; note no star/gal separation yet
                flags1 = srcv1.getFlagForDetection()
                flags2 = srcv2.getFlagForDetection()
                # Get star/gal info here
                isStar = sref1.getFlagForDetection() & measAlg.Flags.STAR

                if (not flags1 & badFlags) and (not flags2 & badFlags) and isStar:
                    ra1, dec1, ra2, dec2 = [x * num.pi / 180.0 for x in (srcv1.getRa(), srcv1.getDec(), srcv2.getRa(), srcv2.getDec())]
                    dRa  = (ra1 - ra2) * abs(num.cos(dec1))
                    dDec = (dec1 - dec2)

                    self.dRa.append(raft, ccd, dRa)
                    self.dDec.append(raft, ccd, dDec)
                    self.x.append(raft, ccd, srcv1.getXAstrom())
                    self.y.append(raft, ccd, srcv1.getYAstrom())

        testSet = self.getTestSet(data, dataId)
        testSet.addMetadata({"Description": self.description})

        self.medErrArcsec = raftCcdData.RaftCcdData(self.detector)
        self.medThetaRad  = raftCcdData.RaftCcdData(self.detector)
    
        for raft, ccd in self.dRa.raftCcdKeys():
            dRa  = self.dRa.get(raft, ccd)
            dDec = self.dDec.get(raft, ccd)
                
            errArcsec = 206265.0*num.sqrt(dRa**2 + dDec**2)
            thetaRad  = num.arctan2(dDec, dRa)
    
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



