import os, re
import math
import copy

import numpy

import lsst.afw.image   as afwImage
import lsst.afw.coord   as afwCoord
import lsst.pex.logging  as pexLog
import lsst.pex.policy  as pexPolicy
import lsst.meas.astrom as measAstrom
import lsst.meas.algorithms.utils as maUtils

class QaDataUtils(object):

    def __init__(self):

        # specify the labels we need
        self.names = [
            "Ra",          
            "Dec",         
            "XAstrom",     
            "YAstrom",     
            "PsfFlux",     
            "PsfFluxErr",  
            "ApFlux",      
            "ApFluxErr",   
            "ModelFlux",   
            "ModelFluxErr",
            "InstFlux",    
            "InstFluxErr", 
            "Ixx",         
            "Iyy",         
            "Ixy",
            "deblend_nchild",
            "FlagBadCentroid", 
            "FlagPixSaturCen", 
            "FlagPixInterpCen",
            "FlagPixEdge",     
            "FlagNegative",    
            "Extendedness",
            ]
        
        self.types = {
            "Ra":           "D",
            "Dec":          "D",
            "XAstrom":      "D",
            "YAstrom":      "D",
            "PsfFlux":      "D",
            "PsfFluxErr":   "D",
            "ApFlux":       "D",
            "ApFluxErr":    "D",
            "ModelFlux":    "D",
            "ModelFluxErr": "D",
            "InstFlux":     "D",
            "InstFluxErr":  "D",
            "Ixx":          "D",
            "Iyy":          "D",
            "Ixy":          "D",
            "deblend_nchild": "I",
            
            # source
            "FlagBadCentroid":  "I",
            "FlagPixSaturCen":  "I",
            "FlagPixInterpCen": "I",
            "FlagPixEdge":      "I",
            "FlagNegative":     "I",

            # icsource
            #["FlagBadCentroid":  ,
            #["FlagPixSaturCen":  ,
            #["FlagPixInterpCen": ,
            #["FlagPixEdge":      ,
            #["FlagNegative":     ,

            "Extendedness":  "D",
            }



        srcNameList = self.getSourceSetNameList()
        mappedNames = zip(*srcNameList)[0]
        missing = []
        for n in self.names:
            if not n in mappedNames:
                missing.append(n)
                
        if len(missing) > 0:
            raise RuntimeError, "Required values not provided by getSourceNameList():\n%s" % ("\n".join(missing))
        



    def getSourceSetNameList(self):
        """Associate Source accessor names to database columns in a list of pairs. """
        return zip(self.names, self.names)


    def getSourceSetAccessors(self):
        """Get a list of all accessor names for Source objects. """
        return zip(*self.getSourceSetNameList())[0]

    def getSourceSetAccessorsAndTypes(self):
        types = []
        for n in self.names:
            types.append(self.types[n])
        return copy.copy(self.names), copy.copy(types)


    def getSourceSetDbNames(self, replacementDict):
        """Get a list of all database names for Source objects. """
        rawList = list(zip(*self.getSourceSetNameList())[1])
        for k,v in replacementDict.items():
            matches = [i for i,x in enumerate(rawList) if x == k]
            arg = matches[0]
            rawList[arg] = v
        return rawList



    def getCalexpNameLookup(self):
        """Associate calexp_md names to database names."""

        nameLookup = {
            'FILTER':           'filterName'           ,
            'RA':                   'ra'               ,
            'DEC':                 'decl'              ,
            'CRPIX1':               'crpix1'           ,
            'CRPIX2':               'crpix2'           ,
            'CRVAL1':               'crval1'           ,
            'CRVAL2':               'crval2'           ,
            'CD1_1':                'cd1_1'            ,
            'CD1_2':                'cd1_2'            ,
            'CD2_1':                'cd2_1'            ,
            'CD2_2':                'cd2_2'            ,
            'FLUXMAG0':             'fluxMag0'         ,
            'FLUXMAG0ERR':        'fluxMag0Sigma'      ,
            'SEEING':                 'fwhm'           ,
            'EXPTIME':              'exptime'          ,
            }

        return nameLookup

                                                                                              
    def getSceNameList(self, dataIdNames, replacements={}):
        """Associate SourceCcdExposure names to database columns in a list of pairs. """
        raise NotImplementedError, "This method should be defined in derived class."


    def getSceDbNames(self, dataIdNames, replacements={}):
        """Get SourceCcdExposure database column names."""
        return zip(*self.getSceNameList(dataIdNames, replacements))[1]




    def getCalibObjects(self, butler, filterName, dataId):
        """
        A version of getCalibObjectsImpl that isn't a class method, for use by other code

        @param useOutputSrc             # use fluxes from the "src", not "icSrc"
        """

        #################################
        # color calibration
        #
        # this stuff stolen from pipette.
        # - numbers from suprimecam.paf, method in Calibrate.py
        #####################
        primaryLookup = {
            #"B" : "g",
            #"V" : "g",
            #"R" : "r",
            #"I" : "i",
            "g" : "g",
            "g" : "g",
            "r" : "r",
            "i" : "i",
            }
        # The below colour terms are from the last page of
        # http://www.naoj.org/staff/nakata/suprime/illustration/colorterm_report_ver3.pdf
        secondaryLookup = {
            "g" : ["r", [-0.00569, -0.0427]],
            "r" : ["g", [0.00261, 0.0304]],
            "i" : ["r", [0.00586, 0.0827, -0.0118]],
            "z" : ['i', [0.000329, 0.0608, 0.0219]],
            "y" : ['i', [0.000329, 0.0608, 0.0219]],
            }
        # figure out which filters we need for the color term
        primaryFilter   = primaryLookup[filterName]      # Primary band for correction
        secondaryFilter = secondaryLookup[primaryFilter][0]  # Secondary band for correction


        log = pexLog.Log.getDefaultLog()
        psources  = butler.get('icSrc', dataId)
        pmatches  = butler.get('icMatch', dataId)
        calexp_md = butler.get('calexp_md', dataId)

        wcs = afwImage.makeWcs(calexp_md)
        W, H = calexp_md.get("NAXIS1"), calexp_md.get("NAXIS2")
        calib = afwImage.Calib(calexp_md)

        matches = pmatches.getSourceMatches()
        setMatchListBlobsNone(matches)
        sources = psources.getSources() ## fails
        setSourceSetBlobsNone(sources)

        anid = pmatches.getSourceMatchMetadata().getInt('ANINDID')

        del psources; del pmatches; del calexp_md # cleanup

        useOutputSrc = False
        if useOutputSrc:
            srcs = butler.get('src', dataId).getSources()
            setSourceSetBlobsNone(srcs)

            import lsst.afw.detection as afwDetect
            pmMatch = afwDetect.matchXy(sources, srcs, 1.0, True)
            for icSrc, src, d in pmMatch:
                icSrc.setPsfFlux(src.getPsfFlux())

        # ref sources
        xc, yc = 0.5*W, 0.5*H
        radec = wcs.pixelToSky(xc, yc)
        ra = radec.getLongitude(afwCoord.DEGREES)
        dec = radec.getLatitude(afwCoord.DEGREES)
        radius = wcs.pixelScale().asDegrees()*math.hypot(xc, yc)*1.1

        pol = pexPolicy.Policy()
        pol.set('matchThreshold', 30)
        solver = measAstrom.createSolver(pol, log)
        idName = 'id'

        X = solver.getCatalogue(ra, dec, radius, primaryFilter, idName, anid)
        refsources = X.refsources
        setSourceSetBlobsNone(refsources)

        inds = X.inds

        referrs, stargal = None, None
        cols = solver.getTagAlongColumns(anid)
        colnames = [c.name for c in cols]

        col = 'starnotgal'
        if col in colnames:
            stargal1 = solver.getTagAlongBool(anid, col, inds)
            stargal = []
            for i in range(len(stargal1)):
                stargal.append(stargal1[i])

        fdict = maUtils.getDetectionFlags()


        # get the data for the secondary
        XX = solver.getCatalogue(ra, dec, radius, secondaryFilter, idName, anid)
        secondaryRefsources = XX.refsources
        setSourceSetBlobsNone(secondaryRefsources)

        polyData = secondaryLookup[primaryFilter][1]   # Polynomial correction
        polyData.reverse()                             # Numpy wants decreasing powers
        polynomial = numpy.poly1d(polyData)


        # We already have the 'primary' magnitudes in the matches
        secondariesDict = dict()
        for s in secondaryRefsources:
            secondariesDict[s.getId()] = (s.getPsfFlux(), s.getPsfFluxErr())
        del secondaryRefsources


        keepref = []
        keepi = []
        for i in xrange(len(refsources)):
            ra, dec = refsources[i].getRa(), refsources[i].getDec() # ra,dec in Rads
            x, y = wcs.skyToPixel(afwCoord.Coord(afwGeom.PointD(numpy.degrees(ra), numpy.degrees(dec))))

            if x < 0 or y < 0 or x > W or y > H:
                continue

            refsources[i].setXAstrom(x)
            refsources[i].setYAstrom(y)
            if stargal[i]:
                refsources[i].setFlagForDetection(refsources[i].getFlagForDetection() | fdict["STAR"])

            # color term
            primaryMag = -2.5*numpy.log10(refsources[i].getPsfFlux())
            secondaryMag = -2.5*numpy.log10(secondariesDict[refsources[i].getId()][0])
            diff = polynomial(primaryMag - secondaryMag)
            refsources[i].setPsfFlux(10.0**(-0.4*(primaryMag+diff)))

            keepref.append(refsources[i])
            keepi.append(i)


        refsources = keepref

        if referrs is not None:
            referrs = [referrs[i] for i in keepi]
        if stargal is not None:
            stargal = [stargal[i] for i in keepi]

        measAstrom.joinMatchList(matches, refsources, first=True, log=log)
        measAstrom.joinMatchList(matches, sources, first=False, log=log)

        return matches, calib, refsources



    def calibFluxError(self, f, df, f0, df0):
        if numpy.isinf(f):
            # issues if f is inf
            return numpy.NaN
        else:
            # or df/f, even if they are finite
            try:
                df/f
            except:
                return numpy.NaN

        if f > 0.0 and f0 > 0.0:
            return (df/f + df0/f0)*f/f0
        else:
            return numpy.NaN

    def atEdge(self, bbox, x, y):

        borderWidth = 18
        x0, y0, x1, y1 = bbox
        imgWidth  = x1 - x0
        imgHeight = y1 - y0

        if x < borderWidth or imgWidth - x < borderWidth:
            return True
        if y < borderWidth or imgHeight - y < borderWidth:
            return True

        return False



