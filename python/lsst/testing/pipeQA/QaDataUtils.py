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



