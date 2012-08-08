import os, re
import numpy

sourceClass = 'python'
if os.environ.has_key('SOURCECLASS'):
    if re.search('afw', os.environ['SOURCECLASS']):
        sourceClass = 'afw'
    else:
        sourceClass = 'python'
        
        
useSource = "afw" #sourceClass
useRefSource = sourceClass

# these mimic the old measAlg.Flags values
# since we're wrapping the Source, we'll use these
STAR          = 0x1
SATUR_CENTER  = 0x2
EDGE          = 0x4
INTERP_CENTER = 0x8

#################################################################
# RefSource
import lsst.afw.table  as afwTab
import QaDataUtils as qaDataUtils

class _RefCatalog(object):
    
    def __init__(self):
        self.schema = afwTab.SourceTable.makeMinimalSchema()

        #setMethods = [x for x in qaDataUtils.getSourceSetAccessors()]
        
        setMethods = ["Ra", "Dec", "PsfFlux", "ApFlux", "ModelFlux", "InstFlux"]

        self.keyDict = {}
        #for sm in setMethods0:
        #    key = self.schema.addField(sm, type="D")
        #    self.keyDict[sm] = key
            
        self.setKeys = []
        for sm in setMethods:
            if sm == 'Id':
                continue
            key = self.schema.addField(sm, type="D")
            self.setKeys.append(key)
            self.keyDict[sm] = key
            setattr(self, sm+"Key", key)
            
        self.table = afwTab.SourceTable.make(self.schema)
        self.catalog = afwTab.SourceCatalog(self.table)


        
class _RefSource(object):

    def __init__(self):
        self.ints = numpy.zeros(2, dtype=numpy.int)
        self.radec = numpy.zeros(2)
        self.flux  = numpy.zeros(4, dtype=numpy.float32)
        #self._id               = 0            
        #self._flagForDetection = 0x0
        #self._ra               = 0.0
        #self._dec              = 0.0
        #self._psfFlux          = 0.0
        #self._apFlux           = 0.0
        #self._modelFlux        = 0.0

    def setPhotometry(self, val):
        pass
    def setAstrometry(self, val):
        pass
    def setShape(self, val):
        pass

    def getId(self):      return self.ints[0]
    def setId(self, val): self.ints[0] = val
    def getFlagForDetection(self):  return self.ints[1]
    def setFlagForDetection(self, val): self.ints[1] = val

    def getRa(self):        return self.radec[0]
    def setRa(self, val):   self.radec[0] = val
    def getDec(self):       return self.radec[1]
    def setDec(self, val):  self.radec[1] = val

    def getPsfFlux(self):      return self.flux[0]
    def setPsfFlux(self, val): self.flux[0] = val
    def getApFlux(self):       return self.flux[1]
    def setApFlux(self, val):  self.flux[1] = val
    def getInstFlux(self):     return self.flux[2]
    def setInstFlux(self, val): self.flux[2] = val
    def getModelFlux(self):    return self.flux[3]
    def setModelFlux(self, val): self.flux[3] = val




if useRefSource == 'afw':
    from lsst.afw.detection import Source as RefSource
else:
    RefSource = _RefSource
    RefCatalog = _RefCatalog


    


################################################
# a local Catalog wrapper
class _Catalog(object):

    def __init__(self):
        self.schema = afwTab.SourceTable.makeMinimalSchema()

        setMethods = [x for x in qaDataUtils.getSourceSetAccessors()]

        #self.schema.addField('Id', type="L")

        self.setKeys = []
        self.keyDict = {}
        for sm in setMethods:
            if sm == 'Id':
                continue
            key = self.schema.addField(sm, type="D")
            self.setKeys.append(key)
            self.keyDict[sm] = key
            setattr(self, sm+"Key", key)

        self.table = afwTab.SourceTable.make(self.schema)
        self.catalog = afwTab.SourceCatalog(self.table)

        
##################################################
# a local Source object
class _Source(object):

    iPsfIxx     = 0
    iPsfIxxErr  = 1
    iPsfIyy     = 2
    iPsfIyyErr  = 3
    iPsfIxy     = 4
    iPsfIxyErr  = 5
    iResolution = 6
    iE1         = 7   
    iE1Err      = 8
    iE2         = 9   
    iE2Err      = 10
    iShear1     = 11
    iShear1Err  = 12
    iShear2     = 13
    iShear2Err  = 14
    iSigma      = 15
    iSigmaErr   = 16


    def __init__(self):
        self.ints              = numpy.zeros(2, dtype=int)
        self.radec             = numpy.zeros(2)
        self.astrom            = numpy.zeros(2, dtype=numpy.float32)
        self.flux              = numpy.zeros(4, dtype=numpy.float32)
        self.fluxErr           = numpy.zeros(4, dtype=numpy.float32)
        self.izz               = numpy.zeros(3, dtype=numpy.float32)
        self.izzErr            = numpy.zeros(3, dtype=numpy.float32)
        self.shape             = numpy.zeros(17, dtype=numpy.float32)


    def setPhotometry(self, val): pass
    def setAstrometry(self, val): pass
    def setShape(self, val):      pass

    def getId(self):            return self.ints[0]
    def setId(self, val):       self.ints[0] = val
    def getFlagForDetection(self):      return self.ints[1] 
    def setFlagForDetection(self, val): self.ints[1] = val


    def getRa(self):            return self.radec[0]
    def setRa(self, val):       self.radec[0] = val
    def getDec(self):           return self.radec[1]
    def setDec(self, val):      self.radec[1] = val

    def getXAstrom(self):       return self.astrom[0]
    def setXAstrom(self, val):  self.astrom[0] = val
    def getYAstrom(self):       return self.astrom[1]
    def setYAstrom(self, val):  self.astrom[1] = val

    def getPsfFlux(self):       return self.flux[0]
    def setPsfFlux(self, val):  self.flux[0] = val
    def getPsfFluxErr(self):      return self.fluxErr[0]
    def setPsfFluxErr(self, val): self.fluxErr[0] = val

    def getApFlux(self):        return self.flux[1]
    def setApFlux(self, val):   self.flux[1] = val
    def getApFluxErr(self):       return self.fluxErr[1]
    def setApFluxErr(self, val):  self.fluxErr[1] = val

    def getInstFlux(self):      return self.flux[2]
    def setInstFlux(self, val): self.flux[2] = val
    def getInstFluxErr(self):     return self.fluxErr[2]
    def setInstFluxErr(self, val): self.fluxErr[2] = val

    def getModelFlux(self):      return self.flux[3]
    def setModelFlux(self, val): self.flux[3] = val
    def getModelFluxErr(self):   return self.fluxErr[3]
    def setModelFluxErr(self, val):  self.fluxErr[3] = val

    def getIxx(self):             return self.izz[0]
    def setIxx(self, val):        self.izz[0] = val
    def getIxxErr(self):          return self.izzErr[0]
    def setIxxErr(self, val):     self.izzErr[0] = val
    def getIyy(self):             return self.izz[1]
    def setIyy(self, val):        self.izz[1] = val
    def getIyyErr(self):          return self.izzErr[1]
    def setIyyErr(self, val):     self.izzErr[1] = val
    def getIxy(self):              return self.izz[2]
    def setIxy(self, val):         self.izz[2] = val
    def getIxyErr(self):           return self.izzErr[2]
    def setIxyErr(self, val):      self.izzErr[2] = val


    def getPsfIxx(self):          return self.shape[Source.iPsfIxx]
    def setPsfIxx(self, val):     self.shape[Source.iPsfIxx] = val
    def getPsfIxxErr(self):       return self.shape[Source.iPsfIxxErr]
    def setPsfIxxErr(self, val):  self.shape[Source.iPsfIxxErr] = val
    def getPsfIyy(self):          return self.shape[Source.iPsfIyy]
    def setPsfIyy(self, val):     self.shape[Source.iPsfIyy] = val
    def getPsfIyyErr(self):       return self.shape[Source.iPsfIyyErr]
    def setPsfIyyErr(self, val):  self.shape[Source.iPsfIyyErr] = val
    def getPsfIxy(self):          return self.shape[Source.iPsfIxy]
    def setPsfIxy(self, val):     self.shape[Source.iPsfIxy] = val
    def getPsfIxyErr(self):       return self.shape[Source.iPsfIxyErr]
    def setPsfIxyErr(self, val):  self.shape[Source.iPsfIxyErr] = val

    def getResolution(self):      return self.shape[Source.iResolution]
    def setResolution(self, val): self.shape[Source.iResolution] = val

    def getE1(self):              return self.shape[Source.iE1]   
    def setE1(self, val):         self.shape[Source.iE1] = val
    def getE1Err(self):           return self.shape[Source.iE1Err]
    def setE1Err(self, val):      self.shape[Source.iE1Err] = val
    def getE2(self):              return self.shape[Source.iE2]   
    def setE2(self, val):         self.shape[Source.iE2] = val
    def getE2Err(self):           return self.shape[Source.iE2Err]
    def setE2Err(self, val):      self.shape[Source.iE2Err] = val

    def getShear1(self):          return self.shape[Source.iShear1]
    def setShear1(self, val):     self.shape[Source.iShear1] = val
    def getShear1Err(self):       return self.shape[Source.iShear1Err]
    def setShear1Err(self, val):  self.shape[Source.iShear1Err] = val
    def getShear2(self):          return self.shape[Source.iShear2]
    def setShear2(self, val):     self.shape[Source.iShear2] = val
    def getShear2Err(self):       return self.shape[Source.iShear2Err]
    def setShear2Err(self, val):  self.shape[Source.iShear2Err] = val

    def getSigma(self):           return self.shape[Source.iSigma]
    def setSigma(self, val):      self.shape[Source.iSigma] = val
    def getSigmaErr(self):        return self.shape[Source.iSigmaErr]
    def setSigmaErr(self, val):   self.shape[Source.iSigmaErr] = val


#     def getPsfIxx(self):             return self._psfIxx
#     def setPsfIxx(self, val):             self._psfIxx = val
#     def getPsfIxxErr(self):             return self._psfIxxErr
#     def setPsfIxxErr(self, val):             self._psfIxxErr = val
#     def getPsfIyy(self):             return self._psfIyy
#     def setPsfIyy(self, val):             self._psfIyy = val
#     def getPsfIyyErr(self):             return self._psfIyyErr
#     def setPsfIyyErr(self, val):             self._psfIyyErr = val
#     def getPsfIxy(self):             return self._psfIxy
#     def setPsfIxy(self, val):             self._psfIxy = val
#     def getPsfIxyErr(self):             return self._psfIxyErr
#     def setPsfIxyErr(self, val):             self._psfIxyErr = val

#     def getResolution(self):             return self._resolution
#     def setResolution(self, val):             self._resolution = val

#     def getE1(self):                return self._e1
#     def setE1(self, val):             self._e1 = val
#     def getE1Err(self):             return self._e1Err
#     def setE1Err(self, val):             self._e1Err = val
#     def getE2(self):                return self._e2
#     def setE2(self, val):             self._e2 = val
#     def getE2Err(self):             return self._e2Err
#     def setE2Err(self, val):             self._e2Err = val

#     def getShear1(self):             return self._shear1
#     def setShear1(self, val):             self._shear1 = val
#     def getShear1Err(self):             return self._shear1Err
#     def setShear1Err(self, val):             self._shear1Err = val
#     def getShear2(self):             return self._shear2
#     def setShear2(self, val):             self._shear2 = val
#     def getShear2Err(self):             return self._shear2Err
#     def setShear2Err(self, val):             self._shear2Err = val

#     def getSigma(self):             return self._sigma
#     def setSigma(self, val):             self._sigma = val
#     def getSigmaErr(self):             return self._sigmaErr
#     def setSigmaErr(self, val):             self._sigmaErr = val


################################################################
# Source
    
if useSource == 'afw':    
    Catalog = _Catalog
    Source = _Source
    
    
else:
    Source = _Source
