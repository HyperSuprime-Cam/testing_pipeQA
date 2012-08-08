import os, re
import numpy

fields = [
    "refObjectId",
    "isStar",
    "ra", "decl",
    "uMag", "gMag", "rMag", "iMag", "zMag", "yMag"
    ]



class SimRefObjectSet(list):
    def push_back(self, sro):
        self.append(sro)


class SimRefObject(object):

    flookup = {
        "u":0, "g": 1, "r":2, "i":3, "z":4, "y":5,
        "B":1, 'V':2, 'R':2, 'I':3,
        }

    def __init__(self, *sroStuff):

        if len(sroStuff) == 10:
            self.refObjectId = sroStuff[0]
            self.isStar = sroStuff[1]
            self.radec = numpy.array(sroStuff[2:4])
            self.mag = numpy.array(sroStuff[4:10], dtype=numpy.float32)

        elif len(sroStuff) == 0:
            self.refObjectId = 0
            self.isStar = 0
            self.radec = numpy.zeros(2, dtype=numpy.float64)
            self.mag = numpy.zeros(6, dtype=numpy.float32)

        else:
            raise NotImplementedError, "Cannot instantiate SimRefObject with" + \
                str(len(sroStuff)) + " args."

        #self.refObjectId = 0
        #self.isStar      = 0
        #"varClass",       
        #self.radec       = numpy.zeros(0)
        #"gLat",           
        #"gLon",           
        #"sedName",        
        #self.mag         = numpy.zeros(5)
        #"muRa",           
        #"muDecl",         
        #"parallax",       
        #"vRad",           
        #"redshift",       
        #"semiMajorBulge", 
        #"semiMinorBulge", 
        #"semiMajorDisk",  
        #"semiMinorDisk",  
        #"uCov",           
        #"gCov",           
        #"rCov",           
        #"iCov",           
        #"zCov",           
        #"yCov",           



    def getId(self): return self.refObjectId
    def setId(self, val): self.refObjectId = val
    def getIsStar(self): return self.isStar
    def setIsStar(self, val): self.isStar = val

    def getRa(self):       return self.radec[0]
    def setRa(self, ra):   self.radec[0] = ra
    def getDecl(self):      return self.radec[1]
    def setDecl(self, dec):      self.radec[1] = dec

    def setMag(self, magNew, filter):
        i = SimRefObject.flookup[filter]
        self.mag[i] = newMag

    def setFlux(self, fluxNew, filter):
        i = SimRefObject.flookup[filter]
        if fluxNew > 0 and not numpy.isnan(fluxNew):
            self.mag[i] = -2.5*numpy.log10(fluxNew)
        else:
            self.mag[i] = numpy.NaN

    def getMag(self, filter):
        return self.mag[SimRefObject.flookup[filter]]

    def getFlux(self, filter):
        return 10.0**(-0.4*self.mag[SimRefObject.flookup[filter]])

    #def getCov(self, filter):
    #    cov = getattr(self, filter+"Cov")
    #    return cov

