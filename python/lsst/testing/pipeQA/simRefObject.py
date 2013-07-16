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

