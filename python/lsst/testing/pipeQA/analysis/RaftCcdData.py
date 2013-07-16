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
import lsst.afw.math as afwMath

class RaftCcdData(object):

    def __init__(self, detector, initValue=0.0):
        self.detector = detector
        self.data = {}
        self.keys = []
        self.reset(initValue)
        
        self.cache = None

    def raftCcdKeys(self):
        keyList = []
        for raft in sorted(self.data.keys()):
            for ccd in sorted(self.data[raft].keys()):
                keyList.append([raft, ccd])
        return keyList
        

    def listKeysAndValues(self):
        kvList = []
        for raft in sorted(self.data.keys()):
            for ccd in sorted(self.data[raft].keys()):
                kvList.append([raft, ccd, self.data[raft][ccd]])
        return kvList


    def reset(self, value=0.0):
        for key, detector in self.detector.items():
            raft = None
            if hasattr(detector, 'getParent'):
                raft = detector.getParent().getId().getName()
            ccd = None
            if hasattr(detector, 'getId'):
                ccd = detector.getId().getName()
            if (not raft is None) and (not ccd is None):
                if not self.data.has_key(raft):
                    self.data[raft] = {}
                if not self.data[raft].has_key(ccd):
                    self.data[raft][ccd] = value

    def set(self, raft, ccd, value):
        self.data[raft][ccd] = value
    def get(self, raft, ccd, default=None):
        if self.data.has_key(raft) and self.data[raft].has_key(ccd):
            return self.data[raft][ccd]
        else:
            return default

    def cacheValues(self, recache=False):
        if (self.cache is None) or recache:
            self.cache = numpy.array([])
            for raft, ccdDict in self.data.items():
                for ccd, value in ccdDict.items():
                    self.cache = numpy.append(self.cache, value)

    def summarize(self, methodName, recache=False, default=numpy.NaN):
        self.cacheValues(recache)
        if methodName == 'median':
            return numpy.median(self.cache)
        if len(self.cache) > 0:
            method = getattr(self.cache, methodName)
            value = method()
        else:
            value = default
        return value
    
        

class RaftCcdVector(RaftCcdData):

    def __init__(self, detector):
        RaftCcdData.__init__(self, detector, initValue=numpy.array([], dtype=numpy.float))

    def xxxlistKeysAndValues(self, methodName=None, nHighest=None, nLowest=None):
        kvList = []
        for raft in sorted(self.data.keys()):
            for ccd in sorted(self.data[raft].keys()):
                # if not method specified, return the list
                if methodName is None:
                    value = self.data[raft][ccd]

                # otherwise reduce the list to a number: eg. mean, std, median
                else:
                    finite = numpy.where( numpy.isfinite(self.data[raft][ccd]) )
                    dtmp = self.data[raft][ccd][finite]
                    if nHighest is not None:
                        dtmp.sort()
                        dtmp = dtmp[-nHighest:]
                    if (not nLowest is None) and (nHighest is None):
                        dtmp.sort()
                        dtmp = dtmp[0:nLowest]

                    value = None
                    if methodName == 'median':
                        value = numpy.median(dtmp)
                    else:
                        method = getattr(dtmp, methodName)
                        value = method()
                kvList.append([raft, ccd, value])
        return kvList


    def listKeysAndValues(self, methodName=None, nHighest=None, nLowest=None, limits=None):

        methods = {
            "median" : afwMath.MEDIAN,
            "meanclip" : afwMath.MEANCLIP,
            "stdevclip" : afwMath.STDEVCLIP,
            "mean" : afwMath.MEAN,
            "stdev" : afwMath.STDEV,
            }
        
        kvList = []
        for raft in sorted(self.data.keys()):
            for ccd in sorted(self.data[raft].keys()):
                dtmp = self.data[raft][ccd].copy()
                if nHighest is not None:
                    dtmp.sort()
                    dtmp = dtmp[-nHighest:]
                if (not nLowest is None) and (nHighest is None):
                    dtmp.sort()
                    dtmp = dtmp[0:nLowest]
                    
                if limits is not None:
                    lo, hi = limits
                    w = numpy.where( (dtmp > lo) & (dtmp < hi) )
                    dtmp = dtmp[w]
                    
                stat = afwMath.makeStatistics(dtmp, afwMath.NPOINT | methods[methodName])
                value = stat.getValue(methods[methodName])
                n = stat.getValue(afwMath.NPOINT)
                kvList.append([raft, ccd, value, n])
                
        return kvList

        
    def reset(self, initValue=numpy.array([])):
        RaftCcdData.reset(self, initValue)

    def append(self, raft, ccd, value):
        self.data[raft][ccd] = numpy.append(self.data[raft][ccd], value)

