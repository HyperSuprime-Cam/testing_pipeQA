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
import cPickle as pickle
import eups
import datetime

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

import lsst.testing.pipeQA.TestCode as testCode
import lsst.testing.pipeQA.figures as qaFig


class QaAnalysisConfig(pexConfig.Config):
    pass


class QaAnalysisTask(pipeBase.Task):
    """Baseclass for analysis classes."""
    ConfigClass  = QaAnalysisConfig
    _DefaultName = "qaAnalysis"


    def __init__(self, testLabel=None, useCache=False, wwwCache=True, summaryProcessing="delay",
                 lazyPlot='sensor', *args, **kwargs):
        """
        @param testLabel   A name for this kind of analysis test.
        """
        pipeBase.Task.__init__(self, *args, **kwargs)

        self.testSets  = {}
        self.testLabel = testLabel

        # if we're not going to use the cached values
        # we'll have to clean the output directory on our first call
        self.useCache     = useCache
        self.clean        = not useCache
        self.wwwCache     = wwwCache
        self.summaryProcessing = summaryProcessing

        self.summOpt = {
            'delay'     : 'delay',  # make summary figures after the final CCD is processed
            'summOnly'  : 'summOnly',   # only make summary figures
            'none'      : 'none'    # don't make summary figures
            }
        

        options = ['none', 'sensor', 'all']
        if not lazyPlot in options:
            raise ValueError, "lazyPlot must be: "+ ",".join(options) + " You said: "+lazyPlot
        
        self.lazyPlot  = lazyPlot


    def combineOutputs(self, data, dataId, label=None):
        ts = self.getTestSet(data, dataId, label, noSuffix=True)
        ts.accrete()
        ts.updateCounts()
        #pass

        
    def getTestSet(self, data, dataId, label=None, noSuffix=False):
        """Get a TestSet object in the correct group.

        @param data    a QaData object
        @param dataId  a dataId dictionary
        @param label   a label for particular TestSet
        """

        dataIdStd = data.cameraInfo.dataIdCameraToStandard(dataId)
        group = dataIdStd['visit']

        raftName, ccdName = data.cameraInfo.getRaftAndSensorNames(dataId)
        ccdName = data.cameraInfo.getDetectorName(raftName, ccdName)
        
        filter = data.getFilterBySensor(dataId)
        # all sensors have the same filter, so just grab one
        key = filter.keys()[0]
        filterName = '?'
        if filter[key] is not None:
            filterName = filter[key].getName()


        ###########################################################
        # collect specific pieces of metadata for the summary panel
        summaryInfo = data.getSummaryDataBySensor(dataId)[key]


        # DATE_OBS
        dateObs = summaryInfo.get('DATE_OBS', None)
        if dateObs is None:
            dateObs = 'unknown'
        else:
            dateObs = dateObs.strftime("%Y-%m-%d")

        dateObs += ' (MJD: %.6f)' % (summaryInfo.get("MJD", 0.0))

        # RADEC
        raDec   = summaryInfo.get('RA', 'unk') + " " + summaryInfo.get('DEC', 'unk')

        # ALTAZ
        alt = summaryInfo.get('ALT', None)
        az  = str(summaryInfo.get('AZ', None))
        if alt:
            altAz   = str(alt) + " " + az + " (airmass: %.2f)" % (summaryInfo.get("AIRMASS", 0.0))
        else:
            altAz = None

        # HST
        hst = summaryInfo.get("HST", None)
        #if hst:
        #    hst = hst.strftime("%H:%M:%S")
        #else:
        #    hst = None

        # Misc. values (we won't list these unless they have a value)
        getOrIgnoreList = 'OBJECT', "EXPTIME", 'INSROT', 'FOCUSZ', 'ADCPOS', 'PA'


        if label is not None:
            label = self.__class__.__name__ + "."+ label
        else:
            label = self.__class__.__name__

        tsIdLabel = "visit-filter"
        groupId = str(group) + '-' + filterName
        tsId = str(group) + '-' + ccdName + '-' + filterName
        if data.cameraInfo.name == 'sdss':
            tsId = group

        if noSuffix:
            return testCode.TestSet(label, group=groupId, clean=self.clean, wwwCache=self.wwwCache, sqliteSuffix="")
        else:
            sqliteSuffix = ccdName

            
        if not self.testSets.has_key(tsId):
            self.testSets[tsId] = testCode.TestSet(label, group=groupId, clean=self.clean,
                                                   wwwCache=self.wwwCache, sqliteSuffix=sqliteSuffix)
            
            self.testSets[tsId].addMetadata('dataset', data.getDataName())
            self.testSets[tsId].addMetadata(tsIdLabel, groupId)

            # we'll always show these, even if showing 'None'
            self.testSets[tsId].addMetadata('DATE_OBS', dateObs)
            self.testSets[tsId].addMetadata('RaDec', raDec)
            self.testSets[tsId].addMetadata('AltAz', altAz)

            # only show these if we have them
            for k in getOrIgnoreList:
                v = summaryInfo.get(k, None)
                if v:
                    self.testSets[tsId].addMetadata(k, v)
                


            # version info
            eupsPkgs = sorted(["testing_pipeQA", "testing_displayQA", "astrometry_net_data"])
            eupsVersions = []
            for pkg in eupsPkgs:
                eupsVersions.append(eups.getSetupVersion(pkg))
            eupsVersStr = "\n".join([" ".join(pair) for pair in zip(eupsPkgs, eupsVersions)])
            self.testSets[tsId].addMetadata("EUPS_qa-"+ccdName, eupsVersStr)
            
            if hasattr(data, 'coaddTable') and not data.coaddTable is None:
                self.testSets[tsId].addMetadata('coaddTable', data.coaddTable)
            if hasattr(data, 'useForced'):
                self.testSets[tsId].addMetadata('forced', "True" if data.useForced else "False")

            key = data._dataIdToString(dataId, defineFully=True)
            sqlCache = data.sqlCache['match'].get(key, "")
            if sqlCache:
                self.testSets[tsId].addMetadata("SQL_match-"+ccdName, sqlCache)
            sqlCache = data.sqlCache['src'].get(key, "")
            if sqlCache:
                self.testSets[tsId].addMetadata("SQL_src-"+ccdName ,  sqlCache)
                
        return self.testSets[tsId]


    def __str__(self):
        testLabel = ""
        if self.testLabel is not None:
            testLabel = "."+self.testLabel
        return self.__class__.__name__ + testLabel
    

    ##########################################
    # pure virtual methods
    #########################################
    def free(self):
        """Method to free attributes to minimize memory consumption."""
        pass
    
    def test(self):
        """Method to perform tests. """
        return []

    def plot(self):
        """Method to make figures."""
        return []

