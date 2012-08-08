import os, sys, re
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import numpy
import numpy.ma as numpyMa

import matplotlib.figure as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
from matplotlib import font_manager as fm
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.patches import Ellipse
from matplotlib import cm
from matplotlib import colors

import colorsys

import QaFigureUtils as qaFigUtils
from QaFigure import QaFigure

class FpaQaFigure(QaFigure):

    def __init__(self, cameraInfo, data=None, map=None):
        """
        @param cameraInfo  CameraInfo object for the device whose focal plane area we're representing.
        @param data        Data to display
        @param map         Map areas to use.
        """

        QaFigure.__init__(self, size=(3.6,3.3))
        self.cameraInfo = cameraInfo
        self.camera = self.cameraInfo.camera
        self.centers, self.rectangles, self.raftBoundaries, self.ccdBoundaries = \
                      qaFigUtils.cameraToRectangles(self.camera)

        self.idByName = {}
        self.raftCcdByAreaLabel = {}

        self.data = {}
        self.reset()
        self.reset(data=self.map)
        self.min = None
        self.max = None
            
        # Fill the data/map values if they were provided.
        if not data is None:
            if not self.validate():
                raise Exception("Data did not pass validation.")
            for raft, ccdDict in data.items():
                for ccd, value in ccdDict.items():
                    self.data[raft][ccd] = data[raft][ccd]
                    
        if not map is None:
            for raft, ccdDict in data.items():
                for ccd, value in ccdDict.items():
                    self.map[raft][ccd] = map[raft][ccd]

    def getArray(self):
        array = []
        if not self.data is None:
            if not self.validate():
                raise Exception("Data did not pass validation.")
            for raft, ccdDict in self.data.items():
                for ccd, value in ccdDict.items():
                    if not value is None and numpy.isfinite(value):
                        array.append(value)
        return numpy.array(array)
        
                    
    def reset(self, data=None):
        """Set all values in data dictionary to None."""
        
        if data is None:
            data = self.data
            
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            data[rlabel] = {}
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                self.idByName[clabel] = ccd.getId()
                data[rlabel][clabel] = None
                areaLabel = self.getAreaLabel(rlabel, clabel)
                self.raftCcdByAreaLabel[areaLabel] = [rlabel, clabel]
                
                
    def validate(self):
        # Since we establish the structure of data in __init__, it
        # should always be valid.  Unless someone mucks with self.data
        # directly...
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            if not self.data.has_key(rlabel):
                return False
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                if not self.data[rlabel].has_key(clabel):
                    return False
        return True

            

    ###################################################################
    # map related methods
    ###################################################################

    def getAreaLabel(self, raft, ccd):
        """Get the area label to use for this raft,ccd."""
        id = self.idByName[ccd]
        name = re.sub("\s+", "_", id.getName())
        serial = "%04d" % (id.getSerial())
        return name + "--"+serial


    def setMapInfo(self):
        """Establish map areas for any defined CCDs"""
        
        for r in self.camera:
            raft = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            for c in raft:
                ccd = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                info = self.map[rlabel][clabel]
                
                if ((not info is None) and
                    self.ccdBoundaries.has_key(clabel) and
                    (not self.ccdBoundaries[clabel] is None)):
                    bound = self.ccdBoundaries[clabel]
                    x0, x1 = bound[0]
                    y0, y1 = bound[1]
                    label = self.getAreaLabel(rlabel, clabel)
                    self.addMapArea(label, [x0, y0, x1, y1], info)

                    
        
    ###################################################################
    # plot related methods
    ###################################################################


    def plotRaftBoundaries(self, sp, boundaryColors):
        for b in self.raftBoundaries:
            sp.plot(b[0], b[1], '%s-' % (boundaryColors), lw=3)
    def plotCcdBoundaries(self, sp):
        for b in self.ccdBoundaries.values():
            x0, x1 = b[0]
            y0, y1 = b[1]
            x = [x0, x0, x1, x1, x0]
            y = [y0, y1, y1, y0, y0]
            sp.plot(numpy.array(x), numpy.array(y), 'k-', lw=1.0)
    def markMissingCcds(self, sp, missingCcds):
        for b in missingCcds.values():
            x0, x1 = b[0]
            y0, y1 = b[1]
            x = [x0, x1, x0, x1]
            y = [y0, y1, y1, y0]
            sp.plot(numpy.array(x), numpy.array(y), 'k-', lw=0.5)
            
    def markFailedCcds(self, sp, failedCcds, cmap, vlimits):

        for label, value in failedCcds.items():
            b = self.ccdBoundaries[label]
            x0, x1 = b[0]
            y0, y1 = b[1]
            x = x0 + 0.2*(x1-x0)
            y = y0 + 0.8*(y1-y0)
            text = "F"

            bgColor = cmap((value-vlimits[0])/(vlimits[1]-vlimits[0]))
            r, g, b, alph = bgColor
            lum = 0.25*r + 0.7*g + 0.05*b
            clr     = "w"        if  lum < 0.5  else "k"
            sp.text(x, y, text, color=clr, horizontalalignment="left", verticalalignment="top",
                    fontsize=8)


    def labelSensors(self, sp):
        for r in self.rectangles.values():
            label = r.get_label()
            bbox  = r.get_bbox()
            xplot = 0.5 * (bbox.x0 + bbox.x1)
            yplot = bbox.y0 + 0.8*(bbox.y1 - bbox.y0)
            sp.text(xplot, yplot, label, horizontalalignment='center', fontsize = 8)
        
    def adjustTickLabels(self, sp, cb=None):
        if not cb is None:
            for tic in cb.ax.get_yticklabels():
                tic.set_size("xx-small")
        for tic in sp.get_xticklabels():
            tic.set_size("x-small")
            tic.set_rotation(22)
        for tic in sp.get_yticklabels():
            tic.set_size("x-small")
            tic.set_rotation(45)
        
    def getFpaLimits(self):
        x, y = [], []
        for r in self.rectangles.values():
            bbox  = r.get_bbox()
            x += [bbox.x0, bbox.x1]
            y += [bbox.y0, bbox.y1]
        x, y = numpy.array(x), numpy.array(y)
        return x.min(), y.min(), x.max(), y.max()
            

    def getDataArray(self):

        arr = []
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                arr.append(self.data[rlabel][clabel])
        return numpy.array(arr)
    

    def makeFigure(self, 
                   borderPix = 100,
                   boundaryColors = 'r', doLabel = False, showUndefined=False,
                   vlimits=None, cmap="jet", title=None,
                   cmapOver=None, cmapUnder=None,
                   failLimits=None, failColor='k',
                   ):
        """Make the figure.

        @param borderPix        width of border in pixels
        @param boundaryColors   matplotlib color specifier for border
        @param doLabel          add ccd labels to the figure
        @param showUndefined    show sensors even if their value is None
        @param vlimits          [low, high] limits for colormap
        @param cmap             string name of an matplotlib.cm colormap
        @param title            title of the figure
        @param cmapOver         Color to use if above cmap high limit.
        @param cmapUnder        Color to use if below cmap low limit.
        @param failLimits       Limits to mark failure.
        @param failColor        Color to use to mark failed sensors.
        """

        
        if vlimits is None:
            arr = self.getDataArray()
            vlimits = [arr.min(), arr.max()]

        if failLimits is None:
            failLimits = vlimits

        histWidth = 0.05
        useHist = True
        if useHist:
            left, right, bottom, top = 0.195, 0.95-2.0*histWidth, 0.15, 0.9
        else:
            left, right, bottom, top = 0.195, 0.95, 0.15, 0.9
        self.fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)
        
        sp     = self.fig.gca()

        values = []  # needs to be synchronized with self.rectangles
        patches = []
        allValues = []
        missingCcds = {}
        failedCcds = {}
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                value = self.data[rlabel][clabel]
                allValues.append(value)
                #if (not value is None) or (showUndefined):
                if value is None:
                    value = numpy.NaN
                    missingCcds[clabel] = self.ccdBoundaries[clabel]
                if value < failLimits[0] or value > failLimits[1]:
                    failedCcds[clabel] = value
                #if value > failLimits[1]:
                #    failedCcds[clabel] = 1
                values.append(value)
                patches.append(self.rectangles[clabel])


        if not vlimits is None:
            norm = colors.Normalize(vmin=vlimits[0], vmax=vlimits[1], clip=False)
        else:
            norm = colors.Normalize()
            
        if len(patches) == 0:
            patches = self.rectangles.values()
            values = allValues

        cmap = getattr(cm, cmap)
        cmap.set_bad('k', 0.2)
        if not cmapOver is None:
            cmap.set_over(cmapOver, 1.0)
        if not cmapUnder is None:
            cmap.set_under(cmapUnder, 1.0)
        p = PatchCollection(patches, norm=norm, cmap=cmap)
        value_array = numpy.array(values)
        masked_value_array = numpyMa.masked_where(numpy.isnan(value_array), value_array)
        p.set_array(masked_value_array)
        cb = self.fig.colorbar(p)
        sp.add_collection(p)

        ##############################
        # put a histogram on the side
        if useHist:
            nbins = 30
            histValues = numpy.array(values)
            finiteValues = histValues[numpy.where(numpy.isfinite(histValues))[0]]
            left = right+0.05
            axH = self.fig.add_axes([left, bottom, histWidth, top-bottom])
            binwid = float(vlimits[1] - vlimits[0])/nbins
            if binwid > 0.0 and len(finiteValues) > 0:
                eps = 1.0e-4*binwid
                #print finiteValues, nbins, vlimits
                nu, bu, pu = axH.hist(finiteValues, bins=nbins,
                                      range=[vlimits[0]-eps,vlimits[1]+eps],
                                      orientation='horizontal', color='#aaaaaa', fill=True)
                #nu = numpy.array([1.0])
                axH.set_ylim(vlimits)
                xmax = 1.2*nu.max()
                if xmax <= 0.0:
                    xmax = 1.0
                axH.set_xlim([0, xmax])

            axH.set_xlabel('')
            axH.set_ylabel('')
            axH.set_xticklabels([])
            axH.set_yticklabels([])
            self.fig.sca(sp)

        self.plotRaftBoundaries(sp, boundaryColors)
        self.plotCcdBoundaries(sp)
        self.markMissingCcds(sp, missingCcds)
        self.markFailedCcds(sp, failedCcds, cmap, vlimits)
        if self.cameraInfo.doLabel or doLabel:
            self.labelSensors(sp)

        if not title is None:
            self.fig.text(0.5, 0.94, title, horizontalalignment="center", fontsize=10)
            #sp.set_title(title, fontsize=12)
        sp.set_xlabel("Focal Plane X", fontsize = 9)
        sp.set_ylabel("Focal Plane Y", fontsize = 9)

        self.adjustTickLabels(sp, cb)

        x0, y0, x1, y1 = self.getFpaLimits()
        sp.set_xlim((x0 - borderPix, x1 + borderPix))
        sp.set_ylim((y0 - borderPix, y1 + borderPix))

        self.setMapInfo()
        


class VectorFpaQaFigure(FpaQaFigure):

    def __init__(self, cameraInfo, data=None, map=None):
        FpaQaFigure.__init__(self, cameraInfo, data=data, map=map)


    def getDataArray(self):

        arr = []
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                arr.append(self.data[rlabel][clabel][2])
        return numpy.array(arr)
    
    def makeFigure(self, 
                   borderPix = 100,
                   boundaryColors = 'r', doLabel = False, showUndefined=False,
                   vlimits=None, cmap="jet", title=None,
                   cmapOver=None, cmapUnder=None,
                   failLimits=None, failColor='k',
                   ):
        """Make the figure.

        @param borderPix        width of border in pixels
        @param boundaryColors   matplotlib color specifier for border
        @param doLabel          add ccd labels to the figure
        @param showUndefined    show sensors even if their value is None
        @param vlimits          [low, high] limits for colormap
        @param cmap             string name of an matplotlib.cm colormap
        @param title            title of the figure
        @param cmapOver         Color to use if above cmap high limit.
        @param cmapUnder        Color to use if below cmap low limit.
        @param failLimits       Limits to mark failure.
        @param failColor        Color to use to mark failed sensors.
        """

        if vlimits is None:
            arr = self.getDataArray()
            vlimits = [arr.min(), arr.max()]
        if failLimits is None:
            failLimits = vlimits

        histWidth = 0.05
        useHist = True
        if useHist:
            left, right, bottom, top = 0.195, 0.95-2.0*histWidth, 0.15, 0.9
        else:
            left, right, bottom, top = 0.195, 0.95, 0.15, 0.9
            
        self.fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top)

        colorValues = []  # needs to be synchronized with self.rectangles
        patches = []
        allValues = []

        radiansWrtX = {}
        lenInPix = {}
        colorScalar = {}
        haveColors = False
        
        missingCcds = {}
        failedCcds = {}
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                values = self.data[rlabel][clabel]
                defaultLen = 1500.0
                if isinstance(values, list):
                    if len(values) == 3:
                        radiansWrtXtmp, lenInPixtmp, colorScalartmp = values
                    elif len(values) == 2:
                        radiansWrtXtmp, lenInPixtmp = values
                        colorScalartmp = 0.0
                    else:
                        raise Exception("values for Vector must be float or [radians, lenInPix, [colorFloat]].")
                    if lenInPixtmp is None:
                        lenInPixtmp = defaultLen
                        
                else:
                    if not values is None:
                        values = float(values)
                        defaultLen = 1500.0  # this works for most ccds, but should change to eg. ccdwidth/2
                    else:
                        values = 0.0
                        defaultLen = 0.0
                    radiansWrtXtmp, lenInPixtmp, colorScalartmp = values, defaultLen, None

                if lenInPixtmp < 0:
                    lenInPixtmp = 0
                #print clabel, radiansWrtXtmp, lenInPixtmp, colorScalartmp
                
                #print clabel, radiansWrtXtmp, lenInPixtmp, colorScalartmp
                lenInPix[clabel]    = lenInPixtmp
                allValues.append(radiansWrtXtmp)
                if (not radiansWrtXtmp is None): # or (showUndefined):
                    if not colorScalartmp is None:
                        colorValues.append(colorScalartmp)
                        patches.append(self.rectangles[clabel])
                        colorScalar[clabel] = colorScalartmp
                        haveColors = True
                        if colorScalartmp < failLimits[0] or colorScalartmp > failLimits[1]:
                            failedCcds[clabel] = colorScalartmp
                        #if colorScalartmp > failLimits[1]:
                        #    failedCcds[clabel] = 1
                    else:
                        colorValues.append(numpy.NaN)
                        patches.append(self.rectangles[clabel])
                        missingCcds[clabel] = self.ccdBoundaries[clabel]
                    radiansWrtX[clabel] = radiansWrtXtmp

                else:
                    colorValues.append(numpy.NaN)
                    patches.append(self.rectangles[clabel])
                    missingCcds[clabel] = self.ccdBoundaries[clabel]
                    

        if not vlimits is None:
            norm = colors.Normalize(vmin=vlimits[0], vmax=vlimits[1], clip=False)
        else:
            norm = colors.Normalize()


        sp     = self.fig.gca() #add_axes([left, bottom, right-left, top-bottom]) #gca()

        ##############################
        # put a histogram on the side
        if useHist:
            nbins = 30
            histValues = numpy.array(colorValues)
            finiteValues = histValues[numpy.where(numpy.isfinite(histValues))[0]]
            left = right+0.05
            axH = self.fig.add_axes([left, bottom, histWidth, top-bottom])
            binwid = float(vlimits[1] - vlimits[0])/nbins
            if binwid > 0.0 and len(finiteValues) > 0:
                eps = 1.0e-4*binwid
                nu, bu, pu = axH.hist(finiteValues, bins=nbins,
                                      range=[vlimits[0]-eps,vlimits[1]+eps],
                                      orientation='horizontal', color='#aaaaaa', fill=True)
                axH.set_ylim(vlimits)
                axH.set_xlim([0, 1.2*nu.max()])
            axH.set_xlabel('')
            axH.set_ylabel('')
            axH.set_xticklabels([])
            axH.set_yticklabels([])
            self.fig.sca(sp)
        
        if len(patches) > 0:

            cmap = getattr(cm, cmap)
            cmap.set_bad('k', 0.2)
            if not cmapOver is None:
                cmap.set_over(cmapOver, 1.0)
            if not cmapUnder is None:
                cmap.set_under(cmapUnder, 1.0)

            p = PatchCollection(patches, norm=norm, cmap=cmap)
            value_array = numpy.array(colorValues)
            masked_value_array = numpyMa.masked_where(numpy.isnan(value_array), value_array)
            p.set_array(masked_value_array)
            if haveColors or (not vlimits is None):
                cb = self.fig.colorbar(p)
            sp.add_collection(p)


        for label, angle in radiansWrtX.items():
            xc, yc = self.centers[label]
            arrowLen = lenInPix[label]
            x = xc - 0.5*arrowLen*numpy.cos(angle)
            y = yc - 0.5*arrowLen*numpy.sin(angle)
            dx = arrowLen*numpy.cos(angle)
            dy = arrowLen*numpy.sin(angle)
            if numpy.abs(dx) < 1.0e-15 or numpy.abs(dy) < 1.0e-15:
                continue
            sp.arrow(x, y, dx, dy) #, ec="k", lw=3)

        self.plotRaftBoundaries(sp, boundaryColors)
        self.plotCcdBoundaries(sp)
        self.markMissingCcds(sp, missingCcds)
        self.markFailedCcds(sp, failedCcds, cmap, vlimits)
        if self.cameraInfo.doLabel or doLabel:
            self.labelSensors(sp)

        if not title is None:
            self.fig.text(0.5, 0.94, title, horizontalalignment="center", fontsize=10)
            #sp.set_title(title, fontsize=12)
        sp.set_xlabel("Focal Plane X", fontsize = 9)
        sp.set_ylabel("Focal Plane Y", fontsize = 9)

        if (not haveColors) and (vlimits is None):
            cb = None
        self.adjustTickLabels(sp, cb)

        x0, y0, x1, y1 = self.getFpaLimits()
        sp.set_xlim((x0 - borderPix, x1 + borderPix))
        sp.set_ylim((y0 - borderPix, y1 + borderPix))



        
        self.setMapInfo()

