import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.cameraGeom.utils as cameraGeomUtils
import numpy

import pylab
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.patches import Ellipse


import QaFigureUtils as qaFigUtils


class QaFig(object):

    def __init__(self):
        self.fig         = pylab.figure()

    def reset(self):
        self.fig.clf()
	

    def validate(self):
	pass
    
    def makeFigure(self):
        # override
        pass

    def getFigure(self):
	return self.fig

    def savefig(self, path, **kwargs):
	self.fig.savefig(path, **kwargs)



class FpaQaFigure(QaFig):

    def __init__(self, camera, data=None):
        QaFig.__init__(self)
	self.camera = camera
        self.rectangles, self.boundaries = qaFigUtils.cameraToRectangles(self.camera)

        # To be filled in by child class
	if not data is None:
	    if not self.validate():
		raise Exception("Data did not pass validation.")
	else:
	    self.data = {}
	    self.reset()
	
    def reset(self):
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            self.data[rlabel] = {}
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                self.data[rlabel][clabel] = 0.0
		
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


    def makeFigure(self, 
                   DPI = 100., size = (1024, 1024), borderPix = 100,
                   boundaryColors = 'r', doLabel = False):

        if not self.validate():
            Trace("lsst.testing.pipeQA.FpaFigure", 1, "Invalid Data")
            return None

        self.fig.set_size_inches(size[0] / DPI, size[1] / DPI)
        
        sp     = self.fig.gca()
        values = []  # needs to be synchronized with self.rectangles
        for r in self.camera:
            raft   = cameraGeom.cast_Raft(r)
            rlabel = raft.getId().getName()
            for c in raft:
                ccd    = cameraGeom.cast_Ccd(c)
                clabel = ccd.getId().getName()
                values.append(self.data[rlabel][clabel])

        p = PatchCollection(self.rectangles)
        p.set_array(numpy.array(values))
        cb = self.fig.colorbar(p)
        sp.add_collection(p)

        for b in self.boundaries:
            sp.plot(b[0], b[1], '%s-' % (boundaryColors), lw=3)

        if doLabel:
            for r in self.rectangles:
                label = r.get_label()
                bbox  = r.get_bbox()
                xplot = 0.5 * (bbox.x0 + bbox.x1)
                yplot = bbox.y1 - size[1]//2
                sp.text(xplot, yplot, label, horizontalalignment='center', fontsize = 8, weight = 'bold')

        sp.set_xlabel("Focal Plane X", fontsize = 20, weight = 'bold')
        sp.set_ylabel("Focal Plane Y", fontsize = 20, weight = 'bold')

        xmin = +1e10
        ymin = +1e10
        xmax = -1e10
        ymax = -1e10
        for r in self.rectangles:
            bbox  = r.get_bbox()
            
            if (bbox.x0 < xmin):
                xmin = bbox.x0
            if (bbox.x1 > xmax):
                xmax = bbox.x1
            if (bbox.y0 < ymin):
                ymin = bbox.y0
            if (bbox.y1 > ymax):
                ymax = bbox.y1
        sp.set_xlim((xmin - borderPix, xmax + borderPix))
        sp.set_ylim((ymin - borderPix, ymax + borderPix))
    