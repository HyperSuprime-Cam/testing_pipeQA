import os
import glob, re, shelve
import numpy

from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.patches import Ellipse
import matplotlib.cm as cm


def unshelveGlob(filename, testSet=None, flag='c'):

    if testSet:
        filename = os.path.join(testSet.wwwDir, filename)
    
    data = {}
    n = 0

    fsub = re.sub("-all", "*", filename)
    for f in glob.glob(fsub):
        if re.search("-all", f):
            continue
        shelf = shelve.open(f+".shelve", flag=flag)
        for k,v in shelf.items():
            if not data.has_key(k):
                data[k] = v
            else:
                if isinstance(v, numpy.ndarray):
                    data[k] = numpy.append(data[k], v)
        shelf.close()
        n += 1
    isSummary = fsub != filename
    return data, isSummary
    


def binDistrib(x, y, dy, binSizeX = 0.5, minPts = 2):
    bx  = []
    bs  = []
    by  = []
    bdy = []

    bins  = num.arange(min(x) + binSizeX, max(x), binSizeX)
    for i in range(len(bins) - 1):
        idx = num.where((x >= bins[i]) & (x < bins[i+1]))
        if len(idx[0]) < minPts:
            continue

        # find median x-values in this bin
        xmed = num.median(x[idx])
        bx.append(xmed)

        # find median y-values in this bin
        ymed = num.median(y[idx])
        by.append(ymed)

        # find width about mean y using sigIQR
        bs.append(sigIQR(y[idx]))

        # find median dy-values in this bin
        dymed = num.median(dy[idx])
        bdy.append(dymed)

    return num.array(bx), num.array(by), num.array(bs), num.array(bdy)


def plotSparseContour(sp, x, y, binSizeX, binSizeY, minCont = 500, nCont = 7):
    idx   = num.isfinite(x)
    x     = x[idx]
    y     = y[idx]

    xmin  = x.min()
    xmax  = x.max()
    ymin  = y.min()
    ymax  = y.max()
    xbins = num.arange(xmin, xmax, binSizeX)
    ybins = num.arange(ymin, ymax, binSizeY)

    cdata = num.zeros((len(ybins), len(xbins)))
    for i in range(len(x)):
        xidx = (x[i] - xmin) // binSizeX
        yidx = (y[i] - ymin) // binSizeY
        cdata[yidx][xidx] += 1


    cs    = num.arange(minCont, cdata.max(), (cdata.max() - minCont) // nCont).astype(num.int)
    c     = sp.contour(cdata, cs, origin='lower', linewidths=1, extent=(xmin,xmax,ymin,ymax))
    outer = c.collections[0]._paths

    xp = []
    yp = []
    for i in range(len(x)):
        found = [o.contains_point((x[i], y[i])) for o in outer]
        if not (True in found):
            xp.append(x[i])
            yp.append(y[i])
    sp.plot(xp, yp, 'r.', ms = 1)








def qaSetp(ticklabels, **kwargs):
    for ticklabel in ticklabels:
        for k,v in kwargs.items():
            methodName = "set_" + str(k)
            if hasattr(ticklabel, methodName):
                method = getattr(ticklabel, methodName)
                method(v)







#get precentile levels for contour plot
#uses a 'fill the bath tub method'
##################################
def getLevels(hist_data,percentile_list=[0.5]):
    """
    returns levels from histogramed data for the pecentiles listed
    the method simply 'fills up' the contours until only the desired
    percentage of the data is above the water mark. therefore ~30% of 
    the data will be contained in the first contour. This is assuming
    well-behaved data, but multiple peaks are ok
    hist_data--regularly grided histogram data (1d or 2d)
    percentile list--list of desired percentiles, needs to be in
    """

    #sort the list in ascending order
    percentile_list = sorted(percentile_list)
    
    #levels
    levels = numpy.zeros(len(percentile_list))

    #Normed histdata
    norm = numpy.sum(numpy.asarray(hist_data).flatten())
    nHist = sorted(numpy.asarray(hist_data/norm).flatten())
    cSum = numpy.cumsum(nHist)

    indices = []
    for p in percentile_list:
        x = numpy.where(cSum > p)[0]
        indices.append(x[0])

    levels = nHist[numpy.array(indices)]
    if isinstance(levels, numpy.float64):
        levels = numpy.array([levels])

    #renormalize to the hist data level
    levels = norm*levels
    return levels


#make a contour plot
def make_densityContour(axes,x,y,xlims=None,ylims=None,bins=(50,50),
                        log=False, color='g',levels=3,normed=True,
                        percentiles=False,lw=1.0,ls='solid'):
    if xlims==None:
        xlims=(x.min(),x.max())
    if ylims==None:
        ylims=(y.min(),y.max())
    x_p = x[(x>=xlims[0]) & (x<=xlims[1]) & (y>=ylims[0]) & (y<=ylims[1])]
    y_p = y[(x>=xlims[0]) & (x<=xlims[1]) & (y>=ylims[0]) & (y<=ylims[1])]

    if len(bins)==2:
        bins = (numpy.linspace(xlims[0],xlims[1],num=bins[0]),
                numpy.linspace(ylims[0],ylims[1],num=bins[1]))

    hist_xy,xedges,yedges = numpy.histogram2d(x_p,y_p,bins=bins,normed=normed)

    #if percentiles is set, get the levels from the percentiles
    if percentiles:
        if isinstance(levels,list):
            levels = getLevels(hist_xy,levels)
        else:
            levels = getLevels(hist_xy)

    if log:
        new = axes.contour( numpy.log10(numpy.flipud(numpy.rot90(hist_xy)) + 0.1), numpy.log10(levels + 0.1),
                            extent=[xedges[0], xedges[-1], 
                                    yedges[0], yedges[-1]],
                            colors=color,linestyles=ls,linewidths=lw)
        
    else:
        new = axes.contour( xedges[:-1], yedges[:-1], numpy.rot90(hist_xy), levels,
                            extent=[xedges[0], 
                                    xedges[-1], yedges[0], yedges[-1]],
                            colors=color,linestyles=ls,linewidths=lw)
        
    return new



#make a density plot and return the patches to make a colorbar
def make_densityplot(axes,x,y,xlims=None,
                     ylims=None,bins=(50,50),log=False):
    if xlims==None:
        xlims=(x.min(),x.max())
    if ylims==None:
        ylims=(y.min(),y.max())
    x_p = x[(x>=xlims[0])&(x<=xlims[1])&
            (y>=ylims[0])&(y<=ylims[1])]
    y_p = y[(x>=xlims[0])&(x<=xlims[1])&
            (y>=ylims[0])&(y<=ylims[1])]

    if len(bins)==2:
        bins = (numpy.linspace(xlims[0],xlims[1],num=bins[0]),
                numpy.linspace(ylims[0],ylims[1],num=bins[1]))

    hist_xy,xedges,yedges = numpy.histogram2d(x_p,y_p,bins=bins)

    if log:
        return axes.imshow( numpy.log10(numpy.rot90(hist_xy)+1.),
                            extent=[xedges[0], xedges[-1], 
                                    yedges[0], yedges[-1]],
                            aspect='auto',
                            interpolation='nearest',cmap=cm.gray_r)
    else:
        return axes.imshow( numpy.rot90(hist_xy),
                            extent=[xedges[0], xedges[-1], 
                                    yedges[0], yedges[-1]],
                            aspect='auto',interpolation='nearest',
                            cmap=cm.gray_r)
