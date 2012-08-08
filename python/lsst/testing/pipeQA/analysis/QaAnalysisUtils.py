import sys, os, re, copy
import numpy

import lsst.afw.math        as afwMath
import lsst.meas.algorithms as measAlg
import lsst.testing.pipeQA.source as pqaSource

# NOTE: please replace this with (s.getFlagForDetection() & pqaSource.STAR)
#       when we eventually start setting it.
def isStarMoment(ss):
    """Quick and dirty isStar() based on comparison to psf ixx/yy/xy"""

    vixx, vixy, viyy = [], [], []
    for s0 in ss:

        s = s0
        if isinstance(s, list):
            sref, s, d = s0
            
        ixx, ixy, iyy    = s.getIxx(), s.getIyy(), s.getIxy()

        vixx.append(ixx)
        vixy.append(ixy)
        viyy.append(iyy)

    sxx = afwMath.makeStatistics(vixx, afwMath.MEANCLIP | afwMath.STDEVCLIP)
    sxy = afwMath.makeStatistics(vixy, afwMath.MEANCLIP | afwMath.STDEVCLIP)
    syy = afwMath.makeStatistics(viyy, afwMath.MEANCLIP | afwMath.STDEVCLIP)

    for s0 in ss:

        s = s0
        if isinstance(s, list):
            sref, s, d = s0
            
        xxOk = (sxx.getValue(afwMath.MEANCLIP) - s.getIxx())/sxx.getValue(afwMath.STDEVCLIP) < 3.0
        xyOk = (sxy.getValue(afwMath.MEANCLIP) - s.getIxy())/sxx.getValue(afwMath.STDEVCLIP) < 3.0
        yyOk = (syy.getValue(afwMath.MEANCLIP) - s.getIyy())/sxx.getValue(afwMath.STDEVCLIP) < 3.0

        if xxOk and xyOk and yyOk:
            s.setFlagForDetection(s.getFlagForDetection() | pqaSource.STAR)



def isStarDeltaMag(ss):

    for s0 in ss:

        s = s0
        if isinstance(s, list):
            sref, s, d = s0
        f_psf, f_mod, f_inst = s.getPsfFlux(), s.getModelFlux(), s.getInstFlux()

        # allow either inst or mod fluxes to discriminate
        for f in [f_mod, f_inst]:
            if f_psf > 0 and f > 0:
                m_psf, m = -2.5*numpy.log10(f_psf), -2.5*numpy.log10(f)

                if abs(m_psf - m) < 0.1:
                    s.setFlagForDetection(s.getFlagForDetection() | pqaSource.STAR)


def isStar(ss):
    return isStarDeltaMag(ss)



def lineFit(x, y, dy=None):
    """A standard linear least squares line fitter with errors and chi2."""
    
    N = len(x)
    if N < 2:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    no_err = False
    if dy is None:
        no_err = True
        dy = numpy.ones(N)
        
    var = dy**2
    S  = (1.0/var).sum()
    Sx = (x/var).sum()
    Sy = (y/var).sum()
    Sxx = ((x**2)/var).sum()
    Sxy = ((x*y)/var).sum()
    
    if (no_err):
        Stt = ((x - Sx/S)**2).sum()
        Sty = ((x - Sx/S)*y).sum()
    
    Delta = S*Sxx - Sx**2
    if no_err:
        bb = (1.0/Stt)*Sty
        aa = (Sy - Sx*bb)/S
        var_aa = (1.0/S)*(1.0 + Sx**2/(S*Stt))
        var_bb = (1.0/Stt)
    else:
        bb = ( S*Sxy - Sx*Sy ) / Delta
        aa = ( Sxx*Sy - Sx*Sxy ) / Delta
        var_aa = Sxx/Delta
        var_bb = S / Delta
        
    rms_aa = numpy.sqrt(var_aa)
    rms_bb = numpy.sqrt(var_bb)

    # coefficient of correlation
    if no_err:
        cov_ab = (-Sx/(S*Stt))
        r_ab   = (cov_ab/(rms_aa*rms_bb))
    else:
        cov_ab = (-Sx/Delta)
        r_ab   = -Sx/numpy.sqrt(S*Sxx)

    # get chi_squared
    X2 = (((y - aa - bb*x) / dy)**2).sum()

    return aa, rms_aa, bb, rms_bb, r_ab, X2
  


def robustPolyFit(x, y, order, nbin=3, sigma=3.0, niter=1):

    xNew, yNew = copy.copy(x), copy.copy(y)
    
    # bin and take medians in each bin
    epsilon = 1.0e-6
    xmin, xmax = xNew.min() - epsilon, xNew.max() + epsilon

    step = (xmax-xmin)/(nbin)
    xMeds, yMeds, yErrs = [], [], []
    for i in range(nbin):
        w = numpy.where( (xNew > xmin + i*step) & (xNew <= xmin + (i+1)*step) )
        if len(xNew[w]) == 0:
            continue
        xMeds.append(numpy.median(xNew[w]))
        yMeds.append(numpy.median(yNew[w]))
        yErr = numpy.std(yNew[w])/numpy.sqrt(len(yNew[w]))
        yErrs.append(yErr)
        
    # use these new coords to fit the line ... with sigma clipping if requested
    xNew, yNew, dyNew = numpy.array(xMeds), numpy.array(yMeds), numpy.array(yErrs)

    # if there's only one point in a bin, dy=0 ... use the average error instead
    w0 = numpy.where(dyNew == 0)[0]
    wnot0 = numpy.where(dyNew > 0)[0]
    if len(w0) > 0:
        if len(wnot0) > 0:
            meanError = numpy.mean(dyNew[wnot0])
        # if *all* bins have a single point, then all dy are zero
        # take a blind guess and use stdev of all values we got
        else:
            meanError = numpy.std(yNew)

        # last ... if we have errors of zero, use 1.0
        # sounds silly, but if we fit a line to a difference value
        # and it's the difference of identical values
        # (eg. instFlux and modFlux are loaded with the same value)
        # then the dy could be zero
        if meanError == 0.0:
            meanError = 1.0
        dyNew[w0] = meanError

    
    for i in range(niter):

        #ret = numpy.polyfit(xNew, yNew, order)
        #p = ret[0:2]
        a, da, b, db, rab, x2 = lineFit(xNew, yNew, dyNew)
        
        residuals = yNew - numpy.polyval((a, b), xNew)

        if i == 0:
            mean = numpy.median(residuals)
        else:
            mean = numpy.mean(residuals)
            
        std = numpy.std(residuals)

        if niter > 1:
            w = numpy.where( (numpy.abs(residuals - mean)/std) < sigma )
            xNew = xNew[w]
            yNew = yNew[w]
            dyNew = dyNew[w]
            
    return b, db, a, da



def dictToList(d, withDelete=False):
    out = numpy.array([])
    for k,v in d.items():
        out = numpy.append(out, numpy.array(v))
        if withDelete:
            del d[k]
    return out
