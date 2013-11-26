"""
ins.detection
=============

Functions to perform automatic detection of events on continuous data.


"""

from numpy import *
from scipy import stats, signal, optimize

def lfdr(x, nbins=50, dc=1, doplot=True):
    """
    Computes density statistics on `x`, where a "null" hypothesis distribution H0,
    here a Gaussian distribution, is fit to the center of the data, and this
    provides an estimation of how likely each bin comes from the Gaussian
    distribution, a measure called the local false discovery rate.

    Decimating the data by specifying `dc` as an integer greater than 1 
    can accelerate the density estimation and safer on large datasets 
    where said estimation takes more time.

    Parameters
    ----------
    x : array
        Data to analyze
    dc : int
        Decimate data, accelerating analysis
    nbins : int or None, optional
        Number of points at which to evaluate the density

    Returns
    -------
    xb : array
        Points at which the densities are evaluated
    f : array
        Density of `x`
    cf : array
        Estimated "center density" of H0
    fdr : array
        Estimated false discovery rate 
    llx : array
        Log lfdr for each datapoint in `x`

    """
    
    k = stats.gaussian_kde(x[::dc])
    xb = r_[x.min() : x.max() : 1j*nbins]
    f = k(xb)
    f /= f.sum()
    
    # TODO use argsort to avoid evaluating f on linspace twice
    dxb = interp(r_[0.0:1.0:1j*nbins], cumsum(f), xb + (xb[1] - xb[0])/2.0)
    xb = unique(r_[xb, dxb])
    xb.sort()                
    f = k(xb)
    f /= f.sum()

    # estimate center sub-density
    def err(par):
        mu, sigma, alo, ahi = par
        sl = slice(argmin(abs(xb - alo)), argmin(abs(xb - ahi)))
        f0 = stats.norm.pdf(xb[sl], loc=mu, scale=sigma)
        f1 = f[sl]
        f0 = f0/f0.max()*f1.max()
        return sum((f0 - f1)**2)/sum(f1**2) - f1.sum()
    
    mu0, sig0 = x.mean(), x.std()
    mu, sig, _, _ = optimize.fmin(err, (mu0, sig0, mu0-sig0, mu0+sig0), disp=0)
    cf = stats.norm.pdf(xb, loc=mu, scale=sig)
    cf = cf/cf.max()*f.max()
        
    # compute lfdr & transform data to log-lfdr
    fdr = clip(cf/f, 0.0, 1.0)
    llx = interp(x, xb, log(fdr))

    if doplot:
        import pylab as pl
        pl.semilogy(xb, f, 'k')
        pl.semilogy(xb, cf, 'k--')
        pl.semilogy(xb, fdr, 'k.')
        pl.grid(True)
        pl.ylim([1e-5, 1.0])
        plevels=[0.2, 0.1, 0.05, 0.01, 0.001]
        pl.yticks(plevels, map(str, plevels))
        
    return xb, f, cf, fdr, llx


def detect(t, y, flo=15., fhi=40., q=1e-5, minel=0.05, maxel=0.5, rf=0.5, info=False):
    """
    Perform detection on continuous signal `y` as function of time `t`. 

    This function filters the signal, estimates the LFDR (see `lfdr`), and
    identifies where `y` remains "unlikely" (LFDR < `q`), for a period of
    time longer than `minel` and shorter than `maxel`, with a refractory
    period of `rf`.


    Parameters
    ----------

    t : array
        Sequence of times corresponding to the samples in `y`, 
            expected to be equally spaced.
    y : array
        Continuously sampled signal on which to perform detection
    flo : float
        Low frequency of the pass band
    fhi : float
        High frequency of the pass band
    q : float
        Threshold on LFDR
    minel : float
        Minimum event length
    maxel : float
        Maximum event length
    rf : float
        Refractory period of events
    info : bool
        Whether to return all information in dictionary or just 
            event times

    Returns
    -------

    ret : dict or array
        An array of event peak times, or if `info` is `True`, a
            dictionary of the local variables.

    """

    fs = 1/(t[1] - t[0])
    
    # filter the temporally diff
    b, a = signal.butter(3, [2*flo/fs, 2*fhi/fs], 'pass')
    fy = signal.filtfilt(b, a, diff(y))
    
    # analyze distribution
    xb, f, cf, fdr, llx = lfdr(fy, doplot=False)
    
    # fdr tarnsform hilbert 
    hy = abs(signal.hilbert(fy))
    llh = interp(hy, xb, log(fdr))
    
    # generate events    
    ev = c_[isfinite(llh), llh < log(q)].all(axis=1)
    ev[~isfinite(llh)] = True
    e0, = argwhere(c_[~ev[:-1], ev[1:]].all(axis=1)).T
    e1, = argwhere(c_[ev[:-1], ~ev[1:]].all(axis=1)).T
    
    if len(e0) == 0 or len(e1) == 0:
        return locals() if info else []
    
    # boundaries
    if e1[0] < e0[0]:
        e0 = r_[0, e0]
    if e0[-1] > e1[-1]:
        e1 = r_[e1, llh.shape[0]-1]
    
    # compute event length and long-enough mask
    el = diff(c_[e0, e1])[:, 0]/fs
    le = c_[el > minel, el < maxel].all(axis=1)
    
    # pull out remaining peaks & align them
    peaks = []
    for i, (e0i, e1i) in enumerate(zip(e0[le], e1[le])):
        hi = hy[e0i:e1i]
        peaks.append((hi.max(), t[e0i] + argmax(hi)/fs))
    
    if len(peaks) == 0:
        return locals() if info else []
    
    # mask refractory period
    ph, pt = array(peaks).T
    nonmask = []
    for i, (phi, pti) in enumerate(peaks):
        mask = c_[pt > pti - rf, pt < pti + rf].all(axis=1)
        if (phi >= ph[mask]).all():
            nonmask.append(i)
    
    nfpeak = array(peaks)[array(nonmask), 1]
    
    return locals() if info else nfpeak


def review_detection(t, y, info):

    peaks = info['nfpeak']
    mark = lambda t: axvline(t, color='r', alpha=0.5)

    fig = figure()
    ax = subplot(311)
    plot(ica.t, ica.ys[0], 'k')
    xlim([0, 10]), grid(True), xticks(r_[:int(ica.t.max())], [])
    map(mark, peaks)

    subplot(312, sharex=ax)
    plot(ica.t[:-1], info['fy'], 'k')
    plot(ica.t[:-1], info['hy'], 'k')
    xlim([0, 10]), grid(True), xticks(r_[:int(ica.t.max())], [])
    map(mark, peaks)

    subplot(313, sharex=ax)
    plot(ica.t[:-1], info['llh'], 'k')
    xt = r_[:int(ica.t.max())]
    grid(True), xticks(xt, map(str, xt)), xlim([0, 10])
    map(mark, peaks)

    tight_layout()
