
"""
iidc.detect
===========


"""

import functools

import numpy as np
from scipy import stats, signal, optimize

from . import stat, util


def lfdr_detect(fs, y, q=1e-3, minel=0.0, maxel=0.5, rf=0.5, ht=False, dc=1, info=False):
    """
    Perform detection on continuous signal `y`.

    This function filters the signal, estimates the LFDR (see `lfdr`), and
    identifies where `y` remains "unlikely" (LFDR < `q`), for a period of
    time longer than `minel` and shorter than `maxel`, with a refractory
    period of `rf`.
    
    `y` may be an array of signals, i.e. shape == (nchan, ntime) in which 
    case detection is performed on each channel, and a list is returned.


    Parameters
    ----------

    fs : float
        Sampling frequency of `y`
    y : array
        Continuously sampled signal on which to perform detection
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

    fs = float(fs)
    
    t = np.r_[:len(y)]/fs
    
    # analyze distribution
    xb, f, cf, fdr, llx = stat.lfdr(y, dc=dc, doplot=False)
    
    # fdr tarnsform hilbert 
    if ht:
        hy = abs(signal.hilbert(y))
        llh = np.interp(hy, xb, np.log(fdr))
    else:
        hy, llh = y, llx
    
    # generate events    
    ev = np.c_[np.isfinite(llh), llh < np.log(q)].all(axis=1)
    ev[~np.isfinite(llh)] = True
    e0, = np.argwhere(np.c_[~ev[:-1],  ev[1:]].all(axis=1)).T
    e1, = np.argwhere(np.c_[ ev[:-1], ~ev[1:]].all(axis=1)).T
    
    if len(e0) == 0 or len(e1) == 0:
        return locals() if info else []
    
    # boundaries
    if e1[0] < e0[0]:
        e0 = np.r_[0, e0]
    if e0[-1] > e1[-1]:
        e1 = np.r_[e1, llh.shape[0]-1]
    
    # compute event length and long-enough mask
    el = np.diff(np.c_[e0, e1])[:, 0]/fs
    le = np.c_[el > minel, el < maxel].all(axis=1)
    
    # pull out remaining peaks & align them
    peaks = []
    for i, (e0i, e1i) in enumerate(zip(e0[le], e1[le])):
        hi = hy[e0i:e1i]
        peaks.append((hi.max(), t[e0i] + np.argmax(hi)/fs))
    
    if len(peaks) == 0:
        return locals() if info else []
    
    # mask refractory period
    ph, pt = np.array(peaks).T
    nonmask = []
    for i, (phi, pti) in enumerate(peaks):
        mask = np.c_[pt > pti - rf, pt < pti + rf].all(axis=1)
        if (phi >= ph[mask]).all():
            nonmask.append(i)
    
    nfpeak = np.array(peaks)[np.array(nonmask), 1]
    
    return locals() if info else nfpeak


def _lfdr_detect_single(arg):
    fs, yi, kwds = arg
    try:
        out = lfdr_detect(fs, yi, **kwds)
    except Exception as out:
        pass
    return out


def batch_lfdr_detect(fs, y, n_jobs=1, **kwds):
    """
    Perform lfdr detection on each row of `y` using `n_jobs` processes. `kwds`
    will be provided to the detection function. 

    """

    jobs = [(fs, yi, kwds) for yi in y]
    with util.mpool(n_jobs) as p:
        results = p.map(_lfdr_detect_single, jobs)
    return results
 
