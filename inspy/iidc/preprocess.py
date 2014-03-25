# -*- coding: utf-8 -*-

import numpy as np

from scipy import signal

def _filter_signal(b, a, data):
    
    for i, y in enumerate(data):
        data[i] = signal.filtfilt(b, a, y)
        
    return data    

def diff_bp(fs, data, bp=[15.0, 45.0], copy=True):
    """
    Perform a differentiation and then bandpass on signals in data, assuming
    data.shape = (nchan, ntime)
    
    By default, the incoming array is copied before processed, otherwise, 
    use copy=False.
    
    """
    
    if copy:
        data = data.copy()
        
    data = np.diff(data)
    
    lo, hi = bp
    b, a = signal.butter(3, [2*lo/fs, 2*hi/fs], "bandpass")

    return _filter_signal(b, a, data)

def hp(fs, data, hp=2.0, copy=True):
    """
    Perform a high pass on signals in data, assuming data.shape = (nchan, ntime)
    
    By default, the incoming array is copied before processed, otherwise, 
    use copy=False.
    
    """
    
    if copy:
        data = data.copy()
    
    b, a = signal.butter(3, [2*hp/fs], "highpass")

    return _filter_signal(b, a, data)



    
    

    
    
