
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.widgets as mw

def quick(t, y, tlo=0., thi=5., ylo=0, yhi=10, labels=[]):
    tsl = np.c_[t>tlo, t<thi].all(axis=1)
    y_ = y[ylo:yhi, tsl]
    pl.plot(t[tsl], y_.T/y_.ptp(axis=1) + np.r_[:len(y_)], 'k')
    pl.grid(True)

