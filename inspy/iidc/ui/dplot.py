
"""
Density plots

"""

import pyqtgraph as pg
import numpy as np

def cross_density(data, max_chan=20):
    """
    Compute histograms between all pairs of signals and plot 
    as a matrix of cross-signal densities. 

    """

    gv = pg.GraphicsWindow()
    hd = data

    if len(hd) > max_nchan:
        raise ValueError('too many channels (%d > %d)' % (len(hd), max_nchan))

    for di in hd:
        for dj in hd:
            H, xb, yb = np.histogram2d(di, dj, 100)
            im_ij = pg.ImageItem()
            im_ij.setImage(np.log(H+1))
            pij = gv.addPlot()
            pij.addItem(im_ij)
        gv.nextRow()
    gv.show()
    return gv


def pasta(ep, labels, ch_labels):
    pw = pg.GraphicsWindow()
    for i, ch in enumerate(ch_labels):
        epi = ep[..., i]
        for j, lab in enumerate(unique(labels)):
            epij = epi[labels==lab]
            T = tile(r_[:epij.shape[1]], (epij.shape[0], 1)).ravel()
            Y = epij.ravel()
            if j==0:
                pw.addLabel(ch_labels[i])
            q5, q95 = percentile(Y, [0.1, 99.9])
            ny = epij.shape[0]/3
            H, _, _ = np.histogram2d(T, clip(Y, q5, q95), (epij.shape[1], ny))
            ii = pg.ImageItem()
            ii.setImage(log(H+1))
            pyi = pw.addPlot()
            pyi.addItem(ii)
            lo, hi = epij.min(), epij.max()
            pens = [pg.mkPen('r', width=i) for i in [1, 2, 4, 2, 1]]
            for mode, pen in zip(percentile(epij, [5, 25, 50, 75, 95], axis=0), pens):
                mode *= ny/(q95 - q5)
                mode += ny*(epij.mean() - q5)/(q95 - q5)
                modeline = pg.PlotDataItem(y=mode, pen=pen)
                pyi.addItem(modeline)
            pyi.addLine(y=(epij.mean() - q5)*ny/(q95 - q5), pen=pg.mkPen('k', style=pg.QtCore.Qt.DashLine, width=2))
            pyi.showAxis('bottom', False)
            pyi.showAxis('left', False)
        pw.nextRow()
    pw.show()

