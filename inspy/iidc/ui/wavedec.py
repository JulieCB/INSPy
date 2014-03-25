
"""
Simple wavelet decomposition viewer components

TODO Don't compute dec all at once but online as Signal view changes
TODO Drop-down menu to choose wavelet form
TODO normalization ?? 

"""

import time
import numpy as np
import pywt
import pyqtgraph as pg

def wtrect(n,h):
    rect = pg.QtCore.QRect()
    rect.setTop(h)
    rect.setBottom(h)
    rect.setLeft(0)
    rect.setRight(n)
    return rect

class Signal(pg.PlotItem):
    def __init__(self, fs, data, parent=None):
        pg.PlotItem.__init__(self, parent)
        self.time = np.r_[:data.size]*1.0/fs
        #HACK
        self.time = np.r_[:int(data.size*1.0/fs+1):1j*data.size]*1.0
        self.data = data
        self.pdi = pg.PlotDataItem(x=self.time, y=self.data)
        self.addItem(self.pdi)
        self.getAxis('left').setWidth(80)

class WaveDec(pg.PlotItem):
    def __init__(self, fs, data, wtype='db1', minfreq=2.0, maxfreq=150.0, parent=None):
        pg.PlotItem.__init__(self, parent)
        self.endtime = (len(data)-1)*1.0/fs
        self.endtime = int(data.size*1.0/fs)
        self.data = data
        self.wtype = wtype
        self.compute_wavedec(wtype)
        self.ticks = []
        #self.wfs = [fs/(2**(len(self.cD)-i)) for i in len(self.cD)]
        for i, cdi in enumerate(self.cD):
            fi = fs/(2**(len(self.cD)-i))
            im = pg.ImageItem(cdi[:, np.newaxis]*cdi.size)#/(2**(len(self.cD)-i)))
            im.setRect(wtrect(self.endtime, len(self.ticks)))
            self.addItem(im)
            self.ticks.append((i+0.5, '%.0f Hz' % (fi, ) if fi > 1 else '%.0f s' % (1/fi,)))
        self.getAxis('left').setTicks([self.ticks])
        self.getAxis('left').setWidth(80)

    def compute_wavedec(self, wtype):
        self.wavelet = pywt.Wavelet(wtype)
        tic = time.time()
        cs = pywt.wavedec(self.data, self.wavelet)
        print 'wavedec required %.2f s' % (time.time() - tic, )
        self.cA, self.cD = cs[0], cs[1:]


class WaveView(pg.GraphicsWindow):

    def __init__(self, fs, signal, wtype='db1', minfreq=2.0, parent=None):
        pg.GraphicsWindow.__init__(self, parent)
        self.fs = fs
        self.signal = signal
        self.p_ts = Signal(fs, signal)
        self.p_ts.showGrid(x=1, y=1)
        self.addItem(self.p_ts)
        self.nextRow()
        self.p_wd = WaveDec(fs, signal, wtype=wtype, minfreq=minfreq)
        self.p_wd.setXLink(self.p_ts)
        self.addItem(self.p_wd)
        self.show()


