
"""
Provides a widget for viewing a set of time series.

"""

import numpy as np
import pyqtgraph as pg

class PagingViewBox(pg.ViewBox):
    def mouseMoveEvent(self, ev):
        print ev
    def mouseReleaseEvent(self, ev):
        print self.viewRange() 
        pg.ViewBox.mouseReleaseEvent(self, ev)

class MinutesSecondsAxis(pg.AxisItem):
    def tickStrings(self, values, scale, spacing):
        out = []
        s = 1
        m = 60*s
        for value in map(int, values):
            nm = value / m
            ns = value % m
            out.append('%dm%ds' % (nm, ns) if nm else '%ds' % (ns, ))
        return out

class TimeSeriesWidget(pg.PlotWidget):

    def __init__(self, fs, data, **kwds):
        labels=kwds.pop('labels', [])
        self.yax = pg.AxisItem(orientation='left')
        self.xax = MinutesSecondsAxis(orientation='bottom')
        self.yax.setTicks([[], [(i, l) for i, l in enumerate(labels)]])
        kwds['axisItems'] = {'left': self.yax, 'bottom': self.xax}
        show = kwds.pop('show', True)

        self.pvb = PagingViewBox()
        kwds['viewBox'] = self.pvb

        pg.PlotWidget.__init__(self, **kwds)

        self.setDownsampling(auto=True, mode='peak')
        self.setClipToView(True)
        self.setXRange(0, 10)
        self._fs = fs
        self._data = data
        self._dstd = data.std()
        self.t = np.r_[:data.shape[1]]*1.0/fs
        self.scale = 0.5
        self.lines = []
        for _ in xrange(len(data)):
            pdi = pg.PlotDataItem()
            self.addItem(pdi)
            self.lines.append(pdi)
        self.update_lines()
        self.showGrid(x=True, y=True)
        if show:
            self.show()

    def update_data(self, data):
        self._data = data
        self.update_lines()

    def update_lines(self):
        for i, (line, y) in enumerate(zip(self.lines, self._data)):
            line.setData(self.t, y/self._dstd*self.scale + i, antialias=True)

    def wheelEvent(self, ev):
        self.scale *= 1 + 0.1*(1 if ev.delta()>0 else -1)
        self.update_lines()

    dd = 0.8

    def keyReleaseEvent(self, ev):

        koi = 'Left Right Up Down Minus Plus Home End'
        l, r, u, d, m, p, h, e = [getattr(pg.QtCore.Qt, 'Key_'+s) for s in koi.split(' ')]
        ek = ev.key()

        dx, dy = 0, 0
        vb = self.getPlotItem().getViewBox()
        (xmin, xmax), (ymin, ymax) = vb.viewRange()

        # just arrow keys
        if ek in (l, r):
            dx = -1 if ek==l else 1
        elif ek in (u, d):
            if ev.modifiers() and pg.QtCore.Qt.Key_Control:
                self.scale *= 1 + 0.25*(1 if ek==u else -1)
                return
            dy = -1 if ek==d else 1
        elif ek in (h, e):
            dx = -xmin if ek==h else (self.t[-1] - xmax)
            
        # shift + arrow
        if ev.modifiers() and pg.QtCore.Qt.Key_Shift:
            if ek in (l, r):
                dx = dx*(xmax - xmin)*self.dd
            elif ek in (u, d):
                dy = dy*(ymax - ymin)*self.dd

        vb.translateBy((dx, dy))

def compare(fs, labels, *data, **kwds):
    win = pg.QtGui.QSplitter()
    win.setOrientation(pg.QtCore.Qt.Vertical)
    kwds['show'] = False
    first_ts = None
    for i, datum in enumerate(data):
        ts = TimeSeriesWidget(fs, datum, labels=labels, **kwds)
        win.addWidget(ts)
        if i==0:
            first_ts = ts
        else:
            ts.getPlotItem().setXLink(first_ts.getPlotItem())
    win.show()
    return win

class TimeSeriesWork(object):
    """
    Use parameter tree to detail channels and quick apply filter sets

    """

    def __init__(self):

        mw = pg.QtGui.QWidget()
        lay = pg.QtGui.QGridLayout()
        mw.setLayout(lay)
        lay.addWidget(tsw, 0, 0)
        clay = pg.QtGui.QHBoxLayout()
        lay.addLayout(clay, 1, 0)
        ctrls = [
            ('lb_low', pg.QtGui.QLabel('Low freq (Hz)')),
            ('sb_low' , pg.QtGui.QSpinBox()),
            ('lb_high', pg.QtGui.QLabel('Low freq (Hz)')),
            ('sb_high', pg.QtGui.QSpinBox()),
            ('pb_apply', pg.QtGui.QPushButton('Apply Band-Pass')),
        ]
        for _, v in ctrls:
            clay.addWidget(v)
        ctrls_d = {k:v for k, v in ctrls}
        ctrls_d['sb_low'].setValue(15.0)
        ctrls_d['sb_high'].setValue(45.0)
        mw.show()



