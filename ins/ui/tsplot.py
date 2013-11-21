
"""
TSPlot
======

Provides class for efficient plotting of large time series data.

.. todo:: per signal scaling
.. todo:: integration with detection
.. todo:: handle showing markers from file  (~?)
.. todo:: marker annotations (easier)
.. todo:: bipolar (move to general util)
.. todo:: ETPlot.frommarkers(tsp) regroups data per chan around marker
.. todo:: absolute time as well as 
.. todo:: using blit to speed  up vertical cursor on supported backends

"""

import time

from pylab import *
from numpy import *
from matplotlib.widgets import Cursor
from matplotlib.ticker import FuncFormatter
from scipy import signal

from spectrum.lpc import lpc

from ins import util

__all__ = ['time_format_minor_units', 'TSPlot']

time_format_minor_units = True

@FuncFormatter
def TimeFormatter(s, _):
    s, ms = divmod(s, 1)
    m, s = divmod(s, 60)
    h, m = divmod(m, 60)
    d, h = divmod(h, 24)
    w, d = divmod(d, 7)

    mu = time_format_minor_units

    if w:
        return '%dw%d' % (w, d) if mu else '%d w' % (w, )
    elif d:
        return '%dd%02d' % (d, h) if mu else '%d d' % (d, )
    elif h:
        return '%dh%02d' % (h, m) if mu else '%d h' % (h, )
    elif m:
        return '%dm%02d' % (m, s) if mu else '%d m' % (m, )
    else:
        return '%ds%03d' % (s, int(1000*ms)) if mu else '%d s' % (s, )


class TSPlot(object):
    """
    Matplotlib plot to handle time series in a useful way, interactively.

    """

    # TODO burden should be on signal sets to implement a to_tsplot method
    @classmethod
    def fromeeg(cls, eeg):
        obj = cls(eeg.data, eeg.labels, eeg.srate)
        obj.eeg = eeg
        return obj

    @classmethod
    def fromeegmat(cls, filename):
        import eegdata, scipy.io
        return cls.fromeeg(eegdata.EEG(scipy.io.loadmat(filename)['EEG']))

    @classmethod
    def frompicked(cls, self):
        if self.picked_lines:
            labels = []
            signals = empty((len(self.picked_lines), self.ys.shape[1]), self.ys.dtype)
            ii = 0
            for i, l in enumerate(self.signal_lines):
                if l in self.picked_lines:
                    labels.append(self.labels[i])
                    signals[ii, :] = self.ys[i]
                    ii += 1
            return TSPlot(signals, labels, self.fs, scale=self.scale)
        else:
            return self

    @classmethod
    def fromdata(cls, set):
        return cls(set.data.values.T, list(set.chan['label'].values), set.fs)

    @classmethod
    def fromtsp(cls, self):
        return cls(self.ys, self.labels, self.fs, scale=self.scale)

    def topicked(self):
        return TSPlot.frompicked(self)

    whiten = False
    persignorm = True
    bandpass = (1.0, 300.0)

    def __init__(self, ys, labels=None, fs=None, period=None, t0=0.0, ax=None, scale=None, markers=[]):

        self.ys = ys
        self.nchan = ys.shape[0]
        self.nt = ys.shape[-1]
        self.labels = labels if labels is not None else map(str, range(self.nchan))
        self.markers = []
        if fs is not None:
            self.fs = fs
        elif period is not None: 
            self.fs = 1.0/period
        else:
            self.fs = 1.0
        self.ts = t0 + r_[:ys.shape[-1]]*1.0/self.fs

        self.scale = scale or 1.0/ys.std()
        self.rescale = 1.3

        self.picked_lines = []
        self.silenced_lines = []

        # setup figues, axes & lines
        if ax is None:
            self.fig = figure()
            self.ax  = subplot(111)
        else:
            self.fig = ax.figure
            self.ax = ax
        sl = slice(0, self.fs)
        tsl = self.ts[sl]
        ysl = self.ys[:, sl]
        self.signal_lines = []
        for i, ysli in enumerate(ysl):
            line, = self.ax.plot(tsl, ysli*self.scale + i, 'k', picker=4)
            self.signal_lines.append(line)
        self.ax.set_yticks(r_[:self.nchan]) 
        self.ax.set_yticklabels(self.labels)
        self.ax.grid(True)
        self.cursor = Cursor(self.ax, color='y')
        self.cursor.horizOn = False
        self.ax.set_ylim(0, 40 if self.nchan > 40 else self.nchan)
        self.ax.set_xlim(0, 1.0)
        self.ax.xaxis.set_major_formatter(TimeFormatter)
        self.render()
        self.fig.tight_layout()

        # callbacks
        for name in dir(self):
            if name.startswith('cb_'):
                self.fig.canvas.mpl_connect(name[3:], getattr(self, name))

    def cb_resize_event(self, *args):
        self.fig.tight_layout()
        self.fig.canvas.draw_idle()

    fastrender = False
    def render(self):

        time = self.ts
        ys = self.ys

        # forward & inverse data-display transforms
        ftr = self.ax.transData.transform
        itr = self.ax.transData.inverted().transform 

        # obtain axis limits in data space
        t0, t1, y0, y1 = self.ax.axis()

        # obtain time limits in display space
        (p0, py), (p1, _) = ftr([(t0, y0), (t1, y0)])

        # limits are outer integer pixels, inclusive
        p0 = int(floor(p0))
        p1 = int(ceil(p1))

        # obtain pixel border locations in data space
        ts = itr([(px, py) for px in r_[p0:p1+1]])[:, 0]

        # throw away what's not in data
        ts = ts[ts>time[0]]
        ts = ts[ts<time[-1]]

        # indices of time points in view
        i0 = where(time < ts[ 0])[0]
        i1 = where(time > ts[-1])[0]
        i0 = i0[-1] if len(i0) else 0
        i1 = i1[ 0] if len(i1) else len(time)

        # slice data source by what's in view
        tsl = time[i0:i1]
        ysl = ys[:, i0:i1].copy()

        # mark points if visible
        marker = (3,3,3) if (p1 - p0) > 10*(i1 - i0) else None

        # preprocessing needs to be done before decimation
        bplo, bphi = self.bandpass
        bpb, bpa = signal.butter(2, [2.0*bplo/self.fs, 2.0*bphi/self.fs], 'pass')
        for i, yi in enumerate(ysl):
            if i > y0 and i < y1:
                if self.whiten:
                    dy = diff(yi)
                    dy = r_[dy[0], dy]
                    yi = signal.filtfilt(lpc(dy, 10)[0], [1.0], dy)
                ysl[i] = signal.filtfilt(bpb, bpa, yi)
                if self.persignorm:
                    ysl[i] = (ysl[i] - ysl[i].mean())/ysl[i].var()

        # if we have at least two points per pixel, we build data
        # from min & max per pixel
        ds = (i1 - i0)/(p1 - p0)
        if ds:
            if self.fastrender:
                tsl = tsl[::ds]
                ysl = ysl[:, ::ds]
            else:
                ii1 = ysl.shape[1]/ds*ds
                tsl_ = tsl[:ii1]
                ysl_ = ysl[:, :ii1].reshape((ysl.shape[0], -1, ds))
                ysl_ = array([ysl_.min(axis=-1), ysl_.max(axis=-1)])
                ysl_ = transpose(ysl_, (1, 2, 0)).reshape((ysl.shape[0], -1))

                tsl = tsl_[0] + r_[:ysl_.shape[1]]*0.5*(tsl_[ds] - tsl_[0])
                ysl = ysl_

        ## this method takes too long!
        if False and 2*(p1 - p0) < (i1 - i0):
            tsl_ = zeros((             2 * (p1 - p0), ), float32)
            ysl_ = zeros((ys.shape[0], 2 * (p1 - p0)), float32)
            # for each pixel
            for it in range(len(ts) - 1):
                # find data index inside pixel
                ii0 = where(tsl > ts[it  ])[0]
                ii1 = where(tsl < ts[it+1])[0]
                ii0 = ii0[ 0] if len(ii0) else 0
                ii1 = ii1[-1] if len(ii1) else len(tsl)-1
                # slice data
                sl = slice(ii0, ii1)
                _t0, _t1 = tsl[ii0], tsl[ii1]
                dt = _t1 - _t0
                ysl__ = ysl[:, sl]
                # always store t, and y, if we have data
                tsl_[it+0] = _t0 + 1.0*dt/4
                tsl_[it+1] = _t0 + 3.0*dt/4
                if ysl__.shape[1]:
                    ysl_[:, it+0] = ysl__.min(axis=1)
                    ysl_[:, it+1] = ysl__.max(axis=1)
                else:
                    break

            tsl = tsl_[1:it-1]
            ysl = ysl_[:, 1:it-1]

        for offset, (line, datai) in enumerate(zip(self.signal_lines, ysl)):
            line.set_marker(marker)
            line.set_color('r' if line in self.picked_lines else 'k')
            if offset > y0 and offset < y1:
                if line in self.silenced_lines:
                    line.set_linestyle('--')
                    line.set_xdata(tsl[[0, -1]])
                    line.set_ydata(array([offset, offset]))
                else:
                    line.set_linestyle('-')
                    line.set_xdata(tsl)
                    line.set_ydata(datai*self.scale + offset)
            else:
                line.set_xdata(array([]))
                line.set_ydata(array([]))

    def cb_scroll_event(self, ev):
        rescale = self.rescale if ev.button == 'up' else 1.0/self.rescale
        self.scale *= rescale
        for offset, line in enumerate(self.signal_lines):
            line.set_ydata((line.get_ydata() - offset) * rescale + offset)
        self.fig.canvas.draw_idle()


    _doing_cb_button_release_event = False

    def add_marker(self, time, **kwds):
        if type(time) in (int, long, int32, int64):
            time = self.ts[time]
        axvline = self.ax.axvline(time, color=kwds.get('color', 'r'), alpha=0.5, picker=4)
        mark = {'time': time, 'axvline': axvline}
        mark.update(kwds)
        self.markers.append(mark)
        self.fig.canvas.draw_idle()

    def add_markers(self, times):
        [self.add_marker(t) for t in times]

    def add_marker_set(self, mset):
        """
        Add markers to each signal depending on tags given.

        .. todo:: Set up a standard marker set interface.
        .. todo:: Set up standard interface for signal/channel
        .. todo:: This doesn't adhere to the add_marker thing

        """

        ys = []
        for tag, time in zip(mset.tags, mset.times):
            for i, ch in enumerate(self.labels):
                if util.chan_intersect(tag, ch):
                    ys.append((time, i-0.5, i+0.5))
        ys = array(ys)
        x, y0, y1 = ys.T
        self.ax.vlines(x, y0, y1, 'r')


    def remove_marker(self, key, val):
        marks = self.markers[:]
        for mark in marks:
            if mark[key] == val:
                mark['axvline'].remove()
                self.markers.remove(mark)
        self.fig.canvas.draw_idle()

    def dispatch_mouse_event(self, prefix, ev, mev=None):
        mev = mev or ev
        b = mev.button
        key = mev.key or ''
        name = '%s_%s_%s' % (prefix, key.replace('+', '_'), b) # e.g. pick__1, pick_ctrl_shift_E_3
        fn = getattr(self, name, None)
        if fn:
            print 'dispatching %s' % (name,)
            return fn(ev)
 
    def cb_button_press_event(self, ev):
        self.dispatch_mouse_event('press', ev)

    def cb_pick_event(self, ev):
        self.dispatch_mouse_event('pick', ev, ev.mouseevent)

    def cb_button_release_event(self, ev, clock=time.time):
        if not self._doing_cb_button_release_event:
            self._doing_cb_button_release_event = True
            self.render()
            self.fig.canvas.draw()
            self._doing_cb_button_release_event = False

    def pick_control_1(self, ev):
        if ev.artist in self.picked_lines:
            self.picked_lines.remove(ev.artist)
        elif ev.artist in self.signal_lines:
                self.picked_lines.append(ev.artist)

    def pick_control_3(self, ev):
        if ev.artist in self.silenced_lines:
            self.silenced_lines.remove(ev.artist)
        elif ev.artist in self.signal_lines:
                self.silenced_lines.append(ev.artist)
    
    def press_shift_1(self, ev):
        self.add_marker(ev.xdata)

    def pick_shift_3(self, ev):
        self.remove_marker('axvline', ev.artist)

    def pick_A_1(self, ev):
        mev = ev.mouseevent
        self.annotxy = mev.xdata, mev.ydata

    def annotmark(self, text, xy=None):
        self.ax.annotate(text, 
            xy=xy or self.annotxy, 
            xytext=(15, 10), 
            bbox={'boxstyle': 'round,pad=0.2', 
                  'fc': 'yellow', 
                  'alpha':0.5}, 
            textcoords='offset points')
        self.fig.canvas.draw_idle()

    moveprop = 0.25
    def cb_key_release_event(self, ev):
        if ev.key in 'mn':
            xlo, xhi = self.ax.get_xlim()
            dx = (xhi - xlo)*self.moveprop*(1 if ev.key=='m' else -1)
            self.ax.set_xlim([xlo + dx, xhi + dx])
            self.cb_button_release_event(ev)

    def view_all(self):
        self.ax.axis([self.ts[0], self.ts[-1], -1, self.ys.shape[0]+1])

    def pick_marked_signals(self):
        """
        .. todo:: Implement me

        """

        raise NotImplementedError

    @property
    def picked_mask(self):
        return array([l in self.picked_lines for l in self.signal_lines])

    @property
    def picked_labels(self):
        return list(array(self.labels)[self.picked_mask])

    @property
    def silenced_mask(self):
        return array([l in self.silenced_lines for l in self.signal_lines])

    @property
    def silenced_labels(self):
        return list(array(self.labels)[self.silenced_mask])






