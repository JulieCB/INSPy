
"""
Boring classes that make life fun elsewhere.

"""

import re
import numpy as np
from pyqtgraph import QtGui, QtCore

class AutoButtonWidget(QtGui.QWidget):
    """
    Generates and adds QPushButtons to self.layout() based on and connects 
    clicked signal to methods with names starting with 'pb_'.

    If no layout is set yet, a QVBoxLayout is used.

    """
    
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent=parent)

        if self.layout() is None:
            self.setLayout(QtGui.QVBoxLayout())

        for k in dir(self):
            if k.startswith('pb_'):
                parts = k[3:].split('_')
                parts[0] = parts[0].title()
                title = ' '.join(parts)
                pb = QtGui.QPushButton(title)
                self.layout().addWidget(pb)
                pb.clicked.connect(getattr(self, k))
                setattr(self, '_'+k, pb)


class ChannelSelect(QtGui.QWidget):
    """
    For each channel type, 
    """
    def __init__(self, *args, **kwds):
        QtGui.QWidget.__init__(self, None)
        self.setLayout(QtGui.QHBoxLayout())
        self.selection = {}
        self.callbacks = {}
        self.default_state = kwds.pop('default_state',
                QtCore.Qt.CheckState.Unchecked)
        self.nrow = kwds.pop('nrow', 10)
        doshow = kwds.pop('show', True)
        if len(kwds) == 0 and len(args) > 0:
            for i, arg in enumerate(args):
                kwds['Channel set %d' % (i + 1,)] = arg
        for k, v in kwds.iteritems():
            self.setup_chan_set(k, v)
        if doshow:
            self.show()
    def make_checkbox_callback(self, cb, label):
        def _():
            self.selection[label] = cb.checkState()
        return _
    def setup_chan_set(self, ch_type, ch_labels):
        gb = QtGui.QGroupBox(ch_type)
        lay = QtGui.QVBoxLayout()
        gb.setLayout(lay)
        # setup check box matrix
        l_cb = QtGui.QGridLayout()
        lay.addLayout(l_cb)
        cbs = []
        for i, l in enumerate(ch_labels):
            cb = QtGui.QCheckBox(l)
            cb.setCheckState(self.default_state)
            self.selection[l] = self.default_state
            self.callbacks[l] = self.make_checkbox_callback(cb, l)
            cb.stateChanged.connect(self.callbacks[l])
            cbs.append(cb)
            l_cb.addWidget(cb, i%self.nrow + 1, i/self.nrow)
        setattr(self, 'gb_'+ch_type, gb)
        self.layout().addWidget(gb)
        # setup all none control
        l_pb = QtGui.QHBoxLayout()
        lay.addLayout(l_pb)
        pb_all = QtGui.QPushButton('All')
        @pb_all.clicked.connect
        def set_all():
            for cb in cbs:
               cb.setCheckState(QtCore.Qt.CheckState.Checked)
        l_pb.addWidget(pb_all)
        pb_none = QtGui.QPushButton('None')
        @pb_none.clicked.connect
        def set_none():
            for cb in cbs:
               cb.setCheckState(QtCore.Qt.CheckState.Unchecked)
        l_pb.addWidget(pb_none)
    @property
    def selected(self):
        return [k for k, v in self.selection.iteritems()
                if v == QtCore.Qt.CheckState.Checked]
    def indices(self, labels):
        return np.array([labels.index(k) for k in self.selected])


def seeg_major_minor_ch_labels(labels):
    tick_pos_per_reg = {}
    ch=[(i+0.5,)+re.match(r'([A-Z][a-zA-Z\']+)(\d+)', lab).groups() 
            for i, lab in enumerate(labels)]

    for p, r, i in ch:
        if r not in tick_pos_per_reg:
            tick_pos_per_reg[r]=[]
        tick_pos_per_reg[r].append(p)

    minor_ticks = [(p, i) for p, r, i in ch]
    major_ticks = [(np.max(tps)+0.5, r+'    ') for r, tps in tick_pos_per_reg.iteritems()]

    return major_ticks, minor_ticks
