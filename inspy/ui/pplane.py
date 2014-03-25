
"""
Phase plane tool. 

Likely, it will be easier to do this with pyqtgraph's parameter trees.

"""

import pyqtgraph as pg
from pyqtgraph import QtCore as qc, QtGui as qg

class ParSlider(qg.QSlider):
    def sizeHint(self):
        return qc.QSize(200, 20)

class ParDel(qg.QPushButton):
    def sizeHint(self):
        return qc.QSize(20, 20)

class Parameter(qg.QWidget):
    "A parameter"

    def __init__(self, parent=None):
        qg.QWidget.__init__(self, parent=parent)
        self._parent = parent
        self.lay = qg.QHBoxLayout()
        self.setLayout(self.lay)
        self.le_name = qg.QLineEdit('par')
        self.sp_lo = pg.SpinBox(value=-1)
        self.sp_hi = pg.SpinBox(value=1)
        self.slider = ParSlider(qc.Qt.Horizontal)
        self.pb_rm = ParDel('x')
        self.pb_rm.clicked.connect(self.cb_rm)
        self.lay.addWidget(self.pb_rm)
        self.lay.addWidget(self.le_name)
        self.lay.addWidget(self.sp_lo)
        self.lay.addWidget(self.slider)
        self.lay.addWidget(self.sp_hi)

    def cb_rm(self):
        self._parent.lay_pars.removeWidget(self)


class Parameters(qg.QWidget):
    "A set of parameters"

    sig_pars_changed = qc.Signal(object)

    def __init__(self, parent=None):
        qg.QWidget.__init__(self, parent=parent)
        self.lay = qg.QVBoxLayout()
        self.setLayout(self.lay)
        self.lay_pars = qg.QVBoxLayout()
        self.lay_ctrl = qg.QHBoxLayout()
        self.pb_add_par = qg.QPushButton('Add parameter')
        self.pb_del_par = qg.QPushButton('Remove parameter')
        self.pb_save_par = qg.QPushButton('Save')
        self.pb_load_par = qg.QPushButton('Load')
        self.pb_add_par.clicked.connect(self.cb_add_par)
        self.lay_ctrl.addWidget(self.pb_add_par)
        self.lay_ctrl.addWidget(self.pb_del_par)
        self.lay_ctrl.addWidget(self.pb_save_par)
        self.lay_ctrl.addWidget(self.pb_load_par)
        self.lay.addLayout(self.lay_ctrl)
        self.lay.addLayout(self.lay_pars)
        self.pars = []

    def cb_add_par(self):
        par = Parameter()
        self.pars.append(par)
        self.lay_pars.addWidget(par)
        par.le_name.textEdited.connect(self.pars_changed)
        par.sp_lo.sigValueChanged.connect(self.pars_changed)
        par.sp_hi.sigValueChanged.connect(self.pars_changed)
        par.slider.sliderMoved.connect(self.pars_changed)

    def pars_changed(self):
        print 'ok'
        self.sig_pars_changed.emit({})


if __name__ == '__main__':

    app = qg.QApplication([])
    par = Parameters()
    par.show()
    app.exec_()

