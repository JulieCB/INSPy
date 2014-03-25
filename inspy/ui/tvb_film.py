import h5py
import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl

from pyqtgraph.Qt import QtCore, QtGui

from scipy import interpolate, io
from matplotlib import cm

from ins.data import tvb
from ins.data.layout import FTLayout
from ins.ui.forward import GridFwd

class SurfaceTime(object):

    def __init__(self, surf_h5, fs, fwd_eeg, fwd_meg):
        self.fs = fs
        self.ntime = surf_h5['surf'].shape[0]
        self.surf_h5 = surf_h5
        self.fwd_eeg = fwd_eeg
        self.fwd_meg = fwd_meg
        print 'calculating surf min max'
        self.surf_min = surf_h5['surf'][2000::10, 0, ::10].min()
        self.surf_max = surf_h5['surf'][2000::10, 0, ::10].max()
        print 'done'
        self.create_widgets()
        self.setup_layout()

    def create_widgets(self):

        # forward solution image plots
        self.meg = pg.ImageView()
        self.meg.setLevels(-8e-5, 8e-5)
        self.eeg = pg.ImageView()
        self.eeg.setLevels(-1e3, 1e3)

        # surface mesh item
        self.setup_surface()

        # time controls
        self.time = pg.PlotWidget()
        self.time.showAxis('left', False)
        self.time.setLabel('bottom', 'time', 's')
        self.il_time = pg.InfiniteLine(movable=True)
        self.il_time.sigPositionChanged.connect(self.time_reset)
        self.time.addItem(self.il_time)

    def setup_surface(self):
        self.surf_vw = gl.GLViewWidget()
        self.surf_vw.setCameraPosition(distance=120)

        self.surf_verts = tvb.cortex_reg13.surf_verts
        self.surf_faces = tvb.cortex_reg13.surf_faces
        self.surf_nv = self.surf_verts.shape[0]
        self.surf_colors = np.zeros((self.surf_nv, 4))

        self.surf_data = gl.MeshData(vertexes=self.surf_verts, faces=self.surf_faces, 
                                     vertexColors=self.surf_colors)
        self.surf_mesh = gl.GLMeshItem(meshdata=self.surf_data)
        #m_surf.setGLOptions('additive')
        self.surf_vw.addItem(self.surf_mesh)

    def setup_layout(self):
        self.fwd_split = QtGui.QSplitter(QtCore.Qt.Vertical)
        self.fwd_split.addWidget(self.eeg)
        self.fwd_split.addWidget(self.meg)
        self.fwd_surf_split = QtGui.QSplitter()
        self.fwd_surf_split.addWidget(self.fwd_split)
        self.fwd_surf_split.addWidget(self.surf_vw)
        self.fwd_surf_split.setSizes([400, 400])
        self.main_split = QtGui.QSplitter(QtCore.Qt.Vertical)
        self.main_split.addWidget(self.fwd_surf_split)
        self.main_split.addWidget(self.time)
        self.main_split.setSizes([400, 50])
        
        # TODO self.main_split set weights?

    def time_reset(self, il):
        self.render_at_time(il.value())

    def render_at_time(self, t):
        ti = int(self.fs * t)
        if ti < 0 or ti >= self.ntime:
            ti = 0 if ti < 0 else self.ntime - 1
        st = self.surf_h5['surf'][ti, 0]
        self.meg.setImage(self.fwd_meg.gridded(st), autoLevels=False, levels=[-8e-5, 8e-5])
        self.eeg.setImage(self.fwd_eeg.gridded(st)[::-1][:, ::-1], autoLevels=False, levels=[-1e3, 1e3])
        st = (st - self.surf_min)/(self.surf_max - self.surf_min)
        self.surf_colors[:] = cm.jet(st)
        self.surf_data.setVertexColors(self.surf_colors)
        self.surf_mesh.meshDataChanged()

    def advance(self):
        vb = self.time.getViewBox()
        (xl, xh), (yl, yh) = vb.viewRange()
        newt = self.il_time.value() + (xh - xl)/1e3
        self.il_time.setValue(newt)
        if (newt - xl)/(xh - xl) > 0.8:
            dt = 0.6*(xh - xl)
            vb.setRange(xRange=(xl + dt, xh + dt))
        
class SurfaceWindow(QtGui.QMainWindow):

    _play_status = False

    def __init__(self, *surfargs):
        super(SurfaceWindow, self).__init__()
        self.surface = SurfaceTime(*surfargs)
        self.setCentralWidget(self.surface.main_split)
        self.play_timer = QtCore.QTimer()
        self.play_timer.timeout.connect(self.maybe_update)
        self.play_timer.setInterval(1000/15)
        self.play_timer.start()
        
    def keyPressEvent(self, ev):
        if ev.key() == 32:
            print 'reversing polarity!'
            self._play_status = not self._play_status

    def maybe_update(self):
        if self._play_status:
            self.surface.advance()

if __name__ == '__main__':

    import sys
    import h5py

    # open h5 file
    h5 = h5py.File(sys.argv[1])
    
    fwd_eeg = GridFwd(FTLayout('elec1010'), tvb.eeg.gains, tvb.eeg.labels)
    fwd_meg = GridFwd(FTLayout('4D248'), tvb.meg.gains, tvb.meg.labels)

    app = QtGui.QApplication([])
    win = SurfaceWindow(h5, 512., fwd_eeg, fwd_meg)
    win.show()
    print """
Instructions:

- drag the time line around to change the time point viewed
- press spacebar to play back the time series
- zoom in/out of timeline to move more/less quickly 
"""
    QtGui.QApplication.instance().exec_()
