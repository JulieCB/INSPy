import sys
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np
from matplotlib import cm


app = QtGui.QApplication([])
w_main = QtGui.QMainWindow()
w_cent = QtGui.QSplitter() #QtGui.QWidget()
w_main.setCentralWidget(w_cent)
#lay = QtGui.QHBoxLayout()
#w_cent.setLayout(lay)

# surface plot
w_surf = gl.GLViewWidget()
w_cent.addWidget(w_surf)
w_surf.setCameraPosition(distance=120)

verts = np.load('/home/duke/tvb-vert.npy')
faces = np.load('/home/duke/tvb-tri.npy')
nc = verts.shape[0]
colors = np.zeros((nc, 4))
colors = cm.hot(np.random.rand(nc))

m_data = gl.MeshData(vertexes=verts, faces=faces, vertexColors=colors)
m_surf = gl.GLMeshItem(meshdata=m_data)
#m_surf.setGLOptions('additive')
w_surf.addItem(m_surf)

# fake data
import h5py
h5 = h5py.File('/home/duke/surfspec.h5')
P = h5['all'][int(sys.argv[1])]
fs = np.r_[:P.shape[-1]*1.0]
lfs = np.log10(fs)

# line plots
split_lines = QtGui.QSplitter(QtCore.Qt.Vertical)
w_cent.addWidget(split_lines)

# spectra
w_spec = pg.PlotWidget()
split_lines.addWidget(w_spec)
lr_select = pg.LinearRegionItem([0., 0.2], bounds=[lfs.min(), lfs.max()], movable=True)
w_spec.addItem(lr_select)

P1, P2, P3 = np.percentile(P, [25, 50, 75], axis=0)
w_spec.plot(np.log10(fs), np.log10(P2), pen=pg.mkPen(width=5))
w_spec.plot(np.log10(fs), np.log10(P1))
w_spec.plot(np.log10(fs), np.log10(P3))

# par plots
split_pars = QtGui.QSplitter()
split_lines.addWidget(split_pars)

p_a = 2.0**np.r_[-19:-9]
p_lg = 2.0**np.r_[-19:-9] # changes slower

w_p_a = pg.PlotWidget()
split_pars.addWidget(w_p_a)
w_p_a.plot(p_a)
il_a = pg.InfiniteLine(pos=0, movable=True, bounds=[0, len(p_a)])
w_p_a.addItem(il_a)

w_p_lg = pg.PlotWidget()
split_pars.addWidget(w_p_lg)
w_p_lg.plot(p_lg)
il_lg = pg.InfiniteLine(pos=0, movable=True, bounds=[0, len(p_lg)])

@il_a.sigDragged.connect
def il_a_dragged(il):
    i = np.round(il.value()).astype(int)
    

@lr_select.sigRegionChanged.connect
def reselect(lr):
    lo, hi = map(lambda v: 10.0**v, lr.getRegion())
    mask = np.c_[fs>lo, fs<hi].all(axis=1)
    mP = np.log10(P[:, mask].mean(axis=1))
    mP = np.clip((mP - mP.mean())/(3*mP.std())+0.5, 0.0, 1.0)
    #mP = (mP - mP.min())/mP.ptp()
    colors[:] = cm.jet(mP)
    m_data.setVertexColors(colors)
    m_surf.meshDataChanged()

lr_select.setRegion((0., 1.))

# show window
w_main.show()

## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
