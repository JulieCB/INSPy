import h5py
import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl

from pyqtgraph.Qt import QtCore, QtGui

from scipy import interpolate, io
from matplotlib import cm

class Layout(object):

    def __init__(self, path):
        self.path = path
        self.read_layout()

    def find_data_order(self, data_labels):
        ix = []
        dl = [l.lower() for l in data_labels]
        for l in self.labels:
            l = l.lower()
            ix.append(dl.index(l) if l in dl else -1)
        return np.r_[:len(self.labels)], np.array(ix)

class FTLayout(Layout):
    def read_layout(self):
        self.labels = []
        self.points = []
        with open(self.path, 'r') as fd:
            for line in fd.readlines():
                _, x, y, _, _, name = line.strip().split()
                self.labels.append(name.lower())
                self.points.append(map(float, (x, y)))
        print self.labels
        self.points = np.array(self.points)
        self.px, self.py = self.points.T

class GridFwd(object):

    def _set_ngrid(self, ngrid):
        if not 2 < ngrid < 500:
            raise AttributeError('unexpected ngrid %d' % (ngrid,))
        self._ngrid = ngrid
        l = self.lay
        self.X, self.Y = np.mgrid[l.px.min() : l.px.max() : 1j*ngrid, l.py.min() : l.py.max() : 1j*ngrid]
        self.XY = np.c_[self.X.flat, self.Y.flat]

    ngrid = property(lambda s: s._ngrid, _set_ngrid)

    def __init__(self, lay, gain, labels, ngrid=100):
        self.lay = lay
        self.gain = gain
        self.labels = labels
        self.lix, self.dix = lay.find_data_order(labels)
        print self.dix.shape, np.isfinite(gain).all(axis=1).shape
        data_mask = self.dix >= 0
        print 'data_mask', data_mask.sum()
        gain_mask = np.zeros(self.dix.shape)
        gain_mask[data_mask] = np.isfinite(gain).all(axis=1)
        print 'gain_mask', gain_mask.sum()
        self.mask = data_mask #np.c_[data_mask, gain_mask].all(axis=1)
        self.ngrid = ngrid

    def gridded(self, surf):
        fwd = np.zeros((self.lay.points.shape[0],))
        fwd_ = self.gain.dot(surf)
        fwd[self.mask] = fwd_[self.dix[self.mask]]
        return interpolate.griddata(self.lay.points, fwd, self.XY, fill_value=0.0).reshape(self.X.shape)

class Viz(object):

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

        self.surf_verts = np.load('/home/duke/tvb-vert.npy')
        self.surf_faces = np.load('/home/duke/tvb-tri.npy')
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
        self.main_split = QtGui.QSplitter(QtCore.Qt.Vertical)
        self.main_split.addWidget(self.fwd_surf_split)
        self.main_split.addWidget(self.time)
        self.main_split.setSizes([0.9, 0.1])
        
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
        

if __name__ == '__main__':

    import sys
    import h5py

    # open h5 file
    h5 = h5py.File(sys.argv[1])
    
    lay_eeg = FTLayout('/home/duke/tools/fieldtrip/template/layout/elec1010.lay')
    lay_meg = FTLayout('/home/duke/tools/fieldtrip/template/layout/4D248.lay')

    meg_labels = io.loadmat('/data/local/tvb-meg-channel-names.mat')['meg_names']
    meg_labels = [l[0] for l in meg_labels[0]]

    eeg_labels = io.loadmat('/data/local/tvb-eeg-channel-names.mat')['eeg_name']
    eeg_labels = [l[0] for l in eeg_labels[0]]

    fwds = io.loadmat('/data/local/tvb-lead-fields.mat')
    
    fwd_eeg = GridFwd(lay_eeg, fwds['eeg'], eeg_labels)
    fwd_meg = GridFwd(lay_meg, fwds['meg'], meg_labels)

    app = QtGui.QApplication([])
    win = QtGui.QMainWindow()
    viz = Viz(h5, 512., fwd_eeg, fwd_meg)
    win.setCentralWidget(viz.main_split)
    print 'showing!'
    win.show()
    QtGui.QApplication.instance().exec_()
