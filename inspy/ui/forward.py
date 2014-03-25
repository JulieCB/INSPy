import numpy as np
from scipy import interpolate

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
        #print self.dix.shape, np.isfinite(gain).all(axis=1).shape
        data_mask = self.dix >= 0
        #print 'data_mask', data_mask.sum()
        #gain_mask = np.zeros(self.dix.shape)
        #gain_mask[data_mask] = np.isfinite(gain).all(axis=1)
        #print 'gain_mask', gain_mask.sum()
        self.mask = data_mask #np.c_[data_mask, gain_mask].all(axis=1)
        self.ngrid = ngrid

    def gridded(self, surf):
        fwd = np.zeros((self.lay.points.shape[0],))
        fwd_ = self.gain.dot(surf)
        fwd[self.mask] = fwd_[self.dix[self.mask]]
        return interpolate.griddata(self.lay.points, fwd, self.XY, fill_value=0.0).reshape(self.X.shape)

