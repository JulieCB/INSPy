"""
ins.data.tvb
============

Provide access to data used in TVB, e.g. surface vertices & triangles

"""

import os.path
here = os.path.dirname(os.path.abspath(__file__))

import numpy as np
from scipy import io

class cortex_reg13:
    surf_verts = np.load(os.path.join(here, 'files', 'tvb-cortex-reg13-vert.npy'))
    surf_faces = np.load(os.path.join(here, 'files', 'tvb-cortex-reg13-tri.npy'))

_gains = io.loadmat(os.path.join(here, 'files', 'tvb-lead-fields.mat'))

class meg:
    _labels_path = os.path.join(here, 'files', 'tvb-meg-channel-names.mat')
    labels = [l[0] for l in io.loadmat(_labels_path)['meg_names'][0]]
    del l
    gains = _gains['meg']

class eeg:
    _labels_path = os.path.join(here, 'files', 'tvb-eeg-channel-names.mat')
    labels = [l[0] for l in io.loadmat(_labels_path)['eeg_name'][0]]
    del l
    gains = _gains['eeg']

