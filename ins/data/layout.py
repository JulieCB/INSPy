
import glob
import os
import os.path

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))

ftlayouts = {}
for absfilename in glob.glob(os.path.join(here, 'files', 'fieldtrip_layouts', '*')):
    name, ext = os.path.splitext(os.path.basename(absfilename))
    ftlayouts[name] = absfilename

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
        path = ftlayouts[self.path] if self.path in ftlayouts else self.path
        with open(path, 'r') as fd:
            for line in fd.readlines():
                _, x, y, _, _, name = line.strip().split()
                self.labels.append(name.lower())
                self.points.append(map(float, (x, y)))
        print self.labels
        self.points = np.array(self.points)
        self.px, self.py = self.points.T

