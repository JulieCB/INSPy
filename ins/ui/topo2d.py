
"""
2-dimensional topographies
==========================


"""

from pylab import *
from numpy import *

from matplotlib.tri import Triangulation

from ins.data.layout import Layout

try:
    from matplotlib.tri import UniformTriRefiner
    def trirefine(tri, q, subdiv=3):
        ref = UniformTriRefiner(tri)
        return ref.refine_field(q, subdiv=subdiv)
except ImportError:
    print "UniformTriRefiner not available"
    def trirefine(tri, q, subdiv=3):
        return tri, q

__all__ = ['topo2d']

def topo2d(layname, q, trimask=None, ax=None, chin=None, cbar=False,
        isolines=True):
    """
    """
    lay = Layout(layname)
    if chin is not None:
        q_, trimask = [], []
        for l in lay.labs:
            if l in chin:
                trimask.append(True)
                q_.append(q[chin.index(l)])
            else:
                trimask.append(False)
        q = array(q_)
        trimask = array(trimask)
    tri = lay.compute_triangulation(trimask)
    rtri, rq = trirefine(tri, q)
    if ax is None:
        figure()
        ax = subplot(111, frameon=False)
    cs = ax.tricontourf(rtri, rq, 50, shading='gouraud')
    if cbar:
        ax.figure.colorbar(cs)
    ax.set_xticks([])
    ax.set_yticks([])
    if isolines:
        ax.tricontour(rtri, rq, [0.0], colors=['0.25'], linewidths=[1.])
        ax.tricontour(rtri, rq, 20, colors=['0.5'], linewidths=[0.5])
    return ax


