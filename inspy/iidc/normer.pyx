
import numpy
cimport numpy as np, cython
from cython.parallel import parallel, prange

cdef extern from "math.h" nogil:
    float sqrt(float x)

@cython.boundscheck(False)
@cython.wraparound(False)
def pairwise_norm(np.ndarray[np.float32_t, ndim=3] vs):
    cdef Py_ssize_t vi, vj, si, sj, ci, ne, ns, nc
    cdef float dotij, affij
    cdef np.ndarray[np.float32_t, ndim=2] aff = numpy.zeros((vs.shape[0], vs.shape[0]), numpy.float32)
    ne, ns, nc = vs.shape[0], vs.shape[1], vs.shape[2]
    with nogil:
        for vi in xrange(ne):
            for vj in xrange(ne):
                if vi < vj:
                    continue
                affij = 0.0
                for si in xrange(ns):
                    for sj in xrange(ns):
                        dotij = 0.0
                        for ci in xrange(nc):
                            dotij += vs[vi, si, ci]*vs[vj, sj, ci]
                        affij += dotij*dotij
                aff[vj, vi] = aff[vi, vj] = sqrt(affij)
    return aff

