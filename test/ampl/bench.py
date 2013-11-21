import multiprocessing as mp, itertools as it, time, numpy as np
from scipy import io

import ins.detect.ampl as detect

if __name__ == "__main__":

    pool = mp.Pool(4)
     
    ngrid = 20
    pars = list(it.product((False, True),
                           ['C2-C3', 'TB3-TB4'], 
                           np.r_[0.1:100:1j*ngrid], 
                           np.r_[0.1:100:1j*ngrid]))
    print len(pars), 'est time', len(pars)*2.5/3600.0, 'h'
    tic = time.time()
    results = pool.map(detect.job, pars)
    print (time.time() - tic)/len(pars), 's per test'

    fd = open('benchdetect.dat', 'w')
    for (wh, ds, lo, hi), (auc, bestp, oflo, ofhi) in zip(pars, results):
        wh = 1 if wh else 0
        print >> fd, wh, ds, lo, hi, auc, bestp, oflo, ofhi
        fd.flush()

    fd.close()

