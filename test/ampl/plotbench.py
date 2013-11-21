from pylab import *
from numpy import *

file = 'benchdetect.dat'
columns = (0, 2, 3, 4, 5, 6, 7)

results = loadtxt(file, usecols=columns)
results = rollaxis(results.T.reshape((-1, 2, 2, 20, 20)), 0, 3)

dsnames = ['C2-3', 'TB3-4']
extent = [0, 100, 100, 0]

for i, whiten in enumerate(results):
    for j, dataset in enumerate(whiten):
        wh, lo, hi, auc, bestp, oflo, ofhi = dataset
        subplot(2, 2, i*2 + j + 1)
        print i, j
        imshow(auc, extent=extent, vmin=0, vmax=1)
        colorbar()
        if i==0:
            title(dsnames[j])
        if not i==results.shape[0]-1:
            xticks([])
        if j==0:
            pass
        else:
            yticks([])

