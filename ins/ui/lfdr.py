
"""
Visualizing local false discovery rates. 

"""

from pylab import *
from numpy import *
from scipy import stats
from ins.util import *
from ins.data import ades, mrk
from ins.viz import tsplot
from ins import algo

reload(algo), reload(tsplot)

seeg = ades.Ades(path=ff.marked.ades[0])

y = seeg.data[4].values

marks = mrk.Mrk(fd.marked.files.seeg006.mrk[0])
marks.timei = int64(marks.times*seeg.fs)


figure(1), clf()

reload(algo)
ax221 = subplot(221)
cla()
ylo, yhi, rest = algo.lfdr(y, details=True)
b, fz, f0z, f1z, fdr = [rest[key] for key in ['bs', 'fz', 'f0z', 'f1z', 'fdr']]
plot(b, fz, 'b')
plot(b, f0z, 'r')
plot(b, f1z, 'g')
gca().set_yscale('log')
autoscale()
grid(True)

# look in an area of 100 ms forward and backward for biggest value
hwin = int(0.100*seeg.fs)
epochs = []
for mi in marks.timei:
    epochs.append(y[mi-hwin : mi+hwin])
epochs = array(epochs)
maxs = abs(epochs).max(axis=1)
imax = argmax(abs(epochs), axis=1)
subplot(222)
plot(epochs.T, 'k', alpha=0.4)
peaks = epochs.flat[r_[:len(imax)]*epochs.shape[1] + imax]
plot(imax, peaks, 'ro', alpha=0.4)
subplot(224)
apochs = epochs.copy()
for i, (ap, peak_idx) in enumerate(zip(apochs, imax)):
    apochs[i] = roll(ap, hwin - peak_idx)
plot(apochs.T, 'k', alpha=0.4);

subplot(223, sharex=ax221), cla()
_c, _b, _ = hist(peaks, 100)
vlines(algo.thresh2side(b, f1z>f0z), 0, max(_c), color='r', linewidth=2)
vlines(algo.thresh2side(b, fdr<0.05), 0, max(_c), color='y', linewidth=2)
grid(True)



figure(2), clf()

sig = seeg.data[4].values
for i, y in enumerate([sig, algo.whiten(sig)]):
    ax = subplot2grid((2, 6), (i, 4))
    ylo, yhi, rest = algo.lfdr(y, details=True, refine=False, nbins=500)
    b, fz, f0z, f1z, fdr = [rest[key] for key in ['bs', 'fz', 'f0z', 'f1z', 'fdr']]
    plot(fz, b, 'b')
    plot(f0z, b, 'r')
    plot(f1z, b, 'g')
    gca().set_xscale('log')
    autoscale()
    grid(True)
    yticks([])
    subplot2grid((2, 6), (i, 5), sharey=ax)
    _c, _b = histogram(peaks, 100)
    plot(_c, _b[:-1], 'ko-')
    hlines(algo.thresh2side(b, f1z>f0z), 0, max(_c), color='r', linewidth=2)
    grid(True)
    yticks([])
    subplot2grid((2, 6), (i, 0), colspan=4, sharey=ax)
    ts = r_[:len(y)]/seeg.fs
    plot(ts, y, 'k')
    vlines(marks.times, y.min(), y.max(), color='r')
    hlines(algo.thresh2side(b, f1z>f0z), 0, ts.max(), color='r', linewidth=3, alpha=0.5)
    hlines(algo.thresh2side(b, fdr<0.05), 0, ts.max(), color='y', linewidth=3, alpha=0.5)
    xlim([0, 10])


