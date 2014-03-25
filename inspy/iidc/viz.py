
"""
Visualizations

"""

import numpy as np
#import pylab as pl

from scipy import stats

from . import cluster, preprocess, util, stat

def browse_signals(fs, data, hp=2.0, newfigure=True, ax=None):
    if ax is None:
        fig = pl.figure()
        ax = fig.add_subplot(111)
    data = preprocess.hp(fs, data, hp=hp)
    t = pl.r_[:data.shape[1]]*2.0/fs
    y = (data.T - data.mean(axis=1))/data.ptp(axis=1) + pl.r_[:len(data)]
    return ax.plot(t, y, 'k')

def embedded_scatter(xi, labels):
    labels = labels.astype(int)
    colors = 'rgbcykv'
    for i, lab in enumerate(np.unique(labels)):
        x, y = xi[:, labels==lab]
        pl.plot(x, y, colors[i]+'x')
    #pl.legend(['class %d' % (lab, ) for lab in np.unique(labels)])

    
def class_resume(fs, data, events, labels, win_pre=-0.1, win_post=0.2, 
                 ch_labels=None, alpha=0.05):
    _, ep = util.extract_windows(fs, data, events)
    pl.figure()
    if ch_labels is None:
        ch_labels = map(str, np.r_[:data.shape[0]])
    outs = []
    for i, lab in enumerate(np.unique(labels)):
        epi = ep[labels==lab]
        _, P = stats.ttest_1samp(epi, 0.0, axis=0)
        PT, _ = stat.fdr(P, alpha=alpha)
        sig = stat.threshold_duration(P.T < 0.05, 0.01 * fs)
        ax = pl.subplot(1, len(np.unique(labels)), i+1)
        pl.imshow(-np.log(P.T*sig)*np.sign(epi.mean(axis=0).T), 
                  cmap='RdBu', 
                  interpolation='none', 
                  aspect='auto', 
                  extent=[win_pre*1e3, win_post*1e3, ep.shape[2], 0])
        ax.set_xticks(np.r_[win_pre:win_post:0.05]*1e3)
        ax.set_yticks(np.r_[:data.shape[0]]-0.5)
        ax.set_yticklabels(ch_labels)
        pl.title('Class %d, %d events' % (lab, len(epi)))
        pl.grid(True)
        pl.xlabel('time (ms)')
        outs.append({'epi': epi, 'P':P, 'PT':PT, 'sig':sig})
    pl.tight_layout()
    return outs

def class_dynamics(fs, data, ev, labels):
    pl.figure()
    ncl = len(np.unique(labels))
    browse_signals(fs, data, ax=pl.subplot(ncl+1, 1, 1))
    pl.xticks([])
    labels = labels.astype(int)
    for i, lab in enumerate(np.unique(labels)):
        pl.subplot(ncl+1, 1, i+2)
        pl.hist(ev[labels==lab], 100)
        pl.xticks([])

def plot2d(data, **kwds):
    """
    Embed data in a 2d manifold and make a scatter plot.

    """

    if 'method' not in kwds:
        fig  = pl.figure(figsize=(12, 5))
        methods = cluster.embed_manifold.methods
        ys = []
        for i, meth in enumerate(methods):
            pl.subplot(1, len(methods), i + 1)
            ys.append(plot2d(data, method=meth))
            pl.title(meth)
        return ys

    Y, mnf = cluster.embed_manifold(data, **kwds)
    x, y = Y.T
    pl.plot(x, y, 'kx')
    return Y




if False:
    elec_ix, elec_nm = util.seeg_yticks(fte)

    figure(figsize=(5*ncl, 10))
    for i in range(ncl):
        ep = epochs[dbscan.labels_==i]
        _, P = stats.ttest_1samp(ep, 0.0, axis=0)
        PT, _ = stat.fdr(P, alpha=0.01)
        sig = stat.threshold_duration(P.T < PT, 0.01 * fs)
        ax = subplot(1, ncl, i+1)
        imshow(log(sig*P.T)*sign(ep.mean(axis=0).T), 
               cmap='RdBu', 
               interpolation='none', 
               aspect='auto', 
               extent=[win_pre*1e3, win_post*1e3, epochs.shape[2], 0])
        ax.set_yticks(elec_ix)
        ax.set_yticklabels(elec_nm if i == 0 else [])
        ax.set_xticks(r_[win_pre:win_post:0.05]*1e3)
        title('Class %d, %d events' % (i+1, len(ep)))
        grid(True)
        xlabel('time (ms)')

    figure(figsize=(5*ncl, 15))
    t = linspace(win_pre, win_post, epochs.shape[1])
    selection = ['A_01', 'B_01', 'GPH_01', 'GPH_03', 'TB_03', 'TP_02', 'OF_05']
    chanlabel = [c.name for c in fte.chan]
    print chanlabel
    selection_ix = array([chanlabel.index(s) for s in selection])
    for i in range(ncl):
        ep = rollaxis(epochs[dbscan.labels_==i][:, :, selection_ix], 2)
        for j, (chanlab, epj) in enumerate(zip(selection, ep)):
            ax = subplot(len(ep), ncl, i + j*ncl + 1)
            if True:
                plot(t, epj.T, 'k', alpha=0.1)
                mu, sd = epj.mean(axis=0), epj.std(axis=0)
                [plot(t, y, c) for y, c in [(mu, 'r'), (mu+sd, 'r--'), (mu-sd, 'r--')]]
            else:
                dplot.dplot(tile(t, (epj.shape[0], 1)), epj, 40)
            yticks([0.0], [chanlab])
            xticks([0.0])
            if j == 0:
                title('Class %d, %d events' % (i + 1, ep.shape[1]))
            if i != 0:
                yticks(yticks()[0], [])
            if j != len(ep)-1:
                xticks(xticks()[0], [])
            grid(True)
    tight_layout(0)
