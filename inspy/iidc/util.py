

"""
Utilities, not organized.

"""

import re
import os
import glob
import inspect
import os.path
import contextlib
import multiprocessing

import numpy as np


def which(program):
    """
    Find executable on path.

    Shamelessly borrowed from Jay, at 
    
        http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python

    """

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def find_file(path, patt):
    """
    Find files matching `patt` on path `path`, e.g.

    >>> find_hdr('/data/local/simultane/*/*gou_mi/', 'hdr.mat')
    [u'/data/local/.../hdr.mat',
     ...]

    The returned list is sorted.

    """

    files = []
    path_ = unicode(glob.glob(path)[0])
    for r, _, fs in os.walk(path_):
        for f in fs:
            if re.match('.*' + patt + '.*', os.path.join(r, f)):
                files.append(os.path.join(r, f))

    return sorted(files)


def bti_import_fix_channel_names(data, pdfname):
    """
    Fix channel names on imported 4D data. We keep most channel's original
    names except for EEG where we keep the chan_label. 

    """

    info_ = _read_bti_header(pdfname, pdfname[:-6] + 'config')
    for i, ch in enumerate(info_['chs']):
        rename = ch['name']
        if rename[0] == 'E':
            rename = ch['chan_label']
        data.ch_names[i] = data.info['chs'][i]['ch_name'] = rename


def seeg_ch_name_split(nm):
    """
    Split an sEEG channel name into its electrode name and index

    >>> seeg_ch_name_split('GPH10')
    ('GPH', 10)

    """

    elec, idx = re.match(r'([A-Za-z]+)(\d+)', nm).groups()
    return elec, int(idx)


def apply_picks_to_proj(raw, picks):
    """
    Work around.

    """

    projs = raw.info['projs']
    projs_ = []
    for proj in projs:
        proj_data = proj['data']
        proj_data['col_names'] = list(np.array(proj_data['col_names'])[picks])
        proj_data['data']      = proj_data['data'][:, picks]
        proj_data['ncol']      = len(proj_data['col_names'])
        proj['data'] = proj_data
        projs_.append(proj)
    raw.info['projs'] = projs_


def raw_sort_channels(raw):
    """
    Sort channels in raw

    """

    # use a leading zero so A02 < A10, not A2 > A10
    def key(name):
        try:
            elec, idx = seeg_ch_name_split(name)
            keyed = '%s%03d' % (elec, idx)
        except:
            keyed = name
        return keyed

    ch_names = sorted(raw.ch_names, key=key)
    chix = np.array(map(raw.ch_names.index, ch_names))

    raw.info['chs'] = [raw.info['chs'][i] for i in chix]
    raw.info['ch_names'] = ch_names
    raw._data = raw._data[chix]

    apply_picks_to_proj(raw, chix)

    return raw


def find_bip_idx(seeg):
    """
    Generate the set of indices for sEEG electrodes

    """

    for i in xrange(len(seeg.ch_names)-1):
        (e1, i1), (e2, i2) = map(seeg_ch_name_split, seeg.ch_names[i:i+2])
        if e1==e2:
            yield i, i+1, e1, i1, i2


def create_bipolar_montage(seeg, copy=True):
    """
    Create a bipolar dataset from `seeg`, assuming channels have already
    been sorted. 

    """

    if copy:
        seeg = seeg.copy()

    bip_data = []
    bip_chan = []

    #   16, 17, 'OF',   2,   3  results in 'OF2-3', _data[16] = _data[16] - _data[17]
    for i1, i2, elec, ei1, ei2 in find_bip_idx(seeg):
        ch_info = seeg.info['chs'][i1].copy()
        ch_info['ch_name'] = '%s%d-%s%d' % (elec, ei1, elec, ei2)
        bip_chan.append(ch_info)
        bip_data.append(seeg._data[i1] - seeg._data[i2])
        
    seeg.info['chs'] = bip_chan
    seeg._data = np.array(bip_data)
    seeg.info['ch_names'] = [c['ch_name'] for c in seeg.info['chs']]
    seeg.info['nchan'] = seeg._data.shape[0] 

    return seeg


def _parse_ch_list(chs):
    if type(chs) in (str, unicode):
        chs = re.split(',| ', chs)    
    return chs


def pick_bip_chan(ch_names, include, exclude=[], delim='-'):
    """
    Find indices of bipolar channels containing and excluding certain monopolar
    or bipolar channel names. 
    
    - ch_names : list
        Channel names as in `raw.ch_names`
    - include : str, unicode, list
        Channels to pick, may be a str/unicode with names separated by commas or spaces. 
    - exclude : str, unicode, list
        Channels to exclude, may be string as with include
        
    """
    
    include = _parse_ch_list(include)
    exclude = _parse_ch_list(exclude)
    
    idx = []
    for i, nm in enumerate(ch_names):
        ch, ref = nm.split(delim)
        inc = ch in include or ref in include
        exc = ch in exclude or ref in exclude
        if inc and not exc:
            idx.append(i)
    
    return array(idx)


def seeg_yticks(run):
    """
    Format channel labels nicely for plotting.

    """

    ix = [0]
    nm = [run.chan[0].name.split('_')[0]]
    for i in range(1, len(run.chan)-1):
        (nm1, _), (nm2, _) = [c.name.split('_') for c in run.chan[i:i+2]]
        if nm1 != nm2:
            ix.append(i+1)
            nm.append(nm2)
    return ix, nm


def extract_windows(fs, ys, peaks, pre=-0.1, post=0.2):
    """
    Extract windows around peaks in ys.

    """

    peaks_kept = []
    epochs = []
    for peaki in peaks:
        ilo = int((peaki + pre)*fs)
        ihi = int((peaki + post)*fs)
        if ilo < 0 or ihi > ys.shape[1]:
            continue
        epoch = ys[:, ilo:ihi].T
        epochs.append(epoch)
        peaks_kept.append(peaki)
    return np.array(peaks_kept), np.array(epochs)

def combine_events(events, lim=0.02, max=5000):
    """
    Combine events not farther apart than `lim`.

    """

    events = np.concatenate([e for e in events if isinstance(e, np.ndarray)])
    events.sort()
    events = events[np.diff(events)>lim]
    return events[::int(len(events)/max + 1)]

@contextlib.contextmanager
def mpool(n_jobs=1, **kwds):
    """
    Thin wrapper around the multiprocessing Pool class to make it a context 
    manager that automatically closes itself on exiting the context. 

    The number of jobs defaults to 1 so user must explicitly provide a number
    or None.

    """

    if n_jobs > 1:
        pool = multiprocessing.Pool(n_jobs, **kwds)
        yield pool
        pool.close()
    else:
        import itertools
        class FakePool(object):
            def imap(self, f, it, chunksize=10):
                return (f(i) for i in it)
            def map(self, f, it, chunksize=10):
                return [f(i) for i in it]
        yield FakePool()
        

def exists(names, namespace=None):
    """
    Tests for a name or space-seperated names in the scope of the caller:

    >>> exists('a')
    False
    >>> a = 3
    >>> exists('a')
    True
    >>> exists('a b')
    False
    >>> b = 4
    >>> exists('a b')
    True

    We could imagine a more general checkpointing mechanism that takes a file
    name and set of variable names. 

    - if variables in caller's scope return False
    - if vars not in scope, but found file, load file, put vars into scope and return False
    - if vars and file not found, print msg, run computation, save to file?
    
    but afterwards, we'd need to save to file..

    """

    names = names.split(' ') if ' ' in names else [names]
    if namespace is None:
        namespace = inspect.currentframe().f_back.f_globals
    existing = namespace.keys()
    for name in names:
        if name not in existing:
            return False
    return True


