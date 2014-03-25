function nfpeaks = ins_detect(t, y, flo, fhi, q, minel, maxel, rf)
%     Perform detection on continuous signal `y` as function of time `t`. 
% 
%     This function filters the signal, estimates the LFDR (see `lfdr`), and
%     identifies where `y` remains "unlikely" (LFDR < `q`), for a period of
%     time longer than `minel` and shorter than `maxel`, with a refractory
%     period of `rf`.
% 
% 
%     Parameters
%     ----------
% 
%     t : array
%         Sequence of times corresponding to the samples in `y`, 
%             expected to be equally spaced.
%     y : array
%         Continuously sampled signal on which to perform detection
%     flo : float
%         Low frequency of the pass band
%     fhi : float
%         High frequency of the pass band
%     q : float
%         Threshold on LFDR
%     minel : float
%         Minimum event length, seconds
%     maxel : float
%         Maximum event length, seconds
%     rf : float
%         Refractory period of events, seconds
% 
%     Returns
%     -------
% 
%     ret : dict or array
%         An array of event peak times, or if `info` is `True`, a
%             dictionary of the local variables.
% 
% defaults: flo=15., fhi=40., q=1e-5, minel=0.05, maxel=0.5, rf=0.5, info=False
% 
% mw 11/25/2013 translation from Python
% 

if nargin < 3, flo=15.; end
if nargin < 4, fhi=40.; end
if nargin < 5, q=1e-5; end
if nargin < 6, minel=0.05; end
if nargin < 7, maxel=0.5; end
if nargin < 8, rf=0.5; end


% fs = 1/(t[1] - t[0])

fs = 1/(t(2) - t(1));

%% filter the temporally diff

% b, a = signal.butter(3, [2*flo/fs, 2*fhi/fs], 'pass')
% fy = signal.filtfilt(b, a, diff(y))

[b, a] = butter(3, [2*flo/fs, 2*fhi/fs], 'bandpass');
fy = filtfilt(b, a, diff(y));

%% analyze distribution

% xb, f, cf, fdr, llx = lfdr(fy, doplot=False)

[xb, f, cf, fdr, llx] = ins_lfdr(fy, 50, 0);

%% fdr tarnsform hilbert 

% hy = abs(signal.hilbert(fy))
% llh = interp(hy, xb, log(fdr))

hy = abs(hilbert(fy));
llh = interp1(xb, log(fdr), hy);

%% generate events    

% ev = c_[isfinite(llh), llh < log(q)].all(axis=1)
% ev[~isfinite(llh)] = True
% e0, = argwhere(c_[~ev[:-1], ev[1:]].all(axis=1)).T
% e1, = argwhere(c_[ev[:-1], ~ev[1:]].all(axis=1)).T

ev = isfinite(llh) & llh < log(q);
ev(~isfinite(llh)) = 1;
e0 = find(~ev(1:end-1) & ev(2:end));
e1 = find(ev(1:end-1) & ~ev(2:end));

% if len(e0) == 0 or len(e1) == 0:
%     return locals() if info else []

if isempty(e0) || isempty(e1)
    return
end

%% boundaries

% if e1[0] < e0[0]:
%     e0 = r_[0, e0]
% if e0[-1] > e1[-1]:
%     e1 = r_[e1, llh.shape[0]-1]

if e1(  1) < e0(  1), e0 = [    0;       e0(:)]; end
if e0(end) > e1(end), e1 = [e1(:); length(llh)]; end

%% compute event length and long-enough mask

% el = diff(c_[e0, e1])[:, 0]/fs
% le = c_[el > minel, el < maxel].all(axis=1)

el = diff([e0(:) e1(:)]')/fs;
le = el > minel & el < maxel;

%% pull out remaining peaks & align them

% peaks = []
% for i, (e0i, e1i) in enumerate(zip(e0[le], e1[le])):
%     hi = hy[e0i:e1i]
%     peaks.append((hi.max(), t[e0i] + argmax(hi)/fs))

peaks = [];
e0_le = e0(le);
e1_le = e1(le);
for i=1:length(e0_le)
    hi = hy(e0_le(i) : e1_le(i));
    [peaks(i).hmax, offset] = max(hi);
    peaks(i).time = t(e0_le(i)) + offset/fs;
end

% if len(peaks) == 0:
%     return locals() if info else []

if isempty(peaks)
    return
end

%% mask refractory period
% ph, pt = array(peaks).T
% nonmask = []
% for i, (phi, pti) in enumerate(peaks):
%     mask = c_[pt > pti - rf, pt < pti + rf].all(axis=1)
%     if (phi >= ph[mask]).all():
%         nonmask.append(i)
% 
% nfpeak = array(peaks)[array(nonmask), 1]
 
ph = [peaks(:).hmax]; pt = [peaks(:).time];
nonmask = [];
for i=1:length(peaks)
    phi = peaks(i).hmax;
    pti = peaks(i).time;
    mask = pt > pti - rf & pt < pti + rf;
    if all(phi >= ph(mask))
        nonmask = [nonmask(:); i];
    end
end

nfpeaks = [peaks(nonmask).time];
