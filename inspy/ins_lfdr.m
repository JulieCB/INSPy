function [xb, f, cf, fdr, llx] = ins_lfdr(x, nbins, doplot, dc)
%
%    Computes density statistics on `x`, where a "null" hypothesis distribution H0,
%    here a Gaussian distribution, is fit to the center of the data, and this
%    provides an estimation of how likely each bin comes from the Gaussian
%    distribution, a measure called the local false discovery rate.
%
%    Decimating the data by specifying `dc` as an integer greater than 1 
%    can accelerate the density estimation and safer on large datasets 
%    where said estimation takes more time.
%
%    Parameters
%    ----------
%    x : array
%        Data to analyze
%    nbins : int or None, optional
%        Number of points at which to evaluate the density
%    doplot : bool
%        If true, make plot of densities and lfdr
%    dc : integer
%        Decimate data, accelerating estimation
%
%    Returns
%    -------
%    xb : array
%        Points at which the densities are evaluated
%    f : array
%        Density of `x`
%    cf : array
%        Estimated "center density" of H0
%    fdr : array
%        Estimated false discovery rate 
%    llx : array
%        Log lfdr for each datapoint in `x`
%
% defaults: nbins=50, doplot=1, dc=1
%
% mw 11/25/2013 translation from python


% def lfdr(x, nbins=50, doplot=True):

if nargin < 2; nbins  = 50; end
if nargin < 3; doplot =  1; end
if nargin < 4;     dc =  1; end

%% obtain preliminary estimate of density

% k = stats.gaussian_kde(x)
% xb = r_[x.min() : x.max() : 1j*nbins]
% f = k(xb)
% f /= f.sum()

xb = linspace(min(x), max(x), nbins);
[f, ~] = ksdensity(x(1:dc:end), xb);
f = f / sum(f);

%% update query points where data is dense

% dxb = interp(r_[0.0:1.0:1j*nbins], cumsum(f), xb + (xb[1] - xb[0])/2.0)
% xb = unique(r_[xb, dxb])
% xb.sort()                
% f = k(xb)
% f /= f.sum()

dxb = interp1(cumsum(f), xb + (xb(2) - xb(1))/2, linspace(0, 1, nbins));
xb = unique([ xb(:); dxb(:) ]);
xb = sort(xb);
f = ksdensity(x(1:dc:end), xb);
fin = isfinite(f);
f = f(fin);
xb = xb(fin);
f = f / sum(f);



%% estimate center sub-density via optimization of error & density

% def err(par):
%     mu, sigma, alo, ahi = par
%     sl = slice(argmin(abs(xb - alo)), argmin(abs(xb - ahi)))
%     f0 = stats.norm.pdf(xb[sl], loc=mu, scale=sigma)
%     f1 = f[sl]
%     f0 = f0/f0.max()*f1.max()
%     return sum((f0 - f1)**2)/sum(f1**2) - f1.sum()

function e = err(par)
    mu = par(1); sigma = par(2); alo = par(3); ahi = par(4);
    [~, ilo] = min(abs(xb - alo));
    [~, ihi] = min(abs(xb - ahi));
    f0 = normpdf(xb(ilo:ihi), mu, sigma);
    f1 = f(ilo:ihi);
    f0 = f0/max(f0) * max(f1);
    e = sum((f0 - f1).^2)/sum(f1.^2) - sum(f1);
end

% mu0, sig0 = x.mean(), x.std()
% mu, sig, _, _ = optimize.fmin(err, (mu0, sig0, mu0-sig0, mu0+sig0), disp=0)
% cf = stats.norm.pdf(xb, loc=mu, scale=sig)
% cf = cf/cf.max()*f.max()

mu0 = mean(x); sig0 = std(x);
opar = fminsearch(@err, double([mu0, sig0, mu0 - sig0, mu0 + sig0]));
mu = opar(1); sig = opar(2);
cf = normpdf(xb, mu, sig);
cf = cf/max(cf) * max(f);
        
%% compute lfdr & transform data to log-lfdr

% fdr = clip(cf/f, 0.0, 1.0)
% llx = interp(x, xb, log(fdr))

fdr = cf ./ f;
fdr(fdr<0.0) = 0.0;
fdr(fdr>1.0) = 1.0;
llx = interp1(xb, log(fdr), x);

% if doplot:
%     import pylab as pl
%     pl.semilogy(xb, f, 'k')
%     pl.semilogy(xb, cf, 'k--')
%     pl.semilogy(xb, fdr, 'k.')
%     pl.grid(True)
%     pl.ylim([1e-5, 1.0])
%     plevels=[0.2, 0.1, 0.05, 0.01, 0.001]
%     pl.yticks(plevels, map(str, plevels))

if doplot
     figure
     semilogy(xb, f, 'k');
     hold on
     semilogy(xb, cf, 'k--');
     semilogy(xb, fdr, 'k.');
     grid on
     ylim([1e-5, 1.0]);
     set(gca, 'ytick', sort([0.2, 0.1, 0.05, 0.01, 0.001]));
end


end
