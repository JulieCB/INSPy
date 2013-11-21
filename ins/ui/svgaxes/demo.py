from pylab import *

from svgaxes import read_axes

sw, sh, axs = read_axes('axes.svg')

figure(figsize=(8, 8*sh/sw))

for k, a in axs.iteritems():
    axes(a)
    plot(randn(30, 30) + r_[:30], 'k')
    title(k)

show()
