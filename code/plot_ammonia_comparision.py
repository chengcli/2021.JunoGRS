#! /usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.io import fits
from pylab import *

case = 'PJ07BOOSTV3'
hdul = fits.open('../%s/profile_PJ07_2.7x20.x168.fits' % case)
nh3_grs = hdul[0].data/2.7*350. # g/kg -> ppm
lat = hdul[3].data
pres = hdul[4].data # bar
pavg = mean(pres, axis = 0)

case = 'PJ1345689V3'
hdul = fits.open('../%s/profile_PJ1345689_2.7x20.x168.fits' % case)
nh3_jets = hdul[0].data/2.7*350. # g/kg -> ppm
lat = hdul[3].data
pres = hdul[4].data # bar
pavg = mean(pres, axis = 0)

X, Y = meshgrid(lat, pavg)

fig, axs = subplots(2, 1, figsize = (8, 12), sharex = True, sharey = True)
subplots_adjust(hspace = 0.08, wspace = 0.08)
ax = axs[0]
h = ax.contourf(X, Y, nh3_jets.T, linspace(40, 360, 17),
  cmap = 'inferno', extend = 'min')

# add colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size = "5%", pad = 0.2, aspect = 20)
#c = colorbar(h, cax = cax, label = 'NH$_3$ mixing ratio (ppmv)')
c = colorbar(h, cax = cax)
c.ax.invert_yaxis()

ax.set_yscale('log')
ax.set_ylim([100., 0.3])
#ax.set_xlabel('Planetocentric latitude', fontsize = 15)
ax.set_ylabel('Pressure (bar)', fontsize = 15)
ax.text(-25.5, 80., '(a)', fontsize = 24)


ax = axs[1]
#nh3_diff = nh3_grs - nh3_jets
#h = ax.contourf(X, Y, nh3_diff.T, linspace(-100, 100, 17),
#  cmap = 'coolwarm', extend = 'both')
h = ax.contourf(X, Y, nh3_grs.T, linspace(40, 360, 17),
  cmap = 'inferno', extend = 'min')

# add colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size = "5%", pad = 0.2, aspect = 20)
#c = colorbar(h, cax = cax, label = 'GRS NH$_3$ mixing ratio anomaly (ppmv)')
c = colorbar(h, cax = cax)
#c.ax.invert_yaxis()

#ax.contour(X, Y, nh3_diff.T, [0.], linewidths = 2, colors = '0.7')
ax.plot([-16., -16.], [pavg.max(), pavg.min()], color = 'w')
ax.plot([-24., -24.], [pavg.max(), pavg.min()], color = 'w')
#colorbar(h, label = 'NH$_3$ mixing ratio anomaly (ppmv)')
#ax.set_yscale('log')
#ax.set_ylim([100., 0.3])
ax.set_xlabel('Planetocentric latitude', fontsize = 15)
ax.set_ylabel('Pressure (bar)', fontsize = 15)
ax.text(-25.5, 80., '(b)', fontsize = 24)

savefig('../figs/ammonia_comparison_F3.png', bbox_inches = 'tight', dpi = 400)
savefig('../figs/ammonia_comparison_F3.pdf', bbox_inches = 'tight')
