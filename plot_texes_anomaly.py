#! /bin/python3
import matplotlib
matplotlib.use('Agg')
from pylab import *
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable

hdul = fits.open('Gemini_TEXES_GRS_2017_product.fits')

hdu = hdul[0]
lat1 = hdu.header['LATS']
lat2 = hdu.header['LATN']
dlat = hdu.header['LATSTEP']
lat = arange(lat1, lat2 + dlat, dlat)

lon1 = hdu.header['LONW']
lon2 = hdu.header['LONE']
dlon = hdu.header['LONSTEP']
lon = arange(lon1, lon2 + dlon, dlon)

hdu = hdul[2]
pres = hdu.data/100.  # pa to mbar

fig, axs = subplots(2, 2, figsize = (12, 10), sharey = True,
  gridspec_kw={'width_ratios':[3,1]})
subplots_adjust(hspace = 0.08, wspace = 0.08)
X, Y = meshgrid(lat, pres)

# temperature
var = hdul[0].data
ax = axs[0,0]
avg = mean(var[:10,:,:], axis = 0)
vara = var[37,:,:] - avg

h = ax.contour(X, Y, vara.T, linspace(-11, -1, 6), colors = 'b')
clabel(h, fontsize = 11, inline = 1, fmt = '%.0f')

h = ax.contour(X, Y, vara.T, linspace(1, 11, 6), colors = 'r')
clabel(h, fontsize = 11, inline = 1, fmt = '%.0f')

ax.plot([-24, -24], [max(pres), min(pres)], 'C1-')
ax.plot([-16, -16], [max(pres), min(pres)], 'C1-')
ax.plot([-20, -20], [max(pres), min(pres)], 'C1--')
ax.set_xlim([-26, -14])
ax.set_ylabel('Pressure (mbar)')
ax.xaxis.tick_top()

ax = axs[0,1]
ax.plot(mean(avg, axis = 0), pres)
ax.set_xlabel('Avg. temperature (K)')
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')

var = hdul[1].data*1.E6  # ppmv
ax = axs[1,0]
avg = mean(var[:10,:,:], axis = 0)
vara = var[37,:,:] - avg
h = ax.contour(X, Y, vara.T, linspace(-41, -1, 5), colors = 'b')
clabel(h, fontsize = 11, inline = 1, fmt = '%d')

h = ax.contour(X, Y, vara.T, linspace(1, 11, 11), colors = 'r')
clabel(h, fontsize = 11, inline = 1, fmt = '%d')

ax.plot([-24, -24], [max(pres), min(pres)], 'C1-')
ax.plot([-16, -16], [max(pres), min(pres)], 'C1-')
ax.plot([-20, -20], [max(pres), min(pres)], 'C1--')
ax.set_xlim([-26, -14])
ax.set_ylim([max(pres), 100.])
ax.set_yscale('log')
ax.set_xlabel('PC latitude')
ax.set_ylabel('Pressure (mbar)')

ax = axs[1,1]
ax.plot(mean(avg, axis = 0), pres)
ax.set_xscale('log')
ax.set_xlim([1.E-4, 1.E3])
ax.set_xlabel('Avg. ammonia (ppmv)')

savefig('TEXES_temp_ammonia_anomaly.png', bbox_inches = 'tight')
