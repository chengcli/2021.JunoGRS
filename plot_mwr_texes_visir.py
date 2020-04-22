#! /bin/python3
import matplotlib
matplotlib.use('Agg')
from pylab import *
from astropy.io import fits
from scipy.io import readsav
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MaxNLocator

def Rad2Temp(imgdata, um):
  imgdata *= 1.e-7  # ergs/s/cm2/sr/cm-1 -> W/cm2/sr/cm-1
  h = 6.626e-34
  c = 2.9979e8
  k = 1.3806e-23
  v = (1e4/um)*100. # m-1
  spec = imgdata*100.*100./100. # W m-2 sr-1 (m-1)-1
  c1 = 2*h*c*c
  c2 = h*c/k
  a = c1*v*v*v/spec
  return c2*v/log(a+1)

# visir cylindrical projection
lon1 = -linspace(-360, 0, 720)  # to west longitude
lat1 = linspace(-90, 90, 360)
X1, Y1 = meshgrid(lon1, lat1)

# Q1, 17.6 um
hdulist = fits.open('VLT-VISIR/data/2017jul13/cmap_cal_wvisir_Q1_2017-07-13T23:25:15.4253_Jupiter_filt.fits')
tb_q1 = hdulist[0].data[:]
tb_q1 = Rad2Temp(tb_q1, 17.6)

# NEII_2, 13 um
hdulist = fits.open('VLT-VISIR/data/2017jul13/cmap_cal_wvisir_NEII_2_2017-07-13T23:03:39.7348_Jupiter_filt.fits')
tb_neii2 = hdulist[0].data[:]
tb_neii2 = Rad2Temp(tb_neii2, 13.)

# SIV_2, 10.7 um
hdulist = fits.open('VLT-VISIR/data/2017jul13/cmap_cal_wvisir_SIV_2_2017-07-13T22:59:13.1116_Jupiter_filt.fits')
tb_siv2 = hdulist[0].data[:]
tb_siv2 = Rad2Temp(tb_siv2, 10.7)

# TEXES inverted data
data = readsav('GeminiTEXES2017/GreatRedSpot/pcentric/temp300mbar.sav')
T300 = data['img_new']

data = readsav('GeminiTEXES2017/GreatRedSpot/pcentric/temp600mbar.sav')
T600 = data['img_new']

data = readsav('GeminiTEXES2017/GreatRedSpot/pcentric/ammoniavmr.sav')
lon2 = data['ulon'] + 40. # shift in longitude because TEXES data were taken three months before MWR PJ7
lat2 = data['ulat']
nh3 = data['img_new']
X2, Y2 = meshgrid(lon2, lat2)

# visir inverted data
#lon2 = linspace(-90, -30, 121)
#lat2 = linspace(-30, 0, 61)
#hdulist = fits.open('VLT-VISIR/GreatRedSpot2017/2017-Jul-13_nh3.fits')
#nh3 = hdulist[0].data[:]
#hdulist = fits.open('VLT-VISIR/GreatRedSpot2017/2017-Jul-13_T_290.fits')
#T290 = hdulist[0].data[:]
#hdulist = fits.open('VLT-VISIR/GreatRedSpot2017/2017-Jul-13_T_500.fits')
#T500 = hdulist[0].data[:]

# mwr data
tb_ch6 = genfromtxt('MWR-PJ7-maps/Tb_nadir_freq5_main_beam_on_planet.txt')
lat_ch6 = genfromtxt('MWR-PJ7-maps/latitude_planetocentric_freq5_main_beam_on_planet.txt')
lon_ch6 = genfromtxt('MWR-PJ7-maps/longitude_planetocentric_freq5_main_beam_on_planet.txt')

ix = where((lat_ch6 < 5) & (lat_ch6 > -30) & (lon_ch6 < -30) & (lon_ch6 > -90.))
clevels = linspace(135, 151, 17)

fig, axs = subplots(3, 2, figsize = (12, 8), sharex = True, sharey = True)
subplots_adjust(hspace = 0.08, wspace = 0.08)

ax = axs[0,1]
h1 = ax.pcolormesh(X1, Y1, tb_q1, cmap = 'gist_heat', vmin = 110, vmax = 120)
ax.text(85, 1, 'VISIR, 17.6 $\mu$m', color = 'w', fontsize = 12)
ax.plot([51.6, 66.4], [4.8, -30], 'w--')
ax.plot([49.1, 57.8], [4.8, -30], 'w--')
# add colorbar
divider = make_axes_locatable(axs[0,1])
cax = divider.append_axes("right", size = "10%", pad = 0., aspect = 20)
colorbar(h1, cax = cax)

ax = axs[1,1]
h2 = ax.pcolormesh(X1, Y1, tb_neii2, cmap = 'gist_heat', vmin = 122, vmax = 133)
ax.text(85, 1, 'VISIR, 13 $\mu$m', color = 'w', fontsize = 12)
ax.plot([51.6, 66.4], [4.8, -30], 'w--')
ax.plot([49.1, 57.8], [4.8, -30], 'w--')
# add colorbar
divider = make_axes_locatable(axs[1,1])
cax = divider.append_axes("right", size = "10%", pad = 0., aspect = 20)
colorbar(h2, cax = cax)


ax = axs[2,1]
h3 = ax.pcolormesh(X1, Y1, tb_siv2, cmap = 'gist_heat', vmin = 124, vmax = 134)
ax.text(85, 1, 'VISIR, 10.8 $\mu$m', color = 'w', fontsize = 12)
ax.plot([51.6, 66.4], [4.8, -30], 'w--')
ax.plot([49.1, 57.8], [4.8, -30], 'w--')
# add colorbar
divider = make_axes_locatable(axs[2,1])
cax = divider.append_axes("right", size = "10%", pad = 0., aspect = 20)
colorbar(h3, cax = cax)
ax.set_xlabel('West longitude')


ax = axs[0,0]
#ax.pcolormesh(X2, Y2, T300, cmap = 'binary_r')
#ax.tricontourf(lon_ch6[ix], lat_ch6[ix], tb_ch6[ix], levels = clevels, cmap = 'RdBu_r')
h4 = ax.contourf(X2, Y2, T300, cmap = 'binary_r', levels = linspace(120, 129, 19), alpha = 0.8)
ax.scatter(-lon_ch6[ix], lat_ch6[ix], s = 6, marker = 's', c = tb_ch6[ix], cmap = 'RdBu_r')
# add colorbar
cax = inset_axes(ax, width = "6%", height = "60%", loc = 3,
  bbox_to_anchor = (0.8, 0.08, 1, 1), bbox_transform=ax.transAxes)
colorbar(h4, cax = cax, ticks = [120, 123, 126, 129])
cax.set_xlabel('K')
cax.xaxis.set_label_position('top')
ax.text(85, 1, 'Temperature, 300 mbar', fontsize = 12)
ax.set_ylabel('PC latitude')

ax = axs[1,0]
#ax.pcolormesh(X2, Y2, T600, cmap = 'binary_r')
#ax.tricontourf(lon_ch6[ix], lat_ch6[ix], tb_ch6[ix], levels = clevels, cmap = 'RdBu_r')
h5 = ax.contourf(X2, Y2, T600, cmap = 'binary_r', levels = linspace(144, 153, 19), alpha = 0.8)
ax.scatter(-lon_ch6[ix], lat_ch6[ix], s = 6, marker = 's', c = tb_ch6[ix], cmap = 'RdBu_r')
# add colorbar
cax = inset_axes(ax, width = "6%", height = "60%", loc = 3,
  bbox_to_anchor = (0.8, 0.08, 1, 1), bbox_transform=ax.transAxes)
colorbar(h5, cax = cax, ticks = [144, 147, 150, 153])
cax.set_xlabel('K')
cax.xaxis.set_label_position('top')
ax.text(85, 1, 'Temperature, 600 mbar', fontsize = 12)
ax.set_ylabel('PC latitude')

ax = axs[2,0]
#ax.pcolormesh(X2, Y2, nh3, cmap = 'binary_r')
#ax.tricontourf(lon_ch6[ix], lat_ch6[ix], tb_ch6[ix], levels = clevels, cmap = 'RdBu_r')
h6 = ax.contourf(X2, Y2, nh3, cmap = 'binary_r', levels = linspace(2, 15, 14), alpha = 0.8)
ax.scatter(-lon_ch6[ix], lat_ch6[ix], s = 6, marker = 's', c = tb_ch6[ix], cmap = 'RdBu_r')
cax = inset_axes(ax, width = "6%", height = "60%", loc = 3,
  bbox_to_anchor = (0.8, 0.08, 1, 1), bbox_transform=ax.transAxes)
colorbar(h6, cax = cax)
cax.set_xlabel('ppmv')
cax.xaxis.set_label_position('top')
ax.text(85, 1, 'Ammonia, 440 mbar', fontsize = 12)
ax.set_ylabel('PC latitude')
ax.set_xlabel('West longitude')

ax.set_ylim([-30., 5.])
#ax.set_xlim([-90., -30.])
ax.set_xlim([90., 30.])

savefig('VISIR_GEMINI_MWR.png', bbox_inches = 'tight', dpi = 600)
savefig('VISIR_GEMINI_MWR.tiff', bbox_inches = 'tight')
