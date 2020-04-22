#! /usr/bin/env python3
#import matplotlib 
#matplotlib.use('Agg')
import argparse
from mwrlib3 import *
from pylab import *

#matplotlib.pyplot.style.use('dark_background')

def Graphic2Centric(re, rp, lat):
  return arctan(tan(lat/180.*pi)*(rp/re)**2)/pi*180.

def Centric2Graphic(re, rp, lat):
  return arctan(tan(lat/180.*pi)*(re/rp)**2)/pi*180.

# centric to graphic
re = 71398.E3
rp = 67040.E3

version = 'v1906xx'
chs = [0,1,2,3,4,5]
mu = cos(pi/4.)
ylims = [[790.,920.],[410.,490.],[290.,360.],[218.,270.],[170.,220.],[134.,155.]]
pad = [10., 5., 5., 2., 2., 1.]
freq = [0.6, 1.25, 2.6, 5.2, 10., 22.]
latmin, latmid, latmax = -26, -20, -14
#latmin, latmid, latmax = -27.1, -20, -7.2

figure(1, figsize = (8, 8))
ax = axes()
tb_anomaly = []
for ch in chs:
  # average
  lat, lon, coeff, cov = ReadAveragedCoefficients(version, 'PJ1345689', ch,
    latmin = -80, latmax = 80.)
  #lat = Centric2Graphic(re, rp, lat)

  # std tb and ld
  nsample = 10000
  std_tb, std_ld = [], []
  for i in range(len(lat)):
    sample = multivariate_normal(coeff[i], cov[i], nsample)
    tb0 = TbFunction(sample[:,0], sample[:,1], sample[:,2], 1.)
    tb1 = TbFunction(sample[:,0], sample[:,1], sample[:,2], mu)
    ld = 100.*(tb0 - tb1)/tb0
    std_tb.append(std(tb0, axis = 0))
    std_ld.append(std(ld, axis = 0))
  std_tb = array(std_tb)
  std_ld = array(std_ld)

  # mean tb and ld
  tb0 = TbFunction(coeff[:,0], coeff[:,1], coeff[:,2], 1.)
  tb1 = TbFunction(coeff[:,0], coeff[:,1], coeff[:,2], mu)
  ld = 100.*(tb0 - tb1)/tb0

  if ch == 0:
    ax.set_xlabel('Planetocentric Latitude', fontsize = 18)

  # PJ7
  lat7, lon7, coeff7, _ = ReadCoefficients(version, 7, ch,
    latmin = -80, latmax = 80.)
  #lat7 = Centric2Graphic(re, rp, lat7)

  tb_anomaly.append(coeff7[:,0] - tb0)


  ix = (lat > -80) & (lat < 80.)
  print(tb0[ix])
  #print coeff7[ix,0]

tb_anomaly = array(tb_anomaly)
print(tb_anomaly.min(), tb_anomaly.max())
X, Y = meshgrid(lat, freq)

ix = (lat > -80) & (lat < 80.)
print('# Planetographic Latitudes:')
for lat in lat[ix]:
  print('%12.2f' % lat,)
print()
print('# Frequencies (GHz):')
for f in freq:
  print('%12.2f' % f,)
print()
print('# GRS Brightness Temperature anomaly:')
for i in range(tb_anomaly[:,ix].shape[0]):
  for j in range(tb_anomaly[:,ix].shape[1]):
    print('%12.2f' % tb_anomaly[:,ix][i,j],)
  print()
print('# Mean Brightness Temperature anomaly:')
for i in range(len(freq)):
  print('%12.2f' % mean(tb_anomaly[i, ix]),)
print()

#ax.fill_between(array([-18, latmax]), freq[0], freq[-1], color = '0.6', alpha = 0.5)
ax.plot([latmid, latmid], [freq[0], freq[-1]], 'C1--', linewidth = 2)
ax.plot([-16, -16], [freq[0], freq[-1]], 'C1-', linewidth = 2)
ax.plot([-24, -24], [freq[0], freq[-1]], 'C1-', linewidth = 2)

h1 = ax.contour(X, Y, tb_anomaly, linspace(-11, -1, 5), colors = 'b', linestyles = '-')
#h = ax.contourf(X, Y, tb_anomaly, linspace(-11, -1, 6), cmap = 'Blues')

h2 = ax.contour(X, Y, tb_anomaly, linspace(1, 26, 11), colors = 'r')
#h = ax.contourf(X, Y, tb_anomaly, linspace(1, 27, 14), cmap = 'Reds')

ax.set_xlim([latmin, latmax])
ax.set_xticks(arange(latmin, latmax + 2, 2))
ax.set_yscale('log')
ax.set_yticks(freq)
ax.set_yticklabels(freq)
ax.set_ylabel('Frequency (GHz)', fontsize = 18)

#clabel(h1, fontsize = 12, inline = 1, fmt = '%.1f', manual = True)
clabel(h1, fontsize = 12, inline = 1, fmt = '%.1f')
#clabel(h2, fontsize = 12, inline = 1, fmt = '%.1f', manual = True)
clabel(h2, fontsize = 12, inline = 1, fmt = '%.1f')

savefig('grs_temperature_anomaly.png', bbox_inches = 'tight', dpi = 600)
#savefig('grs_temperature_anomaly.eps', bbox_inches = 'tight')
#show()
