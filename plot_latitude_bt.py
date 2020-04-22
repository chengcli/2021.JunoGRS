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

# graphic to centric
re = 71398.E3
rp = 67040.E3

version = 'v1906xx'
chs = [0,1,2,3,4,5]
mu = cos(pi/4.)
ylims = [[790.,920.],[410.,490.],[290.,360.],[218.,270.],[170.,220.],[134.,155.]]
pad = [10., 5., 5., 2., 2., 1.]
freq = [0.6, 1.25, 2.6, 5.2, 10., 22.]
uwind = genfromtxt('u_vs_lat.jupiter.cassini')
vort = zeros((len(uwind)-1, 2))

vort[:,0] = Graphic2Centric(re, rp, 0.5*(uwind[1:,0] + uwind[:-1,0]))
vort[:,1] = uwind[:-1,1] - uwind[1:,1]

# three point runing average
for i in range(1, len(vort)-1):
  vort[i,1] = (vort[i-1,1] + vort[i,1] + vort[i+1,1])/3.

fig, axs = subplots(6, 1, figsize = (16, 10), sharex = True)
subplots_adjust(hspace = 0.06)
for ch in chs:
  # average
  lat, lon, coeff, cov = ReadAveragedCoefficients(version, 'PJ1345689', ch)
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

  # reference
  ax = axs[5-ch]
  ax.plot(sin(lat/180.*pi), tb0, 'C0', linewidth = 1)
  ax.fill_between(sin(lat/180.*pi), tb0 - std_tb, tb0 + std_tb, color = 'C2', alpha = 0.5)
  ax.set_ylim(ylims[ch])

  if ch == 0:
    ax.set_xlabel('Planetocentric latitude', fontsize = 18)
  ax.set_ylabel('%g GHz' % freq[ch], fontsize = 12)

  # PJ7
  lat, lon, coeff, cov = ReadCoefficients(version, 7, ch)
  #lat = Centric2Graphic(re, rp, lat)
  ax.plot(sin(lat/180.*pi), coeff[:,0], 'C4', linewidth = 2)

  # GRS range
  #latmin, latmax = -26, -18
  latmin, latmax = -24, -16
  ax.plot([sin(latmin/180.*pi), sin(latmin/180.*pi)],
    [ylims[ch][0]+pad[ch], ylims[ch][1]-pad[ch]], 'C1', linewidth = 1)
  ax.plot([sin(latmax/180.*pi), sin(latmax/180.*pi)], 
    [ylims[ch][0]+pad[ch], ylims[ch][1]-pad[ch]], 'C1', linewidth = 1)

  #latmid = -22
  latmid = -20
  ax.plot([sin(latmid/180.*pi), sin(latmid/180.*pi)],
    [ylims[ch][0]+pad[ch], ylims[ch][1]-pad[ch]], 'C1--', linewidth = 1)

  # shade vorticity
  imin = 0
  while imin < len(vort):
    if vort[imin,1]*vort[imin,0] < 0:
      imin += 1
    else:
      imax = imin
      while imax < len(vort):
        if vort[imax,1]*vort[imax,0] < 0:
          if vort[imax,0] - vort[imin,0] > 2.:
            if ch == 5 and vort[imin,0] > -50 and vort[imax,0] < 50:
              if vort[imin,0] < 0.:
                ax.text(sin((vort[imin,0]-1.)/180.*pi), ylims[ch][1] + 2*pad[ch], '$\otimes$', fontsize = 20)
              else:
                ax.text(sin((vort[imin,0]-1.)/180.*pi), ylims[ch][1] + 2*pad[ch], '$\odot$', fontsize = 20)
              if vort[imax,0] < 0.:
                ax.text(sin((vort[imax,0]-1.)/180.*pi), ylims[ch][1] + 2*pad[ch], '$\odot$', fontsize = 20)
              else:
                ax.text(sin((vort[imax,0]-1.)/180.*pi), ylims[ch][1] + 2*pad[ch], '$\otimes$', fontsize = 20)
            print(vort[imin,0], vort[imax,0])
            ax.fill_between(sin(vort[imin:imax,0]/180.*pi), ylims[ch][0]+pad[ch],
              ylims[ch][1]-pad[ch], alpha = 0.5, color = '0.6', step = 'mid', linewidth = 0)
          imin = imax
          #print vort[imin,0]
          break
        else:
          imax += 1
      if imax == len(vort): break;

ax.set_xticks(sin(arange(-80., 100., 10.)/180.*pi))
ax.set_xticklabels(arange(-80, 100, 10), fontsize = 15)
ax.set_xlim([sin(-pi/3), sin(pi/3)])

savefig('pj7_pclatitude_bt.png', bbox_inches = 'tight', dpi = 600)
#savefig('pj7_pclatitude_bt.eps', bbox_inches = 'tight')
#show()
