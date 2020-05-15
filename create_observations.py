#! /usr/bin/env python3
from pylab import *
from astropy.io import fits
from numpy.random import multivariate_normal
import argparse

def TbFunction(a,b,c,mu):
  return a + b*(1.-mu) + c*(1.-mu)**2

def Graphic2Centric(lat):
  re = 71398.E3
  rp = 67040.E3
  return arctan(tan(lat/180.*pi)*(rp/re)**2)/pi*180.

parser = argparse.ArgumentParser()
parser.add_argument('--pj',
  default = 'PJ1345689',
  choices = ['PJ1345689', 'PJ07'],
  help = 'select perijove'
  )
args = vars(parser.parse_args())

fname = 'Gemini_TEXES_GRS_2017_product.fits'
hdul = fits.open(fname)
temp = hdul[0].data
pres = hdul[2].data/1.E5
ip = where((pres < 1.) & (pres > 0.1))[0][::3]

lat1 = hdul[0].header['LATS']
lat2 = hdul[0].header['LATN']
dlat = hdul[0].header['LATSTEP']
irlat = Graphic2Centric(arange(lat1, lat2 + dlat, dlat))

#print(temp[37,12,ip])
temp = mean(temp[:10,:,:], axis = 0)

fname = 'JunoMWR_coeffs_%s_m24-m16_avg.fits' % args['pj']
hdul = fits.open(fname)

# read 1) brightness temperature coefficients
#      2) covariance matrix
#      3) latitudes
coeffs = hdul[0].data
covs = hdul[1].data
lats = hdul[2].data

for i, lat in enumerate(lats):
  if (lat < irlat.min()) or (lat > irlat.max()):
    continue
  if lat < 0: 
    obsname = 'mwr_%s_m%04.1f.obs' % (args['pj'], abs(lat))
  else: 
    obsname = 'mwr_%s_p%04.1f.obs' % (args['pj'], lat)
  with open(obsname, 'w') as file:
    file.write('# Target brightness temperatures:\n')
    file.write('%10d%10d\n' % (6, 3))
    for j in range(6):
      for k in range(3):
        file.write('%12.4f' % coeffs[i,j,k])
      file.write('\n')

    file.write('# Inverse covariance matrix:\n')
    file.write('%10d%10d%10d\n' % (6, 3, 3))
    for j in range(6):
      # calibration error
      tb0 = TbFunction(*coeffs[i,j,:], 1.)
      covs[i,j,0] += (0.02*(tb0 - 300.))**2
      if j == 0 or j == 5:
        covs[i,j,3] *= 64
        covs[i,j,5] *= 256
      else:
        covs[i,j,3] *= 16
        covs[i,j,5] *= 64
      cov = array([[covs[i,j,0], covs[i,j,1], covs[i,j,2]],
                   [covs[i,j,1], covs[i,j,3], covs[i,j,4]],
                   [covs[i,j,2], covs[i,j,4], covs[i,j,5]]])
      cov = linalg.inv(cov)
      for k1 in range(3):
        for k2 in range(3):
          file.write('%12.4g' % cov[k1,k2])
        file.write('\n')
      file.write('\n')

    file.write('# Additional TP data:\n')
    j = where(irlat < lat)[0][-1]
    # linear interpolation
    f1 = (lat - irlat[j])/(irlat[j+1] - irlat[j])
    f2 = (irlat[j+1] - lat)/(irlat[j+1] - irlat[j])
    file.write('%10d%10d\n' % (len(ip), 3))
    for p in ip:
      file.write('%12.4g%12.4g%12.4g\n' % (pres[p], f2*temp[j,p] + f1*temp[j+1,p], 1.))
  print('Observation file written to %s' % obsname)
