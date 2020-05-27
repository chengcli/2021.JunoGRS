#! /usr/bin/env python3
from pylab import *
from astropy.io import fits
from numpy.random import multivariate_normal
from pyathena.utils import *
from pyfits.gpmcmc import *
import argparse, glob

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
parser.add_argument('--base',
  default = 'None',
  choices = ['None', 'PJ1345689', 'PJ07'],
  help = 'fit relative to base perijove'
  )
parser.add_argument('--case',
  default = '.',
  help = 'case folder for base perijove'
  )
args = vars(parser.parse_args())

if args['pj'] == 'PJ1345689':
  fname = 'Cassini_CIRS_Jupiter_2000.fits'
  hdul = fits.open(fname)
  temp = hdul[0].data
  pres = hdul[2].data
  ip = where((pres < 1.) & (pres > 0.1))[0][::3]
  irlat = Graphic2Centric(hdul[3].data)

else:
  fname = 'Gemini_TEXES_GRS_2017_product.fits'
  hdul = fits.open(fname)
  temp = hdul[0].data
  pres = hdul[2].data/1.E5
  ip = where((pres < 1.) & (pres > 0.1))[0][::3]

  lat1 = hdul[0].header['LATS']
  lat2 = hdul[0].header['LATN']
  dlat = hdul[0].header['LATSTEP']
  irlat = Graphic2Centric(arange(lat1, lat2 + dlat, dlat))
  temp = temp[37,:,:]

# baseline model
if args['base'] != 'None':
  fname = 'JunoMWR_coeffs_%s_m26-m14_avg.fits' % args['base']
  hdul = fits.open(fname)
  coeffs_base = hdul[0].data
  covs_base = hdul[1].data
  # correction for positive curvature
  for i in range(coeffs_base.shape[0]):
    for j in range(coeffs_base.shape[1]):
      if coeffs_base[i,j,2] > 0.:
        coeffs_base[i,j,2] = coeffs_base[i,j,1]/2.

# fitting model
fname = 'JunoMWR_coeffs_%s_m26-m14_avg.fits' % args['pj']
hdul = fits.open(fname)
coeffs = hdul[0].data
covs = hdul[1].data
lats = hdul[2].data
# correction for positive curvature
if args['base'] != 'None':
  for i in range(coeffs.shape[0]):
    for j in range(coeffs.shape[1]):
      if coeffs[i,j,2] > 0:
        coeffs[i,j,2] = coeffs_base[i,j,2]
else:
  for i in range(coeffs.shape[0]):
    for j in range(coeffs.shape[1]):
      if coeffs[i,j,2] > 0.:
        coeffs[i,j,2] = coeffs[i,j,1]/2.

for i, lat in enumerate(lats):
  if lat < 0: 
    obsname = 'mwr_%s_m%04.1f.obs' % (args['pj'], abs(lat))
  else: 
    obsname = 'mwr_%s_p%04.1f.obs' % (args['pj'], lat)

  if args['base'] != 'None':
    obsbase = args['case'] + '/' + re.sub(args['pj'], args['base'], obsname)
    obsfits = re.sub('.obs', '*.fits', obsbase)
    obsfits = glob.glob(obsfits)[0]

    par, val, lnp = find_best_fits(obsfits[:-5])
    par[2*len(par)//3:] *= 1.E3
    val = val.reshape(coeffs_base[i].shape)
    coeffs[i] -= coeffs_base[i] - val

  with open(obsname, 'w') as file:
    if args['base'] == 'None':
      file.write('# Target brightness temperatures:\n')
    else:
      file.write('# Target brightness temperatures relative to %s/%s:\n' % (args['base'],args['case']))
      s = '%12.5g'*len(par)
      file.write('# par:' + s % tuple(par) + '\n')
      file.write('# lnp:%12.5g\n' % lnp)
    file.write('%10d%10d\n' % (6, 3))
    for j in range(6):
      for k in range(3):
        file.write('%12.4f' % coeffs[i,j,k])
      file.write('\n')

    file.write('# Inverse covariance matrix:\n')
    file.write('%10d%10d%10d\n' % (6, 3, 3))

    # add calibration error and increase the error on limb darkening
    # increase covariance if not averaged observation
    if args['pj'] != 'PJ1345689':
      covs[i,:,3] *= 4
      covs[i,:,5] *= 8

    if args['base'] == 'None': 
      covs[i,0,0] += (0.01*coeffs[i,0,0])**2
      covs[i,0,3] += (0.1*coeffs[i,0,1])**2
      covs[i,0,5] += (0.1*coeffs[i,0,2])**2

    for j in range(6):
      # baseline error
      covs[i,j,0] += (0.01*coeffs[i,j,0])**2
      covs[i,j,3] += (0.03*coeffs[i,j,1])**2
      covs[i,j,5] += (0.05*coeffs[i,j,2])**2
      # add calibration error
      if args['base'] == 'None': 
        covs[i,j,0] += (0.01*(coeffs[i,j,0] - 300.))**2

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
    try:
      j = where(irlat < lat)[0][-1]
      # linear interpolation
      f1 = (lat - irlat[j])/(irlat[j+1] - irlat[j])
      f2 = (irlat[j+1] - lat)/(irlat[j+1] - irlat[j])
    except:
      j, f1, f2 = 0, 0., 1.
    file.write('%10d%10d\n' % (len(ip), 3))
    for p in ip:
      file.write('%12.4g%12.4g%12.4g\n' % (pres[p], f2*temp[j,p] + f1*temp[j+1,p], 2.))
  print('Observation file written to %s' % obsname)
