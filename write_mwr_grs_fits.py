#! /usr/bin/env python3
from pylab import *
from astropy.io import fits

version = 'v1906xx'
chs = [0,1,2,3,4,5]
latmin, latmax = -24, -16 # PC latitude
freq = [0.6, 1.25, 2.6, 5.2, 10., 22.] # frequency in GHz

nlat = 13
coeff = zeros((nlat, 6, 3)) # 13 lats, 6 freqs, 3 coeffs
cov = zeros((nlat, 6, 6)) # 13 lats, 6 freq, 3 covariance matrix coeffs
for ch in chs:
  data = genfromtxt('JunoMWR_coeffs_%s_%s_CH%d.txt' % (version, 'PJ1345689', ch+1))
  #data = genfromtxt('JunoMWR_coeffs_%s_%s_CH%d.txt' % (version, 'PJ07', ch+1))
  i1 = searchsorted(data[:,0], latmin, 'left')
  i2 = searchsorted(data[:,0], latmax, 'right')
  lat = data[i1:i2,0]
  lon = data[i1:i2,1]
  coeff[:,ch,:] = data[i1:i2,2:5]
  cov[:,ch,:] = data[i1:i2,5:11]

hdu = fits.PrimaryHDU(coeff)
hdu.header['CREATOR'] = ('Cheng Li', 'file creator')
hdu.header['VAR'] = ('Tb', 'brightness temperature coefficients')
hdu.header['LATS'] = (latmin, 'PC latitude')
hdu.header['LATN'] = (latmax, 'PC latitude')
hdu.header['CH1'] = (freq[0], 'MWR CH1 frequency')
hdu.header['CH2'] = (freq[1], 'MWR CH2 frequency')
hdu.header['CH3'] = (freq[2], 'MWR CH3 frequency')
hdu.header['CH4'] = (freq[3], 'MWR CH4 frequency')
hdu.header['CH5'] = (freq[4], 'MWR CH5 frequency')
hdu.header['CH6'] = (freq[5], 'MWR CH6 frequency')

hdu.header['NAXIS1'] = (3, 'a,b,c coefficients')
hdu.header['NAXIS2'] = (6, 'CH1-6')
hdu.header['NAXIS3'] = (nlat, 'PC latitudes')

hlat = fits.ImageHDU(lat)
hlat.header['VAR'] = ('lat', 'PC latitude coordinates')

hlon = fits.ImageHDU(lon)
hlon.header['VAR'] = ('lon', 'SystemIII longitude coordinates')

hcov = fits.ImageHDU(cov)
hcov.header['VAR'] = ('cov', 'covariance matrix coefficients')

hcov.header['NAXIS1'] = (6, 'covariance matrix coefficients')
hcov.header['NAXIS2'] = (6, 'CH1-6')
hcov.header['NAXIS3'] = (nlat, 'PC latitudes')

hdul = fits.HDUList([hdu, hcov, hlat, hlon])
hdul.writeto('JunoMWR_coeffs_PJ1345689_m24-m16_avg.fits', overwrite = True)
#hdul.writeto('JunoMWR_coeffs_PJ07_m24-m16_avg.fits', overwrite = True)
