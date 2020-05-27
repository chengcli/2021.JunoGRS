#! /usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
from scipy.interpolate import interp1d
from pyathena.utils import *
from pydrum.plots1d import *
from pylab import *
from plot_mcmc import *
from astropy.io import fits
from multiprocessing import Pool

lats, tps, h2os, nh3s, pres, rhos = [], [], [], [], [], []

files = glob.glob('mwr*.fits')
cases = [fname[:-5] for fname in files]
with Pool(3) as p:
  results = p.map(plot_profiles, cases)

for i, v in enumerate(results):
  lats.append(parse_signed(cases[i]))
  tps.append(average(v[0], axis = 0, weights = v[5]))
  h2os.append(average(v[1], axis = 0, weights = v[5]))
  nh3s.append(average(v[2], axis = 0, weights = v[5]))
  pres.append(v[3])
  rhos.append(v[4])

print('making 2d maps ...')
lats = array(lats)
tps = array(tps)
h2os = array(h2os)
nh3s = array(nh3s)/2.7*350. # g/kg -> ppm
pres = array(pres)
rhos = array(rhos)

ix = argsort(lats)
lats = lats[ix]
tps = tps[ix,:]
nh3s = nh3s[ix,:]
pres = pres[ix,:]
rhos = rhos[ix,:]

output = re.sub('_[mp]\d+\.\d*', '', cases[0])
output = re.sub('mwr', 'profile', output)

# interpolation to constant pressure
pavg = mean(pres, axis = 0)

for i in range(len(lats)):
  Tfunc = interp1d(log(pres[i,::-1]), tps[i,::-1], 
    bounds_error = False, fill_value = 'extrapolate')
  tps[i,:] = Tfunc(log(pavg))
  NH3func = interp1d(log(pres[i,::-1]), nh3s[i,::-1], 
    bounds_error = False, fill_value = 'extrapolate')
  nh3s[i,:] = NH3func(log(pavg))

# ammonia map
X, Y = meshgrid(lats, pavg)
figure(1, figsize = (10, 8))
ax = axes()
h = ax.contourf(X, Y, nh3s.T, linspace(0, 400, 21), cmap = 'inferno')
c = colorbar(h)
c.ax.invert_yaxis()
ax.set_yscale('log')
ax.set_ylim([100., 0.3])
ax.set_xlabel('PC latitude', fontsize = 15)
ax.set_ylabel('Pressure (bar)', fontsize = 15)
savefig('%s_ammonia2d.png' % output, bbox_inches = 'tight')
savefig('%s_ammonia2d.pdf' % output, bbox_inches = 'tight')
close()

# temperature anomaly map
figure(1, figsize = (10, 8))
ax = axes()
#h1 = ax.contour(X, Y, tps.T, [-5, -3, -1], colors = 'b')
#clabel(h1, fontsize = 12, inline = 1, fmt = '%.1f')
#h2 = ax.contour(X, Y, tps.T, linspace(1, 21, 11), colors = 'r')
#clabel(h2, fontsize = 12, inline = 1, fmt = '%.1f')
h = ax.contourf(X, Y, tps.T, linspace(-5, 11, 17), cmap = 'OrRd')
#ax.contour(X, Y, tps.T, [0.], colors = '0.7', linewidths = 4)
c = colorbar(h)
ax.set_yscale('log')
ax.set_ylim([100., 0.3])
ax.set_xlabel('PC latitude', fontsize = 15)
ax.set_ylabel('Pressure (bar)', fontsize = 15)
savefig('%s_temp2d.png' % output, bbox_inches = 'tight')
savefig('%s_temp2d.pdf' % output, bbox_inches = 'tight')
close()

# save results to fits file
#result = read_baseline(cases[0])
with Pool(8) as p:
  results = p.map(read_baseline, array(cases)[ix])

T1, P1, Z1 = [], [], []
for v in results:
  T1.append(v[0])
  P1.append(v[1])
  Z1.append(v[2])
T1 = array(T1)
P1 = array(P1)
Z1 = array(Z1)

h1 = fits.PrimaryHDU(nh3s/350.*2.7) # ppm -> g/kg
h1.header['CREATOR'] = ('Cheng Li', 'file creator')
h1.header['VAR'] = ('qNH3', 'ammonia mass mixing ratio (g/kg)')
h1.header['LATS'] = (lats.min(), 'PC latitude south')
h1.header['LATN'] = (lats.max(), 'PC latitude north')
h1.header['PTOP'] = (pres.min(), 'top pressure (bar)')
h1.header['PBOT'] = (pres.max(), 'bottom pressure (bar)')

h2 = fits.ImageHDU(h2os)
h2.header['VAR'] = ('qH2O', 'water mass mixing ratio (g/kg)')

h3 = fits.ImageHDU(tps)
h3.header['VAR'] = ('Tp', 'temperature anomaly (K)')

h4 = fits.ImageHDU(lats)
h4.header['VAR'] = ('lat', 'PC latitude coordinates')

h5 = fits.ImageHDU(pres)
h5.header['VAR'] = ('pres', 'pressure (bar)')

h6 = fits.ImageHDU(rhos)
h6.header['VAR'] = ('rho', 'density (kg/m^3)')

h7 = fits.ImageHDU(T1)
h7.header['VAR'] = ('T1', 'reference temperature (K)')

h8 = fits.ImageHDU(P1)
h8.header['VAR'] = ('P1', 'reference pressure (bar)')

h9 = fits.ImageHDU(Z1)
h9.header['VAR'] = ('Z1', 'reference height (km)')

hdul = fits.HDUList([h1, h2, h3, h4, h5, h6, h7, h8, h9])
hdul.writeto('%s.fits' % output, overwrite = True)
