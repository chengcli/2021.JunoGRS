#! /usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
from scipy.interpolate import interp1d
#from pyathena.utils import *
#from pydrum.plots1d import *
from pylab import *
#from plot_mcmc import *
from astropy.io import fits
#from pyfits.gpmcmc import *
from athena_read import athinput
import glob
#from multiprocessing import Pool

def find_best_fits(case):
  hdul = fits.open('%s.fits' % case)

  nstep, nwalker, ndim = hdul[0].data.shape
  nstep, nwalker, nval = hdul[1].data.shape

  par = hdul[0].data.reshape(-1, ndim)
  val = hdul[1].data.reshape(-1, nval)
  lnp = hdul[2].data.flatten()

  ix = lnp.argmax()

  return par[ix], val[ix], lnp[ix]

def plot_map_pair(ax1, ax2, lats, pres, nh3s, tps):
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
  h1 = ax1.contourf(X, Y, nh3s.T, linspace(40, 360, 17),
    cmap = 'inferno', extend = 'min')
  ax1.plot([-24, -24], [100., 0.3], 'w', linewidth = 1)
  ax1.plot([-16, -16], [100., 0.3], 'w', linewidth = 1)
  ax1.plot([-20, -20], [100., 0.3], 'w--', linewidth = 1)
  ax1.set_yscale('log')
  ax1.set_xlim([min(lats), max(lats)])
  ax1.set_ylim([100., 0.3])
  ax1.tick_params(axis = 'y', labelsize = 18)
  ax1.tick_params(axis = 'x', labelsize = 18, labelbottom = False)

  h2 = ax2.contourf(X, Y, tps.T, linspace(-5, 11, 17), cmap = 'OrRd')
  #h2 = ax2.contourf(X, Y, tps.T, linspace(-10, 10, 11), cmap = 'RdBu_r', extend = 'max')
  ax2.plot([-24, -24], [100., 0.3], 'w', linewidth = 1)
  ax2.plot([-16, -16], [100., 0.3], 'w', linewidth = 1)
  ax2.plot([-20, -20], [100., 0.3], 'w--', linewidth = 1)
  ax2.set_yscale('log')
  ax2.set_xlim([min(lats), max(lats)])
  ax2.set_ylim([100., 0.3])
  ax2.tick_params(axis = 'y', labelsize = 18, labelleft = False)
  ax2.tick_params(axis = 'x', labelsize = 18, labelbottom = False)

  return h1, h2

def get_chi2(folder):
  fname = folder + '/mwr*.fits'
  files = sort(glob.glob(fname))
  print(files)
  cases = [fname[:-5] for fname in files]

  chi2lst = []
  for case in cases:
    par, val, lnp = find_best_fits(case)
    val = val.reshape((-1, 3))

    inp = athinput(case + '.inp')
    obsfile = '../' + inp['problem']['obsfile']
    try:
      target = genfromtxt(obsfile, skip_header = 2, max_rows = 6)
      error = genfromtxt(obsfile, skip_header = 10, max_rows = 18)
    except:
      target = genfromtxt(obsfile, skip_header = 4, max_rows = 6)
      error = genfromtxt(obsfile, skip_header = 12, max_rows = 18)

    error = error.reshape((-1, 3, 3))

    chi2 = 0.
    for i in range(len(target)):
      #chi2 += 0.5*dot(val[i,0] - target[i,0], dot(error[i,0,0], val[i,0] - target[i,0]))
      chi2 += abs(val[i,0] - target[i,0])
    chi2lst.append(chi2/6)
  return chi2lst[::-1]


if __name__ == '__main__':
  fig, axs = subplots(5, 2, figsize = (24, 16),
    gridspec_kw = {'height_ratios': [2, 4, 4, 4, 0.5]})
  subplots_adjust(hspace = 0.08, wspace = 0.04)
  fname = 'profile_PJ07_2.7x20.x168.fits'

  hdul = fits.open('../PJ07BOOSTV3/' + fname)
  nh3s = hdul[0].data/2.7*350.
  tps = hdul[2].data
  lats = hdul[3].data
  pres = hdul[4].data
  plot_map_pair(axs[1,0], axs[1,1], lats, pres, nh3s, tps)
  axs[1,0].text(-25.5, 80., '(a)', fontsize = 24)
  axs[1,1].text(-25.5, 80., '(b)', fontsize = 24)

  hdul = fits.open('../PJ07BOOST.Tstd2/' + fname)
  nh3s = hdul[0].data/2.7*350.
  tps = hdul[2].data
  lats = hdul[3].data
  pres = hdul[4].data
  plot_map_pair(axs[2,0], axs[2,1], lats, pres, nh3s, tps)
  axs[2,0].text(-25.5, 80., '(c)', fontsize = 24)
  axs[2,1].text(-25.5, 80., '(d)', fontsize = 24)

  hdul = fits.open('../PJ07BOOST.Tstd1/' + fname)
  nh3s = hdul[0].data/2.7*350.
  tps = hdul[2].data
  lats = hdul[3].data
  pres = hdul[4].data
  h1, h2 = plot_map_pair(axs[3,0], axs[3,1], lats, pres, nh3s, tps)
  axs[3,0].text(-25.5, 80., '(e)', fontsize = 24)
  axs[3,1].text(-25.5, 80., '(f)', fontsize = 24)

  colorbar(h1, cax = axs[4,0], orientation = 'horizontal')
  colorbar(h2, cax = axs[4,1], orientation = 'horizontal')
  axs[4,0].set_xlabel('NH$_3$ mixing ratio (ppmv)', fontsize = 24)
  axs[4,0].tick_params(axis = 'x', labelsize = 18)
  axs[4,1].set_xlabel("T' (K)", fontsize = 24)
  axs[4,1].tick_params(axis = 'x', labelsize = 18)

  #axs[3,0].set_xlabel('Planetocentric latitude', fontsize = 15)
  #axs[3,1].set_xlabel('Planetocentric latitude', fontsize = 15)
  axs[1,0].set_ylabel('$\sigma(T) = 5 K$\n Pressure (bar)', fontsize = 20)
  axs[2,0].set_ylabel('$\sigma(T) = 2 K$\n Pressure (bar)', fontsize = 20)
  axs[3,0].set_ylabel('$\sigma(T) = 1 K$\n Pressure (bar)', fontsize = 20)

  chi2 = get_chi2('../PJ07BOOSTV3')
  axs[0,0].plot(lats, chi2, 'C0')
  axs[0,1].plot(lats, chi2, 'C0')

  chi2 = get_chi2('../PJ07BOOST.Tstd2')
  axs[0,0].plot(lats, chi2, 'C1')
  axs[0,1].plot(lats, chi2, 'C1')

  chi2 = get_chi2('../PJ07BOOST.Tstd1')
  axs[0,0].plot(lats, chi2, 'C2')
  axs[0,1].plot(lats, chi2, 'C2')

  for i in [0,1]:
    axs[0,i].legend(['5k', '2k', '1k'], fontsize = 16, loc = 2, ncol = 3)
    axs[0,i].set_xlim([min(lats), max(lats)])
    axs[0,i].set_ylim([0., 7.5])
    axs[0,i].xaxis.tick_top()
    axs[0,i].xaxis.set_tick_params(labelsize = 18)
    axs[0,i].yaxis.set_tick_params(labelsize = 18)
    axs[0,i].xaxis.set_label_position('top')
    axs[0,i].set_xlabel('PC Latitude', fontsize = 24)

  axs[0,0].set_ylabel('Residual (K)', fontsize = 20)
  axs[0,1].set_yticklabels([])

  savefig('../figs/PJ07BOOST_result_v2.png', bbox_inches = 'tight')
  savefig('../figs/PJ07BOOST_result_v2.pdf', bbox_inches = 'tight')
  close()
