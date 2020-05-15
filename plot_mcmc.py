#! /usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import argparse 
from scipy.interpolate import interp1d
from pyathena.athena_read import athinput
from pylab import *
from pydrum.plots1d import *
from pyfits.gpmcmc import *

def plot_profiles(case):
  print('processing case %s ...' % case)
  # read mcmc results
  par, val, lnp, msk = read_mcmc_reduced(case)
  T1, p1, z1, theta1, theta1v, h2o1, nh3a, tb1, ptb1 = read_baseline(case)
  T3, nh3b = read_perturbed(case, msk)
  T4, p4, rho4, theta4, theta4v, h2o4, nh3d, tb4 = read_rectified(case, msk)

  # averaged rectified pressure and density
  p4a = mean(p4, axis = 0)
  rho4a = mean(rho4, axis = 0)

  # target and other tp data
  inp = athinput(case + '.inp')
  obsfile = inp['problem']['obsfile']
  target = genfromtxt(obsfile, skip_header = 2, max_rows = 6)
  try:
    tpdata = genfromtxt(obsfile, skip_header = 36, max_rows = 6)
  except:
    tpdata = []

  # interpolation function
  lnpTfunc = interp1d(log(p1[::-1]), T1[::-1], bounds_error = False)
  lnpNH3func = interp1d(log(p1[::-1]), nh3a[::-1], bounds_error = False)

  # plot profiles
  fig, axs = subplots(1, 4, figsize = (14,8), sharey = True)
  subplots_adjust(wspace = 0.08, hspace = 0.08)
  ax = axs[0]
  ax.plot(T1, p1, 'C1')
  ax.scatter(tb1, ptb1, s = 50, ec = 'C4', fc = 'none')
  ax.set_xlabel('T (K)', fontsize = 15)
  ax.set_xlim([100., 1500.])
  ax.set_ylabel('Pressure (bar)', fontsize = 15)
  ax.set_yscale('log')
  ax.set_ylim([1.E3, 0.1])

  ax2 = ax.twiny()
  ax2.plot(h2o1, p1, 'C0')
  ax2.plot(nh3a, p1, 'C2')
  ax2.set_xscale('log')
  ax2.set_xlim([0.1, 100.])
  ax2.set_xlabel('mmr (g/kg)')

  ax = axs[1]
  ax.plot([0., 0.], [1.E3, 0.1], '--', color = '0.7')
  #ax.plot(mean(nh3b, axis = 0) - lnpNH3func(log(p4a)), p4a, 'C2--')
  PlotProfile(ax, nh3d - lnpNH3func(log(p4a)), p4a, 'C2')
  ax.set_xlabel('Perturbed NH3 (g/kg)', fontsize = 15)

  ax = axs[2]
  ax.plot([0., 0.], [1.E3, 0.1], '--', color = '0.7')
  #ax.plot(mean(T3, axis = 0) - lnpTfunc(log(p4a)), p4a, 'C1--')
  PlotProfile(ax, T4 - lnpTfunc(log(p4a)), p4a, 'C1')
  ax.scatter(mean(tb4, axis = 0) - tb1, ptb1, marker = 'x', s = 30, color = 'C3')
  ax.errorbar(target[:,0] - tb1, ptb1, xerr = 0.02*target[:,0], fmt = 'none', color = 'C4', capsize = 5)
  if len(tpdata) > 0:
    ax.errorbar(tpdata[:,1] - lnpTfunc(log(tpdata[:,0])), tpdata[:,0], xerr = tpdata[:,2],
      fmt = 'none', color = 'C6', capsize = 5)
  ax.set_xlabel('Perturbed T (K)', fontsize = 15)

  ax = axs[3]
  ax.plot(theta1, p1, '--', color = '0.7')
  PlotProfile(ax, theta4, p4a, 'C1')
  ax.plot(mean(theta4v, axis = 0), p4a, 'C1--')
  ax.set_xlabel('Potential T (K)', fontsize = 15)
  ax.set_xlim([150., 190.])

  ax2 = ax.twiny()
  PlotProfile(ax2, h2o4, p4a, 'C0')
  ax2.set_xlim([-5., 50.])
  ax2.set_xlabel('mmr (g/kg)')

  savefig('mcmc_%s.png' % case, bbox_inches = 'tight')
  savefig('mcmc_%s.pdf' % case, bbox_inches = 'tight')
  close(fig)

  return T4 - lnpTfunc(log(p4a)), h2o4, nh3d, p4a, rho4a

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-c',
      default = 'none',
      help = 'name of the first netcdf file'
      )
  args = vars(parser.parse_args())
  plot_profiles(args['c'])
