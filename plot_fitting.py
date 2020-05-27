#! /usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import argparse, glob
from pylab import *
from pyathena.utils import *
from pyfits.gpmcmc import *

parser = argparse.ArgumentParser()
parser.add_argument('--case',
  default = '0519a',
  help = 'case folder for base perijove'
  )
args = vars(parser.parse_args())

if __name__ == '__main__':
  lats, targets, errs, vals, lnps = read_all_fits(args['case'], 'mwr')

  fig, axs = subplots(2, 2, figsize = (12, 10), sharex = True)
  subplots_adjust(hspace = 0.08)

  # coeff a
  ax = axs[0,0]
  for i in range(1,6):
    j = 3*i
    ax.plot(lats, vals[:,j], color = 'C%d' % i)
    ax.fill_between(lats, targets[:,j] - errs[:,i,0], targets[:,j] + errs[:,i,0],
      facecolor = '0.7', edgecolor = None)
  ax.set_ylim([500., 100.])
  ax.set_ylabel('C0 (K)')

  # coeff b
  ax = axs[1,0]
  for i in range(1,6):
    j = 3*i+1
    ax.plot(lats, vals[:,j], color = 'C%d' % i)
    ax.fill_between(lats, targets[:,j] - errs[:,i,4], targets[:,j] + errs[:,i,4],
      facecolor = '0.7', edgecolor = None)
  ax.set_ylim([-150., 0.])
  ax.set_xlabel('PC latitude')
  ax.set_ylabel('C1 (K)')

  # coeff a, CH1
  ax = axs[0,1]
  ax.plot(lats, vals[:,0], color = 'C0')
  ax.fill_between(lats, targets[:,0] - errs[:,0,0], targets[:,0] + errs[:,0,0],
    facecolor = '0.7', edgecolor = None)

  # coeff b, CH1
  ax = axs[1,1]
  ax.plot(lats, vals[:,1], color = 'C0')
  ax.fill_between(lats, targets[:,1] - errs[:,0,4], targets[:,1] + errs[:,0,4],
    facecolor = '0.7', edgecolor = None)
  ax.set_xlabel('PC latitude')

  savefig('%s/tbld_fit.png' % args['case'])
  savefig('%s/tbld_fit.pdf' % args['case'])
  close()

  # max lnp
  figure(2, figsize = (10, 8))
  ax = axes()
  ax.plot(lats, lnps, color = 'k')
  savefig('%s/lnp_fit.png' % args['case'])
  savefig('%s/lnp_fit.pdf' % args['case'])
  close()
