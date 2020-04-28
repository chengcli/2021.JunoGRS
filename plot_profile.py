#! /usr/bin/env python3
#import matplotlib
#matplotlib.use('Agg')
import argparse, h5py
from netCDF4 import Dataset
from pydrum.plots1d import PlotProfile, DrawPressureAxis
from scipy.interpolate import interp1d
from glob import os, glob
from pylab import *

choices = []
for x in glob('mwr-*.nc'):
  choices.append(x.replace('mwr-','').replace('-main','').replace('-rad','').replace('.nc',''))
choices = sort(unique(choices))

parser = argparse.ArgumentParser()
parser.add_argument('--case',
    default = 'none',
    choices = choices,
    help = 'name of the first netcdf file'
    )
args = vars(parser.parse_args())

data = Dataset('mwr-%s-main.nc' % args['case'], 'r')
data.set_always_mask(False)
T1 = data['temp'][0,0,0,:]
T3 = data['temp'][0,0,2,:]
T4 = data['temp'][0,0,3,:]

theta1 = data['theta'][0,0,0,:]
theta4 = data['theta'][0,0,3,:]
theta4v = data['thetav'][0,0,3,:]

h2o1 = data['vapor1'][0,0,0,:]*1.E3
h2o4 = data['vapor1'][0,0,3,:]*1.E3

nh3a = data['vapor2'][0,0,0,:]*1.E3
nh3b = data['vapor2'][0,0,1,:]*1.E3
nh3d = data['vapor2'][0,0,3,:]*1.E3

x1 = data['x1'][:]/1.E3
p1 = data['press'][0,0,0,:]/1.E5
p4 = data['press'][0,0,3,:]/1.E5
TPfunc = interp1d(T1[::-1], p1[::-1], bounds_error = False)
lnpTfunc = interp1d(log(p1[::-1]), T1[::-1], bounds_error = False)
lnpNH3func = interp1d(log(p1[::-1]), nh3a[::-1], bounds_error = False)

data = Dataset('mwr-%s-rad.nc' % args['case'], 'r')
tb1, tba = zeros(6), zeros(6)
for b in range(6):
  tb1[b] = float(data['b%dtoa1' % (b+1,)][0,0,0])
  tba[b] = float(data['b%dtoa1' % (b+1,)][0,0,3]) - tb1[b]
ptb1 = TPfunc(tb1)

# profiles
os.system('mkdir -p %s' % args['case'])

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
ax.plot(nh3b - nh3a, p4, 'C2--')
ax.plot(nh3d - lnpNH3func(log(p4)), p4, 'C2-')
ax.set_xlabel('Perturbed NH3 (g/kg)', fontsize = 15)

ax = axs[2]
ax.plot([0., 0.], [1.E3, 0.1], '--', color = '0.7')
ax.plot(T3 - T1, p1, 'C1--')
ax.plot(T4 - lnpTfunc(log(p4)), p4, 'C1-')
ax.errorbar(tba, ptb1, xerr = 0.01*(tb1-300.), fmt = 'o', color = 'C4', capsize = 5)
ax.set_xlabel('Perturbed T (K)', fontsize = 15)

ax = axs[3]
ax.plot(theta1, p1, '--', color = '0.7')
ax.plot(theta4v, p4, 'C1--')
ax.plot(theta4, p4, 'C1')
ax.scatter(175.*ones(6), ptb1, s = 50, ec = 'C4', fc = 'none')
ax.set_xlabel('Potential T (K)', fontsize = 15)
ax.set_xlim([160., 175.])

ax2 = ax.twiny()
ax2.plot(h2o4, p4, 'C0')
ax2.set_xlim([-5., 100.])
ax2.set_xlabel('mmr (g/kg)')


savefig('%s/profile_%s.png' % (args['case'],args['case']), bbox_inches = 'tight')
savefig('%s/profile_%s.eps' % (args['case'],args['case']), bbox_inches = 'tight')
