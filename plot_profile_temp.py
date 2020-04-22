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
T1 = data['temp'][0,0,0,:]
T2 = data['temp'][0,0,1,:] - T1
T3 = data['temp'][0,0,2,:] - T1
T4 = data['temp'][0,0,3,:]

theta1 = data['theta'][0,0,0,:]
theta4 = data['theta'][0,0,3,:]
theta4v = data['thetav'][0,0,3,:]

h2o1 = data['vapor1'][0,0,0,:]*1.E3
h2o4 = data['vapor1'][0,0,3,:]*1.E3

x1 = data['x1'][:]/1.E3
pres = mean(data['press'][0,:,:,:], axis = (0,1))/1.E5
Tfunc = interp1d(T1[::-1], x1[::-1], bounds_error = False)

data = Dataset('mwr-%s-rad.nc' % args['case'], 'r')
tb1, tba = zeros(6), zeros(6)
for b in range(6):
  tb1[b]  = float(data['b%dtoa1' % (b+1,)][0,0,0])
  tba[b] = float(data['b%dtoa1' % (b+1,)][0,0,3]) - tb1[b]
ztb1 = Tfunc(tb1)

#ylims = [-800., 50.]
ylims = [min(x1), max(x1)]

# profiles
os.system('mkdir -p %s' % args['case'])

fig, axs = subplots(1, 4, figsize = (14,8), sharey = True)
subplots_adjust(wspace = 0.08, hspace = 0.08)
ax = axs[0]
ax.plot(T1, x1, 'C1')
ax.scatter(tb1, ztb1, s = 50, ec = 'C4', fc = 'none')
ax.set_ylabel('Height (km)', fontsize = 15)
ax.set_xlabel('T (K)', fontsize = 15)
ax.set_ylim(ylims)

ax2 = ax.twiny()
ax2.plot(h2o1, x1, 'C0')
ax2.set_xscale('log')
ax2.set_xlim([0.1, 100.])
ax2.set_xlabel('H2O mmr (g/kg)')

ax = axs[1]
ax.plot(T2, x1, 'C1')
ax.plot([0., 0.], ylims, 'k--', alpha = 0.5)
ax.set_xlabel('Perturbation T (K)', fontsize = 15)
#ax.set_xlim([-20., 20.])

ax = axs[2]
ax.plot(T2, x1, 'C1--')
ax.plot(T3, x1, 'C1')
ax.scatter(zeros(6), ztb1, s = 50, ec = 'C4', fc = 'none')
ax.errorbar(tba, ztb1, xerr = 0.01*(tb1-300.)+0.2, fmt = 'o', color = 'C4', capsize = 5)
ax.plot([0., 0.], ylims, 'k--', alpha = 0.5)
ax.set_xlabel('Rectified T (K)', fontsize = 15)
#ax.set_xlim([-20., 20.])

ax = axs[3]
ax.plot(theta4v, x1, 'C1--')
ax.plot(theta4, x1, 'C1')
ax.scatter(180.*ones(6), ztb1, s = 50, ec = 'C4', fc = 'none')
ax.set_xlabel('Potential T (K)', fontsize = 15)
ax.set_xlim([150., 180.])

ax2 = ax.twiny()
ax2.plot(h2o4, x1, 'C0')
#ax2.set_xscale('log')
ax2.set_xlim([-5., 100.])
ax2.set_xlabel('H2O mmr (g/kg)')

paxis = [1E3, 300., 100., 50., 20., 10., 5., 2., 1., 0.5, 0.1]
DrawPressureAxis(axs[3], x1, paxis, pres, ylims)

savefig('%s/temp_%s.png' % (args['case'],args['case']), bbox_inches = 'tight')
savefig('%s/temp_%s.eps' % (args['case'],args['case']), bbox_inches = 'tight')
