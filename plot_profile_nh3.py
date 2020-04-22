#! /usr/bin/env python3
#import matplotlib
#matplotlib.use('Agg')
import argparse, h5py
from netCDF4 import Dataset
from pydrum.plots1d import PlotProfile, DrawPressureAxis
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

h2o1 = data['vapor1'][0,0,0,:]*1.E3
h2o4 = data['vapor1'][0,0,3,:]*1.E3

nh3a = data['vapor2'][0,0,0,:]*1.E3
nh3b = data['vapor2'][0,0,1,:]*1.E3 - nh3a
nh3c = data['vapor2'][0,0,2,:]*1.E3 - nh3a
nh3d = data['vapor2'][0,0,3,:]*1.E3

x1 = data['x1'][:]/1.E3
pres = mean(data['press'][0,:,:,:], axis = (0,1))/1.E5

ylims = [-800., 50.]

# profiles
os.system('mkdir -p %s' % args['case'])

fig, axs = subplots(1, 4, figsize = (14,8), sharey = True)
subplots_adjust(wspace = 0.08, hspace = 0.08)
ax = axs[0]
ax.plot(nh3a, x1, 'C2')
#ax.plot(h2o1, x1, 'C0')
ax.set_ylabel('Height (km)', fontsize = 15)
ax.set_ylim(ylims)
ax.set_xlabel('Mass mixing ratio (g/kg)')
ax.set_xscale('log')
ax.set_xlim([0.1, 100.])

ax2 = ax.twiny()
ax2.plot(T1, x1, 'C1', lw = 2)
ax2.set_xlabel('T (K)', fontsize = 15)

ax = axs[1]
ax.plot(nh3b, x1, 'C2')
ax.plot([0., 0.], ylims, 'k--', alpha = 0.5)
ax.set_xlabel('Perturbation (g/kg)', fontsize = 15)

ax = axs[2]
ax.plot(nh3b, x1, 'C2--')
ax.plot(nh3c, x1, 'C2')
ax.plot([0., 0.], ylims, 'k--', alpha = 0.5)
ax.set_xlabel('Rectified (g/kg)')

ax = axs[3]
#ax.plot(h2o4, x1, 'C1')
ax.plot(nh3d, x1, 'C2')
#ax.set_xlabel('(K)', fontsize = 15)
#ax.set_xlim([150., 180.])
ax.set_xlabel('Mass mixing ratio (g/kg)')
#ax.set_xlim([0., 100.])

ax2 = ax.twiny()
ax2.plot(theta4, x1, 'C1', lw = 2)
ax2.set_xlim([150., 180.])
ax2.set_xlabel('Potential T (K)', fontsize = 15)

paxis = [1E3, 300., 100., 50., 20., 10., 5., 2., 1., 0.5, 0.1]
DrawPressureAxis(axs[3], x1, paxis, pres, ylims)

savefig('%s/nh3_%s.png' % (args['case'],args['case']), bbox_inches = 'tight')
savefig('%s/nh3_%s.eps' % (args['case'],args['case']), bbox_inches = 'tight')

# contribution function
#figure(2, figsize = (6, 8))
#ax = axes()

#print(temp)
#PlotProfile(ax, temp, pres, color = 'k')
#ax.set_xlabel('Temperature (K)', fontsize = 15)
#ax.set_ylabel('Pressure (bar)', fontsize = 15)
#ax.set_ylim([1.E4, 0.1])
#ax.set_yscale('log')

#ax2 = ax.twiny()
#ff = h5py.File('contribution_function.h5', 'r')
#cf1 = ff['bellotti']['contribution_function_channel_1']
#cf2 = ff['bellotti']['contribution_function_channel_2']
#cf3 = ff['bellotti']['contribution_function_channel_3']
#cf4 = ff['bellotti']['contribution_function_channel_4']
#cf5 = ff['bellotti']['contribution_function_channel_5']
#cf6 = ff['bellotti']['contribution_function_channel_6']
#pcon = ff['bellotti']['pressure']

#ax2.plot(cf1, pcon, 'k--', color = '0.7', linewidth = 2)
#ax2.plot(cf2, pcon, 'k--', color = '0.7', linewidth = 2)
#ax2.plot(cf3, pcon, 'k--', color = '0.7', linewidth = 2)
#ax2.plot(cf4, pcon, 'k--', color = '0.7', linewidth = 2)
#ax2.plot(cf5, pcon, 'k--', color = '0.7', linewidth = 2)
#ax2.plot(cf6, pcon, 'k--', color = '0.7', linewidth = 2)
#ax2.set_xlabel('Contribution function', fontsize = 15)

#savefig('%s/cf_%s.png' % (args['case'],args['case']), bbox_inches = 'tight')
#savefig('%s/cf_%s.eps' % (args['case'],args['case']), bbox_inches = 'tight')
#close(2)
