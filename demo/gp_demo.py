#! /usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import argparse, os
from netCDF4 import Dataset
from pydrum.plots1d import PlotProfile, DrawPressureAxis
from scipy.interpolate import interp1d
from numpy.random import multivariate_normal
from pylab import *

def ExponentialQuadratic(x1, x2, s1, s2, l):
  n, m = len(x1), len(x2)
  X = zeros((n,m))
  for i in range(n):
    for j in range(m):
      X[i,j] = s1[i]*s2[j]*exp(-(x1[i] - x2[j])**2/(2.*l*l))
  return X

parser = argparse.ArgumentParser()
parser.add_argument('--Tstd',
    default = '4.',
    help = 'standard deviation of temperature'
    )
parser.add_argument('--Length',
    default = '60.',
    help = 'correlation length'
    )
args = vars(parser.parse_args())

case = 'ma166K'

# pressure divides
pdiv = [50, 15, 5, 2, 0.6]

# standard deviation and covariance length
stemp, ltemp = float(args['Tstd']), float(args['Length'])

# reference temperature at 1 bar
T0 = 166.

# locate height given pressure divides
data = Dataset('mwr-%s-main.nc' % case, 'r')
x1 = data['x1'][:]/1.E3
pres = mean(data['press'][0,:,:,:], axis = (0,1))/1.E5
temp = mean(data['temp'][0,:,:,:], axis = (0,1))
lnpfunc = interp1d(log(pres[::-1]), x1[::-1])
Tfunc = interp1d(x1[::-1], temp[::-1])

zdiv, zlev = zeros(7), zeros(6)
zdiv[0] = x1.min()
zdiv[1:-1] = lnpfunc(log(pdiv))
zdiv[-1] = x1.max()
#print(zdiv)

# sample height
for i in range(len(zlev)):
  r = rand()
  zlev[i] = zdiv[i]*(1.-r) + zdiv[i+1]*r
#print(zlev)

# find temperature given the sample height
Tlev = Tfunc(zlev)
#print(Tlev)

# sample pertubation temperature
slev = stemp*Tlev/T0
Clev = ExponentialQuadratic(zlev, zlev, slev, slev, ltemp)
dt = multivariate_normal(mean = zeros(6), cov = Clev, size = 5).T

# calculate mean perturbed profile
s1 = stemp*temp/T0
cov = ExponentialQuadratic(x1, zlev, s1, slev, ltemp)
#print(cov)
ptemp = dot(cov,linalg.solve(Clev, dt))
#print(ptemp)

os.system('mkdir -p gp_demo')

figure(1, figsize = (8, 10))
ax = axes()
ylims = [min(x1), max(x1)]

ax.plot(ptemp, x1)
#ax.set_xlim([160., 200.])
ax.set_xlabel('Perturbed temperature (K)', fontsize = 15)
ax.set_ylabel('Height (km)', fontsize = 15)
ax.set_ylim(ylims)

paxis = [1000., 300., 100., 50., 20., 10., 5., 2., 1., 0.5, 0.1]
DrawPressureAxis(ax, x1, paxis, pres, ylims)
#ax2.set_ylim(ylims)

savefig('gp_demo/pt_%s_s=%d_l=%d.png' % (case,stemp,ltemp), bbox_inches = 'tight')
savefig('gp_demo/pt_%s_s=%d_l=%d.eps' % (case,stemp,ltemp), bbox_inches = 'tight')

# rectify profile
dT = temp[0] - temp[1]
print(dT)
Trec = temp.reshape((len(x1),1)) + ptemp
print(Trec.shape)

for j in range(5):
  dTmax = 5.
  for i in range(1, len(x1)):
    if Trec[i-1,j] - Trec[i,j] > dT:  # potentially unstable
      if dTmax < Trec[i-1,j] - Trec[i,j] - dT: # not enough water
        Trec[i,j] = Trec[i-1,j] - dT - dTmax
        dTmax = 0
      else: # has enough water
        dTmax -= Trec[i-1,j] - Trec[i,j] - dT

# plot rectified profile
figure(2, figsize = (8, 10))
ax = axes()
ylims = [min(x1), max(x1)]

ax.plot(Trec - temp.reshape((len(x1),1)), x1)
#ax.set_xlim([160., 200.])
ax.set_xlabel('Rectified pertubation temperature (K)', fontsize = 15)
ax.set_ylabel('Height (km)', fontsize = 15)
ax.set_ylim(ylims)

paxis = [1000., 300., 100., 50., 20., 10., 5., 2., 1., 0.5, 0.1]
DrawPressureAxis(ax, x1, paxis, pres, ylims)
#ax2.set_ylim(ylims)

savefig('gp_demo/rpt_%s_s=%d_l=%d.png' % (case,stemp,ltemp), bbox_inches = 'tight')
savefig('gp_demo/rpt_%s_s=%d_l=%d.eps' % (case,stemp,ltemp), bbox_inches = 'tight')

close(2)

# potential temperature
R_ov_cp = log(temp[1]/temp[0])/log(pres[1]/pres[0])
print(R_ov_cp)

figure(3, figsize = (8, 10))
ax = axes()
ylims = [min(x1), max(x1)]
ax.plot(Trec*pow(1./pres, R_ov_cp).reshape((len(x1),1)), x1)
ax.set_xlabel('Potential temperature (K)', fontsize = 15)
ax.set_ylabel('Height (km)', fontsize = 15)
ax.set_ylim(ylims)
ax.set_xlim([160, 180])

paxis = [1000., 300., 100., 50., 20., 10., 5., 2., 1., 0.5, 0.1]
DrawPressureAxis(ax, x1, paxis, pres, ylims)

savefig('gp_demo/theta_%s_s=%d_l=%d.png' % (case,stemp,ltemp), bbox_inches = 'tight')
savefig('gp_demo/theta_%s_s=%d_l=%d.eps' % (case,stemp,ltemp), bbox_inches = 'tight')
