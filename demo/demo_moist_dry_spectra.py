#! /usr/bin/env python3
from pylab import *

freq = array([0.6,1.25,2.6,5.2,10.,22.])

ma = genfromtxt('ma166k.obs', skip_header = 2, max_rows = 6)
da = genfromtxt('da166k.obs', skip_header = 2, max_rows = 6)

figure(1, figsize = (8, 8))
ax = axes()
ax.errorbar(freq, ma[:,0], yerr = 0.02*(ma[:,0]-300.), fmt = 'o', color = 'C4', capsize = 5)
ax.plot(freq, da[:,0], 'k--')
ax.set_xlabel('Frequencies (GHz)')
ax.set_ylabel('Tb nadir (K)')
ax.set_xscale('log')

savefig('demo_moist_dry_spectra.png', bbox_inches = 'tight')
