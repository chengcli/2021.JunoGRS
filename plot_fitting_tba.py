#! /usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
from pylab import *
from pyfits.gpmcmc import *

# PJ1345689
lats, targets, errs, vala, lnps = read_all_fits('PJ1345689V3', 'mwr')

# PJ07
lats, targets, errs, val7, lnps = read_all_fits('PJ07BOOSTV3', 'mwr')

# brightness temperature anomaly
freq = [0.6, 1.25, 2.6, 5.2, 10., 22.]
X, Y = meshgrid(lats, freq)
tba = val7[:,::3] - vala[:,::3]

figure(3, figsize = (8, 8))
ax = axes()
latmin, latmid, latmax = -26, -20, -14
ax.plot([latmid, latmid], [freq[0], freq[-1]], 'C1--', linewidth = 2)
ax.plot([-16, -16], [freq[0], freq[-1]], 'C1-', linewidth = 2)
ax.plot([-24, -24], [freq[0], freq[-1]], 'C1-', linewidth = 2)

h1 = ax.contour(X, Y, tba.T, linspace(-11, -1, 5), colors = 'b', linestyles = '-')
h2 = ax.contour(X, Y, tba.T, linspace(1, 26, 11), colors = 'r')

ax.set_xlim([latmin, latmax])
ax.set_xticks(arange(latmin, latmax + 2, 2))
ax.set_yscale('log')
ax.set_yticks(freq)
ax.set_yticklabels(freq)
ax.set_ylabel('Frequency (GHz)', fontsize = 18)
clabel(h1, fontsize = 12, inline = 1, fmt = '%.1f')
clabel(h2, fontsize = 12, inline = 1, fmt = '%.1f')

savefig('GRS_fitted_Tba_BOOSTV3.png', bbox_inches = 'tight')
#savefig('GRS_fitted_Tba_V1.pdf', bbox_inches = 'tight')
