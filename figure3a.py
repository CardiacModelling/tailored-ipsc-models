#!/usr/bin/env python
#
# Python code to generate Figure 3: APs generated from tailored models
#
from __future__ import print_function
import os
import sys
import myokit
import methods
import methods.outward as outward
import numpy as np
import matplotlib
import matplotlib.pyplot as pl

#
# Fit all outward current experiments (or used cached result)
#
force = '--force' in sys.argv
table, cells = outward.fit_all(force=force)

#
# Run simulations with each cell (or use cached data)
#
aps = outward.simulate_tailored_aps(table, cells, force=force)

#
# Plot 2APs
#
#
# Order cells by APD at -65mV
#
tmin = 50 * 1e-3
vapd = -65 * 1e-3
order = np.array([aps['membrane.V', i] for i in xrange(len(outward.ORDER))])
offset = np.where(aps.time() >= tmin)[0][0]
order = order[:, offset:]
order = [np.where(v < vapd)[0][:1] for v in order]
order = [v[0] if v else 2000 for v in order]
order = np.argsort(np.argsort(order))

#
# Create colormap
#
cmap = matplotlib.cm.get_cmap('viridis')
norm = matplotlib.colors.Normalize(0, len(cells))

#
# Plot and save
#
pl.figure()
pl.xlabel('Time [ms]')
pl.ylabel('V [mV]')
pl.ylim(-85, 40)
for i, long_id in enumerate(outward.ORDER):
    short_id = table['short_id'][i]
    pl.plot(aps.time()*1e3, aps['membrane.V', i]*1e3,
        color=cmap(norm(order[i])),
        label=short_id)
pl.tight_layout()

methods.save('figure3a')

print('Figure saved to ' + methods.FIGURES)

#pl.show()
