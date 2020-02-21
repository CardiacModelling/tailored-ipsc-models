#!/usr/bin/env python
#
# Python code to generate Figure S7: APs and corresponding current contributions generated from tailored models
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
aps = outward.simulate_tailored_aps(table, cells, store_currents=True,
    force=force)

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
# Create figure
#
fig, axes = pl.subplots(nrows=6, ncols=1, sharex=True,figsize=(5,7.5))
fig.text(0.5, 0.05, r"Time [ms]", ha='center', va='center')

#
# Plot and save
#
axes[0].set_ylabel('V [mV]')
#pl.ylim(-85, 40)
mint=np.inf
maxt=-np.inf
for i, long_id in enumerate(outward.ORDER):
    short_id = table['short_id'][i]
    axes[0].plot(aps.time()*1e3, aps['membrane.V', i]*1e3,
        color=cmap(norm(order[i])),
        label=short_id)
    mint = min(np.min(aps['membrane.V',i]*1e3),mint)
    maxt = max(np.max(aps['membrane.V',i]*1e3),maxt)
#axes[0].set_yticklabels(np.linspace(int(mint),int(maxt),3))
#axes[0].set_ylim([mint-0.075*mint,maxt+0.075*maxt])
print(np.linspace(int(mint),int(maxt),3))

axes[1].set_ylabel('I$_{Na}$ [A/F]')
#pl.ylim(-85, 40)
for i, long_id in enumerate(outward.ORDER):
    short_id = table['short_id'][i]
    axes[1].plot(aps.time()*1e3, aps['ina.INa', i],
        color=cmap(norm(order[i])),
        label=short_id)

axes[2].set_ylabel('I$_{CaL}$ [A/F]')
#pl.ylim(-85, 40)
for i, long_id in enumerate(outward.ORDER):
    short_id = table['short_id'][i]
    axes[2].plot(aps.time()*1e3, aps['ical.ICaL', i],
        color=cmap(norm(order[i])),
        label=short_id)

axes[3].set_ylabel('I$_{Ks}$ [A/F]')
#pl.ylim(-85, 40)
for i, long_id in enumerate(outward.ORDER):
    short_id = table['short_id'][i]
    axes[3].plot(aps.time()*1e3, aps['iks.IKs', i],
        color=cmap(norm(order[i])),
        label=short_id)

axes[4].set_ylabel('I$_{NaCa}$ [A/F]')
#pl.ylim(-85, 40)
for i, long_id in enumerate(outward.ORDER):
    short_id = table['short_id'][i]
    axes[4].plot(aps.time()*1e3, aps['inaca.INaCa', i],
        color=cmap(norm(order[i])),
        label=short_id)
#pl.tight_layout()

axes[5].set_ylabel('I$_{Kr}$ [A/F]')
#pl.ylim(-85, 40)
for i, long_id in enumerate(outward.ORDER):
    short_id = table['short_id'][i]
    axes[5].plot(aps.time()*1e3, aps['ikr.IKr', i],
        color=cmap(norm(order[i])),
        label=short_id)

fig.subplots_adjust(hspace=0)
#plt.setp([a.get_xticklabels() for a in fig.axes[:-1]], visible=False)
#pl.savefig('iall.png',bbox_inches='tight')

for a in axes:
    pl.setp(a.get_yticklabels()[0], visible=False)    
    pl.setp(a.get_yticklabels()[-1], visible=False)
pl.setp(axes[5].get_yticklabels()[-2], visible=False)

methods.save('figureS7')

print('Figure S7 saved to ' + methods.FIGURES)

#pl.show()
