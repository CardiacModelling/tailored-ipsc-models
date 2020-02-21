#!/usr/bin/env python
#
# Python code to generate Figure 4: APs/APDs of simulations compared to optical
# mapping results
#
from __future__ import print_function
import os
import sys
import methods
import methods.outward as outward
import methods.apd as apd
import numpy as np
import matplotlib
import matplotlib.pyplot as pl

#
# Get debug/caching options
#
force = '--force' in sys.argv
if '--apfigures' in sys.argv:
    apd.DEBUG_FIGURES = True

#
# Point at which the APD is measured
#
if '--apd50' in sys.argv:
    PT = 0.5
else:
    PT = 0.1    

#
# Median or mean
#
MEDIAN = True

#
# Load experimental data
#
exp_time, exp_aps = apd.load_base_aps(force=force)
exp_time -= 50
exp_apds = apd.calculate_base_apds(exp_time, exp_aps, PT, force=force)

#
# Fit all outward current experiments (or used cached result)
#
table, cells = outward.fit_all(force=force)

#
# Run simulations with each cell (or use cached data)
#
sim = outward.simulate_tailored_aps(table, cells, force=force)

# Normalise and calculate apd
n = len(outward.ORDER)
sim_time = sim.time() * 1e3 - 50
sim_aps = []
sim_apds = []
for i, long_id in enumerate(outward.ORDER):
    filename = os.path.join(apd.AP_FIGURES, 'sim-base-' + long_id + '.png')
    sim_aps.append(apd.normalise(sim_time, sim['membrane.V', i],
        filename=filename))
    filename = os.path.join(apd.APD_FIGURES, 'sim-base-' + long_id + '.png')
    sim_apds.append(apd.calculate_apd(sim_time, sim_aps[-1],
        filename=filename)[0])
sim_aps = np.array(sim_aps)
sim_apds = np.array(sim_apds)

# Calculate apds
sim_apds = apd.calculate_tailored_apds(sim_time, sim_aps, PT)

#
# Run simulation with original model (or use cached data)
#
org = outward.simulate_original_aps(force=force)
org_time = org.time()*1000 - 50
org_ap = org['membrane.V']
org_ap = apd.normalise(org_time, org_ap)
org_apd, alternans = apd.calculate_apd(org_time, org_ap, PT)
del(org)

#
# Figure 4a: Plot median and quartiles of APs
#
pl.figure()
pl.xlabel('Time [ms]')
pl.ylabel('V (normalised)')
pl.xlim(-50, 550)
pl.grid(True)

# Plot experimental data
if MEDIAN:
    mid, lo, hi = methods.median_and_percentiles(exp_aps)
else:
    mid, lo, hi = methods.mean_and_std(exp_aps)

methods.fill(exp_time, lo, hi, color='black')
pl.plot(exp_time, mid, color='black', label='Optical data ('
    + ('median' if MEDIAN else 'mean') + ', n=94)')
#for ap in exp_aps:
#    pl.plot(exp_time, ap, color='black', alpha=0.1, lw=1)


# Plot original model AP
pl.plot(org_time, org_ap, '--', label='Original model')

# Plot tailored model APs
n = len(sim_aps)
if MEDIAN:
    mid, lo, hi = methods.median_and_percentiles(sim_aps)
else:
    mid, lo, hi = methods.mean_and_std(sim_aps)
methods.fill(sim_time, lo, hi, color='orange')
pl.plot(sim_time, mid, lw=3, label='Fitted models ('
    + ('median' if MEDIAN else 'mean') + ', n='+str(n)+')')
#for ap in sim_aps:
#    pl.plot(sim_time, ap, color='tab:orange', alpha=0.25, lw=2)

# Show APD line
pl.axhline(PT, ls='--', color='gray', label=str(PT))

# Add legend
pl.legend(loc='upper right')

# Store
methods.save('figure4a' if PT == 0.1 else 'figure104a')

#
# Figure 4b: Plot histogram of normalised APDs
#
pl.figure()
pl.xlabel('APD$_{' + ('90' if PT == 0.1 else '50') + '}$ [ms]')
pl.ylabel('Probability density')
if PT == 0.1:
    pl.xlim(150, 500)
elif PT == 0.5:
    pl.xlim(100, 375)

bins = np.arange(0,600,25)

pl.hist(exp_apds,
    normed=True,
    bins=bins,
    label='Optical data',
    color='k',
    edgecolor='k',
    alpha=0.25,
)

pl.axvline(org_apd, lw=3, label='Original model')

pl.hist(sim_apds,
    normed=True,
    bins=bins,
    label='Tailored models',
    color='tab:orange',
    edgecolor='tab:orange',
    alpha=0.5,
)

pl.legend(loc='upper right')

methods.save('figure4b' if PT == 0.1 else 'figure104b')

print('Figures saved to ' + methods.FIGURES)

#pl.show()







