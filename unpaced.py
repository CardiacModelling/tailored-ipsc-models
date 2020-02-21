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
#matplotlib.use('Agg')
import matplotlib.pyplot as pl

# Don't cache
force = "--force" in sys.argv

#
# Fit all outward current experiments (or used cached result)
#
table, cells = outward.fit_all()

#
# Run simulations with each cell (or use cached data)
#
sim = outward.simulate_tailored_aps(table, cells, force=force, beats=100, cl=1,
 stimulate=False)

# Normalise and calculate apd
n = len(outward.ORDER)
sim_time = sim.time() * 1e3 - 50
sim_aps = []
sim_apds = []
for i, long_id in enumerate(outward.ORDER):
    filename = os.path.join(apd.AP_FIGURES, 'sim-base-' + long_id + '.png')
    #sim_aps.append(apd.normalise(sim_time, sim['membrane.V', i],
    #    filename=filename))
    sim_aps.append(sim['membrane.V', i])
    filename = os.path.join(apd.APD_FIGURES, 'sim-base-' + long_id + '.png')
    #sim_apds.append(apd.calculate_apd(sim_time, sim_aps[-1],
    #    filename=filename)[0])
sim_aps = np.array(sim_aps)
#sim_apds = np.array(sim_apds)


#
# Plot APs
#
pl.figure(figsize=(15,5))
pl.xlabel('Time [ms]')
pl.ylabel('V (normalised)')
#pl.xlim(-50, 550)
pl.xlim(-50,20000)
pl.grid(True)

# Plot tailored model APs
n = len(sim_aps)
for ap in sim_aps:
    try:
        whereap = sim_time[ap>-0.04][0]
    except:
        whereap = 0
    try:
        tempidx = np.arange(len(ap))[ap>-0.04]
        whereend = [tempidx[i+1] for i in range(len(tempidx)-1) if tempidx[i+1]-tempidx[i]>1000][1] - 1000
    except:
        whereend = -1000
    pl.plot((sim_time-whereap)[:whereend], ap[:whereend], alpha=0.5, lw=2)

#pl.legend(loc='upper right')

#pl.show()

methods.save('unpaced')





