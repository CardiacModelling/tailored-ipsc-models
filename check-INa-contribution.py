#!/usr/bin/env python
#
# Python code to generate Table 2: Fits to outward current for all cells
#
from __future__ import print_function
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import methods
import methods.outward as outward

#
# Fit all outward current experiments (or used cached result)
#
force = '--force' in sys.argv
table, cells = outward.fit_all(force=force)

# Time of interest, i.e. 0 - 20 ms
icapacitance = 0 # int(round(outward.CAPACITANCE / outward.DT))
iend = icapacitance + int(round( 20e-3 / outward.DT))
# Load simulated original INa
stime_org, sim_org = outward.simulate_outward_current(force)
INa = sim_org['ina.INa']

# Not scaled
sim_org = np.array([sim_org[current] for current in outward.CURRENTS], copy=True)
sim_org = np.sum(sim_org, axis=0)

# Compare INa
scaled_INa_peak = []
sim_peak = []
#exp_peak = []
for i in table['index']:
    i -= 1
    #print (np.min(table['scale.ina.INa'][i]*INa[:,icapacitance:iend], 1))
    scaled_INa_peak.append(np.min(table['scale.ina.INa'][i]*INa[:,icapacitance:iend], 1))
    cell = cells[i]
    exp = []
    sim = []
    for j in range(10):
        exp.append(cell[str(j)+'.exp'])
        sim.append(cell[str(j)+'.sim'])
    exp = np.array(exp)[:,icapacitance:iend]
    sim = np.array(sim)[:,icapacitance:iend]
    #print (np.min(exp, 1))
    #print (np.min(sim, 1))
    sim_peak.append(np.min(sim, 1))
    #exp_peak.append(np.min(exp, 1))

scaled_INa_peak = np.array(scaled_INa_peak)
sim_peak = np.array(sim_peak)
#exp_peak = np.array(exp_peak)

org_ratio = np.min(INa[:,icapacitance:iend], 1) / np.min(sim_org[:,icapacitance:iend], 1)

#plt.plot(cell['time'],INa[0])
#plt.plot(cell['time'],exp[0])
#plt.plot(cell['time'],sim[0])
#plt.savefig('ina.png')

#
# Show result
#
np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})
print('Printing ratio between peak of simulated scaled INa alone and simulated lumped current')
ratio_sim = scaled_INa_peak/sim_peak
print(ratio_sim)
print('Averaged ratio over all cells: ')
print(np.mean(ratio_sim, 0))
print('Original traces ratio: ')
print(org_ratio)
print('Each column is at different holding voltage')
n,m = ratio_sim.shape
print('Averaged ratio over everything: ', np.mean(ratio_sim.reshape(n*m,1)),' +/- ', np.std(ratio_sim.reshape(n*m,1)))
print('**'*20)
print('Averaged ratio over original: ', np.mean(org_ratio), ' +/- ', np.std(org_ratio))
#print('Printing ratio between peak of simulated scaled INa alone and experimental measure')
#print(scaled_INa_peak/exp_peak)
