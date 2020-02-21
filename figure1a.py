#!/usr/bin/env python
#
# Python code to generate Figure 1a: Sodium current experiments
#
# Will run simulations, read experimental data and store figures in the
# directory 'generated'
#
from __future__ import print_function
import os
import methods
import myokit
import numpy as np
import matplotlib.pyplot as pl

#
# Load model
#
model = myokit.load_model('paci-2013-ventricular.mmt')

#
# Solutions for Sodium protocol
#
def fix_concentration(variable, concentration):
    v = model.get(variable)
    if v.is_state():
        v.demote()
    v.set_rhs(concentration)

fix_concentration('sodium.Nai', 10)
fix_concentration('sodium.Nao', 150)
fix_concentration('potassium.Ki', 110)
fix_concentration('potassium.Ko', 4)
fix_concentration('calcium.Cao', 1.2)
fix_concentration('calcium.Cai', 0)

#
# Room temperature
#
model.get('phys.T').set_rhs(298)

#
# Create Sodium protocol
#

# Holding potential
vhold = -80e-3

# Step potentials
vsteps_mv = np.arange(-60, 60 + 10, 10)
vsteps = np.arange(-60, 60 + 2, 2)*1e-3 # Finer step for simulations

# Time between steps
tpre = 5000e-3

# Time at step
tstep = 20e-3

# Create myokit protocol, get duration
sodium_protocol = myokit.pacing.steptrain(vsteps, vhold, tpre, tstep)
duration = sodium_protocol.characteristic_time()

#
# Clamp membrane potential to protocol
#
vm = model.get('membrane.V')
vm.demote()
vm.set_rhs(-80e-3)
model.get('engine.pace').set_binding(None)
vm.set_binding('pace')

#
# Apply protocol
#
current = 'ina.INa'
#current = 'membrane.i_ion'
sim = myokit.Simulation(model, sodium_protocol)
data = sim.run(duration, log=['engine.time', current], log_interval=1e-4)

# Show raw simulation results
if False:
    pl.figure()
    pl.plot(data.time(), data[current], label=current)
    pl.legend()
    pl.show()

#
# Analyse
#
simulated = []
data = data.fold(tpre + tstep)
data = data.trim_left(tpre)
for i, v in enumerate(vsteps):
    simulated.append(np.min(data[current, i]))
simulated = np.array(simulated)

# Show steps
if False:
    pl.figure()
    for i in xrange(len(vsteps)):
        pl.plot(data.time(), data[current, i])
    pl.show()

#
# Load experimental data
#

#
# Cell capacitances
#
capacitance = []
filename = os.path.join('sodium-experiment', 'sodium-capacitance.csv')
with open(filename, 'r') as f:
    line = f.readline()
    capacitance = [float(x) for x in line.split(',')]
capacitance = np.array(capacitance)

#
# Cell currents
#
real = []
filename = os.path.join('sodium-experiment', 'sodium-current.csv')
with open(filename, 'r') as f:
    for line in f:
        real.append([float(x) for x in line.split(',')])
real = np.array(real)
real /= capacitance # Normalise!
real = real.T
ncells = len(real)

#
# Show experimental results
#
if False:
    pl.figure()
    for cell in real:
        pl.plot(vsteps_mv, cell)
    pl.show()

#
# Calculate median and percentiles
#
median, p25, p75 = methods.median_and_percentiles(real)
mean = np.mean(real, axis=0)

#
# Calculate scaling factor to match peak of mean to simulation
#
scale = np.min(mean) / np.min(simulated)
print('Scaling factor: ' + str(scale))

#
# Create plot
#
pl.figure()
pl.xlabel('V [mV]')
pl.ylabel('Peak I$_{Na}$ [A/F]')
pl.ylim(-235, 9)
methods.fill(vsteps_mv, p25, p75, 'red')
pl.plot(vsteps_mv, mean, 'o-', label='Experiment (mean, n=' + str(ncells) +')',
    color='#B22400')
pl.plot(vsteps*1e3, simulated, label='Simulation')
pl.plot(vsteps*1e3, simulated*scale, label='Fitted simulation')
pl.legend(loc='lower center', fontsize=9, ncol=3, columnspacing=1.3,
    handlelength=1.5, handletextpad=0.5)
pl.tight_layout()

# Store
methods.save('figure1a')
#pl.show()

print('Figure saved to ' + methods.FIGURES)

