#!/usr/bin/env python
#
# Python code to generate Figure S4: Outward current in original model
#
from __future__ import print_function
import numpy as np
import methods
import methods.outward as outward
import matplotlib
import matplotlib.pyplot as pl
import myokit

#
# Simulate original currents
#

# Load model
model = myokit.load_model('paci-2013-ventricular.mmt')

# Set solutions for outward current protocol
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

# Set room temperature
model.get('phys.T').set_rhs(298)

# Create outward current protocol
vhold = -80e-3
vsteps = outward.VSTEPS_mV * 1e-3
tpre = 10 # 10 s
tstep = 0.5 # 500ms
tpost = 0.010
sodium_protocol = myokit.pacing.steptrain(vsteps, vhold, tpre, tstep, tpost)
duration = sodium_protocol.characteristic_time()

# Clamp membrane potential to protocol
vm = model.get('membrane.V')
vm.demote()
vm.set_rhs(-80e-3)
model.get('engine.pace').set_binding(None)
vm.set_binding('pace')

# Apply protocol
sim = myokit.Simulation(model, sodium_protocol)
log_vars = ['engine.time'] + outward.CURRENTS
data = sim.run(duration, log=log_vars, log_interval=outward.DT)
data = data.npview()

# Create overlapping step data
data = data.fold(tpre + tstep + tpost)
data = data.trim_left(tpre - tpost, adjust=True)

# Get time, outward current
time = data.time() * 1e3
total = []
for i, v in enumerate(outward.VSTEPS_mV):
    total_i = np.zeros(time.shape)
    for current in outward.CURRENTS:
        total_i += data[current, i]
    total.append(total_i)

#
# Plot
#
cmap = matplotlib.cm.get_cmap('viridis')
norm = matplotlib.colors.Normalize(0, len(total))
pl.figure()
pl.xlabel('Time [ms]')
pl.ylabel('Outward current [A/F]')
pl.ylim(-6, 2)

for i, current in enumerate(total):
    vm = outward.VSTEPS_mV[i]
    pl.plot(time, current, color=cmap(norm(i)), label=str(vm)+' mV')
pl.legend(loc='lower right', ncol=5, columnspacing=1.2, fontsize=9,
        handlelength=1.5, handletextpad=0.5)
methods.save('figureS4')

print('Figure S4 saved to ' + methods.FIGURES)

#pl.show()
