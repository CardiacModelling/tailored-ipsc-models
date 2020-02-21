#!/usr/bin/env python
#
# Python code to generate Figure S2
#
from __future__ import print_function
import sys
import myokit
import methods
import methods.outward as outward
import numpy as np
import matplotlib
import matplotlib.pyplot as pl

CURRENTS = [
    'ina.INa',
    'ical.ICaL',
    'ik1.IK1',
    'ikr.IKr',
    'iks.IKs',
    'ito.Ito',
    'if.If',
    'inaca.INaCa',
    'inak.INaK',
    'ibna.IbNa',
    ]
    
LABELS = {
    'ina.INa' : 'INa',
    'ical.ICaL' : 'ICaL',
    'ik1.IK1' : 'IK1',
    'ikr.IKr' : 'IKr',
    'iks.IKs' : 'IKs',
    'ito.Ito' : 'Ito',
    'if.If' : 'If',
    'inaca.INaCa' : 'INaCa',
    'inak.INaK' : 'INaK',
    'ibna.IbNa' : 'INa,b',
}

# Load model and protocol
model, protocol, x = myokit.load('paci-2013-ventricular.mmt')

# Pre-pace at 1Hz
s = myokit.Simulation(model, protocol)
s.pre(100)
model.set_state(s.state())
del(s, protocol)

# Solutions for outward current protocol
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

# Room temperature
model.get('phys.T').set_rhs(298)

# Create outward current protocol
vhold = -80e-3
vsteps = outward.VSTEPS_mV * 1e-3
tpre = 10 # 10 s
tstep = 0.5 # 500ms
tpost = 0.01
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
log = sim.run(duration, log=['engine.time'] + CURRENTS,
    log_interval=outward.DT).npview()

# Create overlapping step data
log = log.fold(tpre + tstep + tpost)
log = log.trim_left(tpre - 0.01, adjust=True)

# Create colormap
cmap = matplotlib.cm.get_cmap('viridis')
norm = matplotlib.colors.Normalize(min(vsteps), max(vsteps))

# Create plot
pl.figure(figsize=(methods.FIGSIZE[0]*2.1, methods.FIGSIZE[1]*1.3))
n = len(CURRENTS)
n2 = int(np.ceil(n * 0.5))
for i, current in enumerate(CURRENTS):
    pl.subplot(2, n2, i+1)
    if i >= n * 0.5:
        pl.xlabel('Time [ms]')
    if i % n2 == 0:
        pl.ylabel('I [A/F]')
    pl.title(LABELS[current])
    for i, v in enumerate(vsteps):
        pl.plot(log.time()*1e3, log[current, i], color=cmap(norm(v)), label=str(v*1e3)+' mV')
pl.tight_layout()

pl.legend(bbox_to_anchor=(1.05, 0), loc='lower left', borderaxespad=0.)
#plt.figlegend( lines, labels, loc = 'lower center', ncol=5, labelspacing=0. )


methods.save('figureS2')

print('Figures saved to ' + methods.FIGURES)

#pl.show()

