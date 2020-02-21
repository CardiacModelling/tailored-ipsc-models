#!/usr/bin/env python
#
# Plots two successive APs of the current-clamped CID cells.
#
from __future__ import print_function
import os
import scipy.io
import numpy as np
import matplotlib
import matplotlib.pyplot as pl
import methods

#
# Two traces from several cells
#
filenames = [
    '170120-cell-5-middle.csv',
    '170123-cell-2-middle.csv',
    '170124-cell-2-middle.csv',
    '170130-cell-6-middle.csv',
    '170131-cell-3-middle.csv',
    '170120-cell-3-middle.csv',
    '170130-cell-3-middle.csv',
    ]

#
# Load data
#
cells = []
for filename in filenames:
    filename = os.path.join('current-clamp', filename)
    ts = []
    vs = []
    with open(filename, 'r') as f:
        for line in f:
            t, v = line.split(',')
            ts.append(float(t))
            vs.append(float(v))
    cells.append((np.array(ts), np.array(vs)))

#
# Order cells by APD at -20mV
#
tmin = 100
vapd = -20
i = np.where(cells[0][0] > tmin)[0][0]
cells.sort(key = lambda cell: np.where(cell[1][i:] < vapd)[0][:1])

#
# Create colormap
#
cmap = matplotlib.cm.get_cmap('inferno')
norm = matplotlib.colors.Normalize(0, len(cells))

#
# Plot
#
pl.figure()
pl.ylim(-85, 40)
for k, (ts, vs) in enumerate(cells):   
    pl.plot(ts, vs, color=cmap(norm(k)), lw=0.5)
pl.xlabel('Time [ms]')
pl.ylabel('V [mV]')
pl.tight_layout()
methods.save('figure3b')

print('Figure saved to ' + methods.FIGURES)

#pl.show()
