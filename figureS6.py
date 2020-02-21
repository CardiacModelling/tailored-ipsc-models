#!/usr/bin/env python
#
# Python code to generate Figure S6: Histogram of peak outward current
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
# Get caching options
#
force = '--force' in sys.argv

#
# Load experimental data
#
peaks = outward.load_outward_peaks(force)

#
# Plot histogram of peak current
#
pl.figure()
pl.xlabel('Peak outward current [A/F]')
pl.ylabel('Probability density')
#pl.xlim(125, 575)
#pl.ylim(0, 0.011)

for peak in sorted(peaks):
    print('  ' + str(peak))

bins = np.arange(0,70,5)

pl.hist(peaks,
    normed=True,
    bins=bins,
    label='Peak outward current',
    color='tab:green',
    edgecolor='tab:green',
    alpha=0.5,
)

methods.save('figureS6')

print('Figure S6 saved to ' + methods.FIGURES)

#pl.show()







