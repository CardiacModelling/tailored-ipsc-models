#!/usr/bin/env python
#
# Python code to generate Figure S5: Outward current fit to three cells.
#
from __future__ import print_function
import sys
import methods
import methods.outward as outward
import matplotlib.pyplot as pl

#
# Fit all outward current experiments (or used cached result)
#
force = '--force' in sys.argv
table, cells = outward.fit_all(force)

#
# Plot
#
for i, long_id in enumerate(outward.ORDER):
    print('Plotting fit for cell ' + long_id)
    result = outward.select(table, cells, i)
    outward.plot_experiment(*result)
    methods.save('all-outward-fits/' + long_id + '_exp')
    outward.plot_simulation(*result)
    methods.save('all-outward-fits/' + long_id + '_sim')
    pl.close('all')
    
print('Figures saved to ' + methods.FIGURES + ' all-outward-fits')
