#!/usr/bin/env python
#
# Python code to generate Figure 2: Outward current fit to cell 'b3'
#
from __future__ import print_function
import sys
import methods
import methods.outward as outward

#
# Fit all outward current experiments (or used cached result)
#
force = '--force' in sys.argv
table, cells = outward.fit_all(force)

#
# Find cell by name
#
cell = outward.find('b3')

#
# Fit to data
#
#
# Plot
#
result = outward.select(table, cells, cell)
outward.plot_experiment(*result)
methods.save('figure2a')
outward.plot_simulation(*result)
methods.save('figure2b')

print('Figures 2a and 2b saved to ' + methods.FIGURES)
