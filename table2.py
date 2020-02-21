#!/usr/bin/env python
#
# Python code to generate Table 2: Fits to outward current for all cells
#
from __future__ import print_function
import sys
import methods
import methods.outward as outward

#
# Fit all outward current experiments (or used cached result)
#
force = '--force' in sys.argv
table, cells = outward.fit_all(force=force)

#
# Show result
#
#outward.print_table(table)
outward.tex_table(table)
