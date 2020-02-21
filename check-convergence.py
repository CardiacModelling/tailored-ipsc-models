#!/usr/bin/env python
#
# Python code to generate Table 2: Fits to outward current for all cells
#
from __future__ import print_function
import numpy as np
import sys
import methods
import methods.outward as outward

#
# Fit all outward current experiments 
#
force = True #'--force' in sys.argv
ref_vector = []
fit_tol = 1e-5
cutoff = 1e-10
compare_fields = [
        'scale.iks.IKs',
        'scale.inaca.INaCa',
        'scale.ikr.IKr',
        'scale.ik1.IK1',
        'scale.if.If',
        'scale.ito.Ito'
        ]
max_err = []
for i in range(10):
    table, cells = outward.fit_all(force=force)
    if i == 0:
        for field in compare_fields:
            ref_vector.append(table[field][:])
        ref_vector = np.array(ref_vector)
        ref_vector[ref_vector<cutoff] = 0
        print('*'*28)
        print('Testing convergence within error less than ',fit_tol)
    else:
        fitted_vector = []
        for field in compare_fields:
            fitted_vector.append(table[field][:])
        fitted_vector = np.array(fitted_vector)
        fitted_vector[fitted_vector<cutoff] = 0
        #fitted_vector = np.abs(fitted_vector - ref_vector) # abs diff
        norm = np.linalg.norm(fitted_vector - ref_vector)
        #max_err.append(np.max(fitted_vector))
        #print('Maximum difference = ',max_err[-1])
        #if np.all(fitted_vector < fit_tol):
        if norm < fit_tol:
            print('*'*20,' PASS ','*'*20)
        else:
            raise Exception('check')

print('='*48)
print('ALL PASS!')
#print(max_err)


