#!/usr/bin/env python
#
# Python code to generate Figure 
#
from __future__ import print_function
import os
import sys
import methods
import methods.apd as apd
import methods.outward as outward
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
import matplotlib.gridspec as gs

#
# Get debug/caching options
#
force = '--force' in sys.argv
if '--apfigures' in sys.argv:
    apd.DEBUG_FIGURES = True

#
# Point at which the APD is measured
#
if '--apd50' in sys.argv:
    PT = 0.5
else:
    PT = 0.1    

#
# Detailed simulations for higher resolution?
#
FAST = '--fast' in sys.argv

#
# Show median and percentiles or mean and standard deviation
#
ALL = False
MEDIAN = True

#
# Show results for Original/O'Hara/Tailored models
#
EXPERIMENT  = True
ORIGINAL    = True
OHARA       = True
TAILORED    = True

#
# Fit all outward current experiments (or used cached result)
#
table, cells = outward.fit_all(force=force)

#
# Get baseline experimental APDs
#
t, aps = apd.load_base_aps(force=force)
base_apds = apd.calculate_base_apds(t, aps, PT, force=force)
if MEDIAN:
    base_mid, base_lo, base_hi = methods.median_and_percentiles(base_apds)
else:
    base_mid, base_lo, base_hi = methods.mean_and_std(base_apds)
del(t, aps)

#
# Create figure for each drug
#
drugs = [
    'dofetilide',
    'quinidine',
    'sotalol',
    'sotalol-okada',
    'verapamil',
    'paracetamol',
    ]
    
drug_concentrations = [
    [0, 0.03, 0.1, 0.3, 1],
    [0, 0.01, 0.1, 1, 10],
    [0, 0.3, 3, 30, 300],
    [0, 0.3, 3, 30, 300],
    [0, 0.01, 0.1, 1, 10],
    [0, 0.3, 3, 30, 300],
    ]

cycle_lengths = [
    1.375,
    1.176,
    0.933,
    0.933,
    0.905,
    1,
    ]

if PT == 0.1:
    figures = [
        'figure5a',
        'figure5b',
        'figure5c',
        'figure5c-okada',
        'figure5d',
        'figure5e',
        ]
else:
    figures = [
        'figure105a',
        'figure105b',
        'figure105c',
        'figure105c-okada',
        'figure105d',
        'figure105e',
        ]
    
legends = [
    'upper right',
    'upper right',
    'upper right',
    'upper right',
    'upper right',
    'upper right',
    ]

ncells = len(outward.ORDER)
for idrug, drug in enumerate(drugs):
    print('== ' + drug + ' ==========================')

    # Get concentrations used in experiment
    concentrations = drug_concentrations[idrug]
    
    # Get cycle-length corresponding to spontaneous rate in experiment
    cl = cycle_lengths[idrug]
    
    # Create detailed list of concentrations for simulations
    if FAST:
        detailed = concentrations
    else:
        detailed = [0]
        for i in xrange(2, len(concentrations)):
            detailed += list(np.logspace(
                np.log10(concentrations[i-1]), 
                np.log10(concentrations[i]),
                4,
                endpoint=False))
        detailed += concentrations[-1:]

    #
    # Set up figure with split x-axis
    #
    fig = pl.figure()
    grid = gs.GridSpec(1, 2, width_ratios=[1,4])
    grid.update(wspace=0.05, hspace=0.5)
    ax1 = fig.add_subplot(grid[0])
    ax2 = fig.add_subplot(grid[1], sharey=ax1)
    ax1.grid(True)
    ax2.grid(True)
    
    # Set up left x-axis
    ax1.set_ylabel('APD$_{' + ('90' if PT == 0.1 else '50') + '}$ [ms]')
    ax1.set_xlim(-5e-3, detailed[1] * 0.85)
    ax1.set_xticks([0])
    
    # Set up right x-axis
    ax2.set_xlabel(drug.capitalize() + ' [$\mu$M]' + ' '*14)
    ax2.set_xscale('log')
    ax2.set_xlim(0.9 * detailed[1], 1.25 * concentrations[-1])
    ax2.set_xticks(concentrations[1:])
    ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    # Remove right y-axis on left axes, left one on right
    ax1.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax1.yaxis.tick_left()
    ax2.yaxis.tick_right()
    ax2.tick_params(labelright='off')
    
    # Add diagonal lines indicating a break
    d = .015
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((1-d,1+d), (-d,+d), **kwargs)
    ax1.plot((1-d,1+d),(1-d,1+d), **kwargs)
    kwargs.update(transform=ax2.transAxes)
    ax2.plot((-d/3,+d/3), (1-d,1+d), **kwargs)
    ax2.plot((-d/3,+d/3), (-d,+d), **kwargs)
        
    # Create plotting version of concentration list for ax2
    # Instead of point 1 being at zero (infinite distance, so horizontal graph
    # regardless of the value), set it for roughly equal spacing of the points
    plot_concentrations = list(concentrations)
    plot_concentrations[0] = plot_concentrations[1] / 10
    plot_detailed = list(detailed)
    plot_detailed[0] = plot_detailed[1] / 10
    
    #
    # Load experimental data
    #
    if EXPERIMENT:
        # Load apds
        exp_apds = np.array(apd.load_drug_apds(drug, PT, force=force)).T
        averages = []
        label = 'Optical mapping (n=' + str(exp_apds.shape[0]) + ')'
        for i, concentration in enumerate(concentrations):
            try:
                y = exp_apds[:,i]
            except IndexError:
                # Concentration can be too high for real experiment
                break
            x = [concentration]*len(y)
            ax1.plot(x, y, 'o', color='k', label=label)
            ax2.plot(x, y, 'o', color='k', label=label)
            averages.append(np.median(y) if MEDIAN else np.mean(y))
            label = None
        ax1.plot(concentrations[:len(averages)], averages, 'k--')
        ax2.plot(concentrations[:len(averages)], averages, 'k--')

    #
    # Run simulations with original model
    #
    if ORIGINAL:
        apds = []
        for concentration in detailed:
            print('Running for untailored model')
            
            # Run simulations
            sim = outward.simulate_original_aps(drug=drug,
                concentration=concentration, cl=cl, force=force)

            # Normalise and get apd
            t = sim.time() * 1e3
            filename = None
            if concentration > 0 and drug != 'paracetamol':
                filename = os.path.join(apd.AP_FIGURES, 'sim-original-' + drug
                    + '-' + str(concentration) + '.png')
            ap = apd.normalise(t, sim['membrane.V'], filename=filename)
            if concentration > 0 and drug != 'paracetamol':
                filename = os.path.join(apd.APD_FIGURES, 'sim-original-' + drug
                    + '-' + str(concentration) + '.png')
            ad, alternans = apd.calculate_apd(t, ap, PT, filename)
            
            # Check for alternans
            if alternans or ad*1e-3 > 0.95 * cl:
                print('>>> Alternans and/or APD exceeding CL found in original'
                    +' model for ' + drug + ' at ' + str(concentration) + ' uM')
                break
            
            # Store
            apds.append(ad)
        del(ap, ad)
        for i in xrange(len(detailed) - len(apds)):
            apds.append(float('NaN'))

        # Plot
        ax1.plot(detailed[:2], apds[:2], '-', color='tab:blue')
        ax2.plot(plot_detailed, apds, '-', color='tab:blue',
            label='Original model')

    #
    # Run simulations with O'Hara model
    #
    if OHARA:
        apds = []
        for concentration in detailed:
            print('Running for O\'Hara model, ' + str(concentration))
            
            # Run simulations
            sim = outward.simulate_ohara_aps(drug=drug,
                concentration=concentration, cl=cl, force=force)

            # Normalise and get apd
            t = sim.time() * 1e3
            filename = None
            if concentration > 0 and drug != 'paracetamol':
                filename = os.path.join(apd.AP_FIGURES, 'sim-ohara-' + drug
                    + '-' + str(concentration) + '.png')
            ap = apd.normalise(t, sim['membrane.V'], filename=filename)
            if concentration > 0 and drug != 'paracetamol':
                filename = os.path.join(apd.APD_FIGURES, 'sim-ohara-' + drug
                    + '-' + str(concentration) + '.png')
            ad, alternans = apd.calculate_apd(t, ap, PT, filename)
            
            # Check for alternans
            if alternans or ad*1e-3 > 0.95 * cl:
                print('>>> Alternans and/or APD exceeding CL found in O\'Hara'
                    +' model for ' + drug + ' at ' + str(concentration) +' uM')
                break
            
            # Store
            apds.append(ad)
        del(ap, ad)
        for i in xrange(len(detailed) - len(apds)):
            apds.append(float('NaN'))

        # Plot
        ax1.plot(detailed[:2], apds[:2], '-', color='tab:green')
        ax2.plot(plot_detailed, apds, '-', color='tab:green',
            label='O\'Hara model')

    #
    # Run simulations with tailored models
    #
    if TAILORED:
        apds = []
        stars = []
        nmodels = []
        for concentration in detailed:
            print('Running for drug: ' + drug + ' ' + str(concentration) + ' [uM]')
      
            # Run simulations    
            sim = outward.simulate_tailored_aps(table, cells, drug=drug,
                concentration=concentration, cl=cl, force=force)
                
            # Normalise and get apds
            group = []
            t = sim.time() * 1e3
            for i, long_id in enumerate(outward.ORDER):
                filename = None
                if concentration > 0 and drug != 'paracetamol':
                    filename = os.path.join(apd.AP_FIGURES, 'sim-' + long_id + '-'
                        + drug + '-' + str(concentration) + '.png')
                ap = apd.normalise(t, sim['membrane.V', i], filename=filename)
                if concentration > 0 and drug != 'paracetamol':        
                    filename = os.path.join(apd.APD_FIGURES, 'sim-' + long_id + '-'
                        + drug + '-' + str(concentration) + '.png')
                ad, alternans = apd.calculate_apd(t, ap, PT, filename)
                
                # Check for alternans
                if alternans or ad*1e-3 > 0.95 * cl:
                    print('>>> Alternans and/or APD exceeding CL found in tailored'
                        + ' model for ' + drug + ' at ' + str(concentration)
                        + ' uM, cell ' + long_id)
                    if ALL:
                        group.append(float('nan'))
                    else:
                        stars.append((concentration, ad))
                else:
                    # Store
                    group.append(ad)
                
            # Count models in group
            nmodels.append(len(group))
            
            # Store group
            apds.append(group)
            del(ap, group)
        apds = np.array(apds).T
        
        # Get median and percentiles
        if MEDIAN:
            mid, lo, hi = methods.median_and_percentiles(apds)
        else:
            mid, lo, hi = methods.mean_and_std(apds)

        # Plot
        if ALL:
            for trace in apds:
                ax1.plot(plot_detailed, trace, color='tab:orange', alpha=0.5)
                ax2.plot(plot_detailed, trace, color='tab:orange', alpha=0.5)
        else:
            ax1.fill_between(detailed[:2], lo[:2], hi[:2],
                linewidth = 0, alpha = 0.25, color = 'tab:orange')
            ax2.fill_between(plot_detailed, lo, hi,
                linewidth = 0, alpha = 0.25, color = 'tab:orange')
        
            ax1.plot(detailed[:2], mid[:2], '-', color='tab:orange')
            ax2.plot(plot_detailed, mid, '-', color='tab:orange',
                label='Tailored models (' + ('median' if MEDIAN else 'mean')
                    + ', n=' + str(min(nmodels)) + '-' + str(max(nmodels))
                    + ')')
    
        # Show points with alternans
        if False:
            for star in set(stars):
                x = star[0] * 0.95
                y = mid[detailed.index(star[0])] + 25    
                pl.plot(x, y, '*', color='tab:orange')
    
    # Set y-axis limits
    if PT == 0.5:
        ax1.set_ylim(25, 1000)
        ax2.set_ylim(25, 1000)
    else:
        ax1.set_ylim(100, 1100)
        ax2.set_ylim(100, 1100)

    # Fix legend for quinidine figure
    if drug == 'quinidine':
        pl.legend(bbox_to_anchor=(0.58, 0.99))
    else:
        pl.legend(loc=legends[idrug])
        
    # Store figure
    methods.save(figures[idrug])

print('Figures saved to ' + methods.FIGURES + ' (plus apds)')

#pl.show()


