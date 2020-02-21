#!/usr/bin/env python
#
# Shared methods
#
from __future__ import print_function
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl

FIGURES = 'figures'
RESULTS = 'results'

# Figure sizes: 85mm wide (or 180mm for double)
FIGSIZE = (85 / 254, 65 / 254)
FIGSIZE = (6, 4.5)

# Configure matplotlib
params = {
    'axes.grid' : False,
    'lines.markersize' : 6,
    'font.size' : 14,
    'axes.labelsize' : 14,
    'xtick.labelsize' : 12,
    'ytick.labelsize' : 12,
    'legend.fontsize' : 10,
    'legend.facecolor' : '#eeeeee',
    'legend.edgecolor' : '#eeeeee',
    'legend.framealpha' : 1,
    'figure.figsize': FIGSIZE,
    }
pl.rcParams.update(params)

#
# Shared methods
#

def median_and_percentiles(data):
    """
    Calculates the median and the 25th and 75th percentiles of the given data.
    
    Data must be specified such data `data[0]` is the first recording (e.g.
    cell 1), `data[1]` is the second etc.
    
    Returns a tuple ``(median, p25, p75)``
    """
    data = np.asarray(data)
    if (len(data.shape) == 1 and not np.isscalar(data[0])):
        # List of lists of unequal length
        median, p25, p75 = [], [], []
        for column in data:
            median.append(np.median(column))
            p25.append(np.percentile(column, 25))
            p75.append(np.percentile(column, 75))
        return np.array(median), np.array(p25), np.array(p75)
    # Equal lenght entries, use fast numpy method
    median = np.median(data, axis=0)
    p25 = np.percentile(data, 25, axis=0)
    p75 = np.percentile(data, 75, axis=0)
    return median, p25, p75

def mean_and_std(data):
    """
    Calculates the mean and standard deviation of the given data.

    Returns a tuple ``(mean, mean-std, mean+std)``.
    """
    data = np.asarray(data)
    if (len(data.shape) == 1 and not np.isscalar(data[0])):
        # List of lists of unequal length
        mean, std = [], []
        for column in data:
            mean.append(np.mean(column))
            std.append(np.std(column))
        mean = np.array(mean)
        std = np.array(std)
    else:
        # Equal lenght entries, use fast numpy method
        mean = np.mean(data, axis=0)
        std = np.std(data, axis=0)
    return mean, mean-std, mean+std

def fill(x, lower, upper, color='red'):
    """
    Add a plot with upper and lower boundaries to the current figure.
    """
    colors = {
        'red' : 'tab:red',
        'orange' : 'tab:orange',
        'black' : 'k',
        }
    color = colors[color]
    pl.fill_between(np.asarray(x), np.asarray(lower), np.asarray(upper),
        linewidth = 0,
        alpha = 0.25,
        color = color,
    )
    
def save(name):
    """
    Saves the current figure.
    """
    base = os.path.join(FIGURES, name)
    pl.savefig(base + '.png', bbox_inches='tight')
    pl.savefig(base + '.pdf', bbox_inches='tight', format='pdf')
    
    





    
