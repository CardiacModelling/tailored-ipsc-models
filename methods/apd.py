#!/usr/bin/env python
#
# Methods to normalise APs and calculate an APD90 variant
#
from __future__ import print_function
from __future__ import division
import os
import glob
import myokit
import methods
import methods.outward as outward
import numpy as np
import scipy
import scipy.signal
import matplotlib
import matplotlib.pyplot as pl

OPTICAL = 'optical'

DRUG_DIR = {
    'dofetilide'    : 'Drug1',
    'quinidine'     : 'Drug2',
    'sotalol'       : 'Drug3',
    'sotalol-okada' : 'Drug3',
    'verapamil'     : 'Drug4',
    'paracetamol'   : 'Drug5',
    }
DRUG_BASE_DIR = {
    'dofetilide'    : 'Group I',
    'quinidine'     : 'Group II',
    'sotalol'       : 'Group III',
    'sotalol-okada' : 'Group III',
    'verapamil'     : 'Group IV',
    'paracetamol'   : 'Group V',
    }

DEBUG_FIGURES = False
AP_FIGURES = os.path.join(methods.FIGURES, 'aps')
APD_FIGURES = os.path.join(methods.FIGURES, 'apds')

def normalise(t, v, use_start_for_baseline=True, correct_time=False,
        filename=None):
    """
    Returns a normalised version of ``v``.
    
    Requires some signal _before_ the AP to be present.
    
    If a filename is given, a debug plot will be stored at the given location.
    """
    # Reliably find point within the upstroke
    v_naive = (v - np.min(v)) / (np.max(v) - np.min(v))
    iup = np.where(v_naive > 0.5)[0][0]
    if iup < 1:
        # Clever stuff didn't work! Possibly depolarised at start of signal
        ipre = 0
        vmin = np.min(v)
    else:
        # Find vmin as the mean of the signal before the upstroke (with a safety
        # margin)
        ipre = int(0.5 * iup)
        vmin = np.mean(v[:ipre])
        
    # Alternative method for files known to have a moving baseline issue
    if not use_start_for_baseline:
        # Set baseline as lower percentile
        vmin = np.percentile(v, 10)

    # Find end of AP (with safety margin)
    e = None
    idown = np.where(v_naive > 0.9)[0][0]
    try:
        idown += np.where(v_naive[idown:] < 0.3)[0][0]
        # Find vmax as 95th percentile during AP
        vmax = np.percentile(v[iup:idown], 95)
    except IndexError as e:
        # Clever stuff didn't work!
        idown = -1
        vmax = np.max(v)

    # Show points etc.
    if DEBUG_FIGURES and filename:
        pl.figure(figsize=(methods.FIGSIZE[0]*2, methods.FIGSIZE[1]*2))
        pl.title('Normalisation')
        pl.suptitle('"AP" region is not used for APD, only for detection of'
            ' vmax')
        pl.plot(t, v)      
        pl.axvline(t[ipre], color='tab:blue', label='Pre-upstroke area')
        pl.axvline(t[iup], color='tab:green', label=' Start AP area')
        pl.axvline(t[idown], color='tab:green', label='End AP area')        
        pl.axhline(vmin, color='tab:orange', label='vmin (mean pre-upstroke)')
        pl.axhline(vmax, color='tab:orange', label='vmax (95th percentile of'
            ' AP area)')
        pl.legend()
        pl.savefig(filename, bbox_inches='tight')
        pl.close()

    # Normalise
    vnorm = (v - vmin) / (vmax - vmin)
    
    if not correct_time:
        return vnorm

    # Some experimental files have the upstroke at 100ms instead of 50ms, this
    # can mess up visualisation and calculation of medians/means (although apd
    # routine works fine).
    # Correct this:
    iup = np.where(vnorm > 0.5)[0][0]
    if t[iup] > 75:
        ichop = np.where(t > 50)[0][0]
        t = t[:-ichop]
        vnorm = vnorm[ichop:]
    
    # Return
    return t, vnorm

def lowpass(t, v):
    """
    Applies a low-pass filter to the given voltage data (with fixed
    parameters)
    """
    f = 1e3 / (t[1] - t[0])
    g = 100
    m = int(f / g)
    n = len(v)
    x = np.zeros(n + 2 * m)
    for i in xrange(2*m+1):
        x[i:i+n] += v
    return x[m:m+n] / float(2*m+1)

def calculate_apd(t, v, pt=0.1, filename=None):
    """
    Calculates an APD.
    
    If a filename is given, a debug plot will be stored at the given location.
    """
    alternans = False
    
    # Detect big upstroke
    imin = np.where(v > 0.8)[0][0]    
    
    # Check for alternans
    if imin == 0:
        # Already depolarised: Sign of alternans!
        alternans = True
    if imin > 0.4 * len(v):
        # Big AP much too late in signal: sign of alternans!
        alternans = True
    
    # Detect end of ap
    vfilter = lowpass(t, v)
    try:
        imax = np.where(np.logical_and(t >= t[imin], vfilter < pt))[0][0]
    except IndexError:
        # APD longer than given signal: Sign of alternans!
        alternans = True
        imax = -1
    apd = t[imax] - t[imin]
    
    # Detect second upstroke
    try:
        # Start of upstroke
        i1 = imax + np.where(v[imax:] > 0.35)[0][0]
        try:
            # Only works if second upstroke is big
            i2 = i1 + np.where(v[i1:] > 0.8)[0][0]
        except IndexError:
            # Second upstroke doesn't get very far: sign of alternans!
            alternans = True
    except IndexError:
        pass
    
    # Store debug output
    if DEBUG_FIGURES and filename:
        pl.figure(figsize=(methods.FIGSIZE[0]*2, methods.FIGSIZE[1]*2))
        pl.title('APD = ' + str(apd) + (' (alternans)' if alternans else ''))
        pl.grid(True)
        pl.plot(t, v)
        pl.plot(t, vfilter)        
        pl.axhline(pt, label=str(pt), color='tab:green')
        pl.axvline(t[imin], label='AP start', color='tab:purple')
        if apd > 0:
            pl.axvline(t[imax], label='AP end', color='tab:red')
        pl.legend(loc='upper right')
        pl.savefig(filename, bbox_inches='tight')
        pl.close()

    return apd, alternans

def load_base_aps(force=False):
    """
    Loads the baseline APs from the optical mapping data.
    """
    cachefile = os.path.join(methods.RESULTS, 'exp-aps-base.csv')
    if not force and os.path.isfile(cachefile):
        print('Using cached baseline optical APs')
        log = myokit.DataLog.load_csv(cachefile).npview()
        time = log.time()
        aps = [log['ap', i] for i in xrange(len(log) - 1)]
        return time, np.array(aps)

    print('Collecting baseline optical APs!')
    files = glob.glob('optical/Group */Well*.csv')
    aps = []
    log = myokit.DataLog()
    log.set_time_key('time')
    for i, filename in enumerate(files):
        data = np.loadtxt(filename, delimiter=',')
        
        # Split into time and voltage
        t = data[:,0]*1000
        v = data[:,1]

        # Normalise
        filename = os.path.join(AP_FIGURES, 'exp-base-' + str(i+1) + '.png')
        t, v = normalise(t, v, filename=filename, correct_time=True)

        # Select bit that's present in all traces
        imin = np.where(t >= 0)[0][0]
        imax = np.where(t >= 575)[0][0]
        t = t[imin:imax]
        v = v[imin:imax]
        
        # Store
        aps.append(v)
        log['ap', i] = v
    
    # Store and return
    log['time'] = t
    log.save_csv(cachefile)    
    return t, np.array(aps)
    
def calculate_base_apds(time, aps, pt=0.1, force=False):
    """
    Calculates APDs for baseline optical experiments.
    """
    cachefile = os.path.join(methods.RESULTS, 'exp-apds-base-' + str(pt)
        + '.csv')
    if not force and os.path.isfile(cachefile):
        print('Using cached baseline optical APDs')
        log = myokit.DataLog.load_csv(cachefile).npview()
        return log['apds']

    print('Collecting baseline optical APDs!')
    apds = []
    for i, ap in enumerate(aps):
        filename = os.path.join(APD_FIGURES, 'exp-base-' + str(i) + '.png')
        apds.append(calculate_apd(time, ap, pt, filename=filename)[0])
    apds = np.array(apds)
    
    # Store and return
    log = myokit.DataLog()
    log['apds'] = apds
    log.save_csv(cachefile)
    return apds

def calculate_tailored_apds(time, aps, pt=0.1, force=False):
    """
    Calculates the tailored simulation aps using the results from
    :meth:`normalise_tailored_simulations`.
    """
    cachefile = os.path.join(methods.RESULTS, 'sim-base-apds-' + str(pt)
        + '.csv')
    if not force and os.path.isfile(cachefile):
        print('Using cached baseline tailored APDs')
        log = myokit.DataLog.load_csv(cachefile).npview()
        return log['apds']

    print('Collecting baseline tailored APDs!')
    apds = []
    for i, long_id in enumerate(outward.ORDER):
        filename = os.path.join(APD_FIGURES, 'sim-base-' + long_id + '.png')
        apds.append(calculate_apd(time, aps[i], pt, filename=filename)[0])
    apds = np.array(apds)
    
    # Store and return
    log = myokit.DataLog()
    log['apds'] = apds
    log.save_csv(cachefile)
    return apds    

def load_drug_apds(drug, pt=0.1, force=False):
    """
    Calculates APDs for optical experiments with drugs.
    """
    # Check drug variable by getting path
    drug = str(drug)
    path = os.path.join(OPTICAL, DRUG_DIR[drug])
    
    # Load cached if possible
    cachefile = os.path.join(methods.RESULTS, 'exp-apds-' + drug + '-' 
        + str(pt) + '.csv')
    if not force and os.path.isfile(cachefile):
        print('Using cached optical APDs for ' + drug)
        log = myokit.DataLog.load_csv(cachefile).npview()
        return [log['apds', i] for i in xrange(len(log))]

    # Get drug base apds
    print('Collecting baseline optical APs for ' + drug)
    files = os.path.join('optical', DRUG_BASE_DIR[drug],'Well*.csv')
    files = glob.glob(files)[:5] # Use five only!
    group = []
    for i, filename in enumerate(files):
        # Load ap
        t, ap = np.loadtxt(filename, delimiter=',').T
        t *= 1e3
        
        # Normalise, calculate apd (figures are already made!)
        t, ap = normalise(t, ap, correct_time=True)
        apd, alternans = calculate_apd(t, ap, pt)
        
        group.append(apd)

    # Get drug applied APDs
    print('Calculating APDs for ' + drug)
    apds = [group]
    group = []
    for i, filename in enumerate(sorted(os.listdir(path))):
        # Load ap
        t, ap = np.loadtxt(os.path.join(path, filename), delimiter=',').T
        t *= 1e3
        
        # Some files need an alternative normalisation method due to a moving
        # baseline!
        use_normal_method = True
        if drug == 'quinidine' and filename in ['Well16.csv', 'Well19.csv']:
            print('>>> Using alternative baseline detection method for '
                + drug + ' file ' + filename)
            use_normal_method = False
        
        # Normalise
        filename = os.path.join(AP_FIGURES, 'exp-' + drug + '-' + str(i+1)
            + '.png')
        t, ap = normalise(t, ap, filename=filename, correct_time=True,
            use_start_for_baseline=use_normal_method)
        
        # Calculate apd
        filename = os.path.join(APD_FIGURES, 'exp-' + drug + '-' + str(i+1)
            + '.png')
        apd, alternans = calculate_apd(t, ap, pt, filename=filename)
        
        group.append(apd)        
        if i % 5 == 4:
            apds.append(group)
            group = []
    if i % 5 != 4:
        raise Exception('Expected drug AP files to be multiple of 5.')
    
    # Store
    log = myokit.DataLog()
    for i, group in enumerate(apds):
        log['apds', i] = group
    log.save_csv(cachefile)

    return apds









