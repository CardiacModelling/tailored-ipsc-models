#!/usr/bin/env python
#
# Methods to fit the outward current experiment
#
from __future__ import print_function
from __future__ import division
import os
import myokit
import myokit.lib.fit as fitlib
import methods
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl

# Some versions of CMA enable matplotlib's interactive mode
interactive = matplotlib.is_interactive()
#import cma
matplotlib.interactive(interactive)
del(interactive)

EXP = 'outward-experiment'
SIM = os.path.join(methods.RESULTS, 'sim-outward-currents.csv')
TABLE = os.path.join(methods.RESULTS, 'sim-outward-current-table.csv')

# Log interval
DT = 0.2e-3

# Voltage steps
VSTEPS_mV = np.arange(-40, 50 + 10, 10)

# Offset ignored when fitting at the start of each trace
CAPACITANCE = 0.01      # 10 ms

#
# List of cells, ordered as in the excel sheet?
#
ORDER = [
    'H3_KIV_IKs',
    'G10_KIV_IKs',
    'G7_KIV_IKs',
    'F8_KIV_IKs',
    'E10_KIV_IKs',
    'E8_KIV_IKs',
    'E7_KIV_IKs',
    'E3_KIV_IKsr',
    'D9_KIV_IKsr',
    'D8_KIV_IKsr',
    'C9_KIV_IKsr',
    'C3_KIV_IKsr',
    'C1_KIV_IKs',
    'B10_KIV_IKs',
    'A10_KIV_IKsr',
    'A9_KIV_IKsr',
    'A8_KIV_IKsr',
    'A3_KIV_IKs',
    'B3_KIV_IKsr',
    'B9_KIV_IKtor',
    'E1_KIV_IKtos',
    'A7_KIV_IKs',
    ]

CURRENTS = [
    'ina.INa',
    'ical.ICaL',
    'ik1.IK1',
    'ikr.IKr',
    'iks.IKs',
    'ito.Ito',
    'if.If',
    'inaca.INaCa',
#    'inak.INaK',
    ]
# Note: IpCa=0 when Cai is fixed to 0
# Note: IbCa=Inf when Cai is fixed to 0
# Note: IbNa (defined as `g * (V - ENa)` becomes very large if we add it to the
#       simulation, and takes over the role of IKr.

# Table for fit all
FIELDS = [
    'index',
    'long_id',
    'short_id',
    'score',
    ]
FIELDS += ['scale.' + current for current in CURRENTS]
FIELDS += ['contribution.' + current for current in CURRENTS]

def count():
    """
    Returns the number of cells.
    """
    return len(order)

def find(string):
    """
    Returns the indice for the given string.
    """
    letter, number = string[0], int(string[1:])
    string = letter.upper() + str(number)
    for index, long_id in enumerate(ORDER):
        short_id = long_id.split('_')[0]
        if short_id == string:
            return index
    raise ValueError('No cell found for: "' + str(string) + '"')

def load(index):
    """
    Loads the outward current data for a single cell.
    """
    # Get long and short id
    long_id = ORDER[index]
    short_id = long_id.split('_')[0]
    letter, number = short_id[0], int(short_id[1:])
    
    # Load capacitance
    '''
    cm = None
    with open(os.path.join(ROOT, 'cm.csv'), 'r') as f:
        for line in f:
            i, cm = line.split(',')
            let, num = i[0], int(i[1:])
            if let==letter and num==number:
                break
            cm = None
    if cm is None:
        raise ValueError('Capacitance not found for ' + long_id)
    cm = float(cm)
    '''
            
    # Load data
    data = []
    with open(os.path.join(EXP, long_id) + '.txt', 'r') as f:
        for line in f:
            data.append([float(x) for x in line.split()])
    data = np.array(data).T
    
    # Separate time and current
    time = np.array(data[0], copy=True)
    data = np.array(data[1:], copy=True)
    
    # Set time to seconds (like the model etc.)
    # And fix offset
    time /= 1000000
    time -= 0.050
    
    # Normalise current
    #data /= cm
    
    return time, data

def trim(time, data):
    """
    Trims the data to just the step part.
    """
    lo, hi = int(round(0.05/DT)), -int(round(0.05/DT))
    return time[lo:hi], data[:,lo:hi]

def plot(time, data):
    """
    Creates a figure of some loaded data.
    """
    n = len(VSTEPS_mV)
    if len(data) != n:
        raise ValueError('Expected data of size ' + str(n))

    # Create colormap
    cmap = matplotlib.cm.get_cmap('viridis')
    norm = matplotlib.colors.Normalize(np.min(VSTEPS_mV), np.max(VSTEPS_mV))

    # Figure out good y-axes limits
    t2, d2 = trim(time, data)
    upper = np.max(np.max(d2[:,1000:]))
    upper = 2 * np.ceil(upper * 0.5)
    lower = np.min(np.min(d2[:,50:]))
    lower = max(-4, 2 * np.floor(lower * 0.5))

    # Create main figure
    fig = pl.figure()
    pl.subplot(1,1,1)
    pl.xlim(-15, 515)
    pl.ylim(lower, upper)
    #pl.ylim(-2, 15)
    for i, trace in enumerate(data):
        pl.plot(time*1000, trace, lw=1, color=cmap(norm(VSTEPS_mV[i])))
    pl.xlabel('Time [ms]')
    pl.ylabel('I$_{outward}$ [A/F]')

    # Add color bar
    fig.subplots_adjust(right=0.875)
    ax = fig.add_axes([0.9, 0.15, 0.015, 0.7])
    cb = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm)
    cb.set_label('V$_{hold}$ [mV]')

def simulate_outward_current(force=False):
    """
    Simulates outward current using the original model.
    
    Won't simulate if the traces are already present. Set `force` to `True` to
    force new simulations.
    """
    def join(log):
        time = np.array(log.time(), copy=True)
        data = {}            
        n = len(VSTEPS_mV)
        for current in CURRENTS:
            data[current] = np.array(
                [log[current, i] for i in xrange(n)], copy=True)
        return time, data
    
    if not force:
        if os.path.isfile(SIM):
            log = myokit.DataLog.load_csv(SIM)
            ok = True
            for current in CURRENTS:
                if '0.' + current not in log:
                    ok = False
                    break
            if ok:
                print('Using cached simulated currents!')
                return join(log)
    print('Running new untailored currents simulation')
    
    # Load model
    model = myokit.load_model('paci-2013-ventricular.mmt')

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
    vsteps = VSTEPS_mV * 1e-3
    tpre = 10 # 10 s
    tstep = 0.5 # 500ms
    sodium_protocol = myokit.pacing.steptrain(vsteps, vhold, tpre, tstep)
    duration = sodium_protocol.characteristic_time()

    # Clamp membrane potential to protocol
    vm = model.get('membrane.V')
    vm.demote()
    vm.set_rhs(-80e-3)
    model.get('engine.pace').set_binding(None)
    vm.set_binding('pace')

    # Apply protocol
    sim = myokit.Simulation(model, sodium_protocol)
    data = sim.run(duration, log=['engine.time'] + CURRENTS, log_interval=DT)

    # Create overlapping step data
    data = data.fold(tpre + tstep)
    data = data.trim_left(tpre, adjust=True)

    # Store data
    data.save_csv(SIM)
    
    # Return data
    return join(data)

def fit(index, force=False):
    """
    Attempts to fit the outward current with some linear combination of
    simulated currents.
    
    Arguments:
    
    ``index``
        The index of the cell to fit
    ``force``
        Set to ``True`` to force new simulations to be made
    
    Returns:
    
    ``score``
        The final error score
    ``scales``
        The scaling factor for each current.
    ``contribution``
        The contribution score for each current
    ``time``
        The times array.
    ``exp``
        The experimental current as a 2d array with indices (repeat, time).
    ``sim``
        The simulated outward current as a 2d array with indices (repeat,
        time).

    """
    # Load experimental data
    etime_org, exp_org = trim(*load(index))
    
    # Load simulated original currents
    stime_org, sim_org = simulate_outward_current(force)
   
    # Chop off capacitance artefacts
    icapacitance = int(round(CAPACITANCE / DT))
    etime = etime_org[icapacitance:]
    stime = stime_org[icapacitance:]
    exp = exp_org[:, icapacitance:]
    sim = {}
    for current in CURRENTS:
        sim[current] = sim_org[current][:, icapacitance:]
    
    # Plot experimental traces with start chopped off
    if False:
        pl.figure()
        for i, trace in enumerate(exp):
            pl.plot(etime, trace)
        pl.show()
    
    # Get 3d blocks of data
    sim = np.array([sim[current] for current in CURRENTS], copy=True)
    sim_org = np.array([sim_org[current] for current in CURRENTS], copy=True)
    
    # Number of currents
    n = len(CURRENTS)
    
    # Number of voltages
    m = len(VSTEPS_mV)
    
    # Number of points per step
    k = len(stime)
    
    # Define score function
    def score(scales):
        scales = np.asarray(scales).reshape((n, 1, 1))
        #error = np.sum(scales * sim, axis=0) - exp
        #error = error.reshape((m, k))
        #return np.sum(np.linalg.norm(error, axis=1))
        #return np.linalg.norm(error)
        return np.sum((np.sum(scales * sim, axis=0) - exp)**2)

    # Set boundaries
    bounds = [(0, 1000) for i in xrange(n)]
    
    # Optimise!
    scales, fx = fitlib.cmaes(score, bounds, parallel=True, verbose=True)
    #scales, fx = fitlib.xnes(score, bounds, max_iter=1000, parallel=True,
    #    verbose=100)

    # Calculate scaled currents and total outward current
    sim_scaled = scales.reshape((n, 1, 1)) * sim_org
    sim_outward = np.sum(sim_scaled, axis=0)
    
    # Calculate contribution scores
    contribution = np.abs(sim_scaled[:,-1,-1])
    contribution /= np.sum(contribution)
    
    # Filter contributions
    contribution[contribution < 1e-6] = 0
    contribution *= 100

    # Print scales and contribution scores
    print('--' + ORDER[index].split('_')[0] + '-' * 40)
    print('Current, Scale, Contribution')
    m1 = max([len(i) for i in CURRENTS]) + 2
    m2 = max([len(i) for i in [str(scale) for scale in scales]]) + 2
    for i, current in enumerate(CURRENTS):
        scale = str(scales[i])
        score = str(contribution[i])
        print(
            current + ' '*(m1 - len(current)) + 
            scale + ' '*(m2 - len(scale)) +
            score
            )
    print('Score: ' + str(fx))    

    # Return lots of data!
    return (
        fx,
        scales,        
        contribution,
        stime_org,
        exp_org,
        sim_outward,
        )

def plot_experiment(score, scales, contributions, time, exp, sim):
    """
    Takes results from fit() and plots the experimental current.
    """
    # Guess limits for y-axis
    ignore = int(0.015 / DT)
    lower = 2 * int(np.floor(0.5 * np.min(np.min(exp[:,ignore:]))))
    upper = 2 * int(np.ceil(0.5 * np.max(np.max(exp[:,ignore:]))))

    # Create colormap
    cmap = matplotlib.cm.get_cmap('viridis')
    norm = matplotlib.colors.Normalize(np.min(VSTEPS_mV), np.max(VSTEPS_mV))

    # Plot
    pl.figure()
    pl.xlabel('Time [ms]')
    pl.ylabel('I$_{outward}$ [A/F]')
    pl.ylim(lower, upper)
    for i, trace in enumerate(exp):
        pl.plot(time*1e3, trace, color=cmap(norm(VSTEPS_mV[i])))
    pl.axvline(CAPACITANCE*1e3, color='r')

def plot_simulation(score, scales, contributions, time, exp, sim):
    """
    Takes results from fit() and plots the simulated current.
    """
    # Guess limits for y-axis
    ignore = int(0.015 / DT)
    lower = 2 * int(np.floor(0.5 * np.min(np.min(exp[:,ignore:]))))
    upper = 2 * int(np.ceil(0.5 * np.max(np.max(exp[:,ignore:]))))

    # Create colormap
    cmap = matplotlib.cm.get_cmap('viridis')
    norm = matplotlib.colors.Normalize(np.min(VSTEPS_mV), np.max(VSTEPS_mV))

    # Plot
    pl.figure()
    pl.xlabel('Time [ms]')
    pl.ylabel('I$_{outward}$ [A/F]')
    pl.ylim(lower, upper)
    for i, trace in enumerate(sim):
        vm = VSTEPS_mV[i]
        pl.plot(time*1e3, trace, color=cmap(norm(vm)), label=str(vm)+' mV')
    pl.legend(loc='lower right', ncol=5, columnspacing=1.2, fontsize=9,
        handlelength=1.5, handletextpad=0.5)

def fit_all(force=False):
    """
    Repeats the outward current fitting process for every cell.
    
    Returns two objects:
    
    1. A dict with information, with keys:
    
        ``index``
            The cell index (1,2,3,...)
        ``name``
            A short name ('b3')
        ``score``
            The final fit score
        ``scale.ina.INa`` etc.
            A scale for each current
        ``contribution.ina.INa`` etc.
            A contribution score for each current
    
    2. A list of dictionaries (ordered by cell) with keys:
    
        ``time``
            The sampling points
        ``0.exp``, ``1.exp``, etc.
            The experimentally recorded outward current, per step
        ``0.sim``, ``1.sim``, etc.
            The simulated outward current, per step

    """
    # Assume we can load everything
    refresh = bool(force)
    
    # Create table structure
    table = {}
    for field in FIELDS:
        table[field] = []

    # Load cached info (assume it's correct)
    if refresh or not os.path.isfile(TABLE):
        refresh = True
    else:
        with open(TABLE, 'r') as f:
            lines = f.readlines()
            header = lines[0].strip().split(',')
            if len(lines) < len(ORDER)+1 or header != FIELDS:
                refresh = True
            else:
                try:
                    for line in lines[1:]:
                        parts = line.strip().split(',')
                        for i, field in enumerate(FIELDS):
                            if field == 'index':
                                table[field].append(int(parts[i]))
                            elif field == 'long_id' or field == 'short_id':
                                table[field].append(parts[i])
                            else:
                                table[field].append(float(parts[i]))
                except ValueError:
                    refresh = True

    # Create structure
    cells = []
    
    def cell_filename(cell):
        return os.path.join(methods.RESULTS, 'sim-outward-current-' + cell +
            '.csv')
    
    # Load cached time series data
    if not refresh:
        for cell in ORDER:
            filename = cell_filename(cell)
            try:
                d = myokit.DataLog.load_csv(filename).npview()
            except Exception:
                refresh = True
                break
            if 'time' not in d:
                refresh = True
                break
            for i in xrange(len(VSTEPS_mV)):
                if str(i) + '.exp' not in d or str(i) + '.sim' not in d:
                    refresh = True
                    break
            if refresh:
                break
            cells.append(d)
    
    # Return cached data if possible
    if not refresh:
        print('Using cached fits!')
        return table, cells

    # Fit all cells!
    print('Fitting all cells')
    for i, cell in enumerate(ORDER):
        # Fit
        force_sim = force and i==0 # If forced refreshing, only simulate once
        score, scales, contributions, time, exp, sim = fit(i, force_sim)
        
        # Append data to table
        table['index'].append(i + 1)
        table['long_id'].append(cell)
        table['short_id'].append(cell.split('_')[0])
        table['score'].append(score)
        for j, current in enumerate(CURRENTS):
            table['scale.' + current].append(scales[j])
            table['contribution.' + current].append(contributions[j])
            
        # Create and store cell data file
        d = myokit.DataLog()
        d.set_time_key('time')
        d['time'] = time
        for j in xrange(len(VSTEPS_mV)):
            d['exp', j] = exp[j]
            d['sim', j] = sim[j]
        d.save_csv(cell_filename(cell))
        cells.append(d)
    
    # Store table
    with open(TABLE, 'w') as f:
        f.write(','.join(FIELDS) + '\n')
        for i, cell in enumerate(ORDER):
            row = [table[field][i] for field in FIELDS]
            f.write(','.join([str(x) for x in row]) + '\n')

    # Return new data
    return table, cells
    
def print_table(table):
    """
    Produces a readable version of tables returned by :meth:`fit_all()`.
    """
    # Create big array of strings
    strings = []
    header = []
    for field in FIELDS:
        if field == 'index':
            header.append('#')
        elif field == 'long_id':
            continue
        elif field == 'short_id':
            header.append('id')
        elif field == 'score':
            header.append(field)
        else:
            parts = field.split('.')
            header.append(parts[0][0] + '_' + parts[2])
    strings.append(header)
    for i, cell in enumerate(ORDER):
        row = []
        for field in FIELDS:
            value = table[field][i]
            if field == 'index':
                row.append(str(value))
            elif field == 'score':
                row.append(str(int(round(value))))
            elif field == 'long_id':
                continue
            elif field == 'short_id':
                row.append(value)
            else:
                if value < 1e-10:
                    row.append('~0')
                else:
                    row.append('{:>1.3g}'.format(value))
        strings.append(row)
    strings = np.array(strings)
    
    # Calculate column widths, padding per cell
    lengths = np.array([[len(x) for x in y] for y in strings])
    widths = np.max(lengths, axis=0)
    padding = widths - lengths + 1
    
    # Print!
    for i, row in enumerate(strings):
        for j, text in enumerate(row):
            print(text + ' '*padding[i][j], end='')
        print('')

def tex_table(table):
    """
    Produces a tex version of tables returned by :meth:`fit_all()`.
    """
    fields = [
        'short_id',
        'scale.inaca.INaCa',
        'contribution.inaca.INaCa',
        'spacer',
        'scale.iks.IKs',
        'contribution.iks.IKs',
        'spacer',
        'scale.ikr.IKr',
        'contribution.ikr.IKr',
        'spacer',
        'scale.ik1.IK1',
        'contribution.ik1.IK1',
        'spacer',
        'scale.if.If',
        'contribution.if.If',
        ]
    
    # Create big array of strings
    strings = []
    for i, cell in enumerate(ORDER):
        row = []
        for field in fields:
            if field == 'short_id':
                value = table[field][i]
                row.append(value + ' ')
            elif field == 'spacer':
                row.append('')                
            else:
                value = table[field][i]
                if value < 1e-10:
                    row.append(' --- ')
                else:
                    row.append(' {:>1.3g} '.format(value))
        strings.append(row)
    strings = np.array(strings)
    
    # Calculate column widths, padding per cell
    lengths = np.array([[len(x) for x in y] for y in strings])
    widths = np.max(lengths, axis=0)
    padding = widths - lengths
    
    # Print!
    for i, row in enumerate(strings):
        for j, text in enumerate(row):
            if j > 0:
                print(' '*padding[i][j-1] + '&', end='')
            print(text, end='')
        print(' \\\\')


def select(table, cells, index):
    """
    Takes the data returned by `fit_all` and returns the data for the cell at
    the given `index` (in the format returned by `fit`).
    """
    score = table['score'][index]
    scales = []
    contributions = []
    for current in CURRENTS:
        scales.append(table['scale.' + current][index])
        contributions.append(table['contribution.' + current][index])
    log = cells[index].npview()
    time = log.time()
    exp = np.array([log['exp', i] for i in xrange(len(VSTEPS_mV))])
    sim = np.array([log['sim', i] for i in xrange(len(VSTEPS_mV))])
    return score, scales, contributions, time, exp, sim
    
def simulate_aps(scales, drug=None, concentration=0, beats=2, cl=1,
        stimulate=True, store_currents=False):
    """
    Generates APs using the given scalings.
    
    Scales must be given as a dictionary `variable: scale`.
    """
    # Load model
    model = myokit.load_model('paci-2013-ventricular.mmt')

    # Room temperature
    model.get('phys.T').set_rhs(298)
    
    # Add drug block scalings
    def hill(x, ic50, h=1):
        return 1 / (1 + (x / ic50)**h)
    def add(current, x, ic50):
        f = hill(x, ic50)
        if current in scales:
            scales[current] *= f
        else:
            scales[current] = f
    
    # Add drug scales
    if drug is None or drug == 'paracetamol':
        pass
    
    elif drug == 'dofetilide':
        add('ikr.IKr', concentration, 5.2e-3)
        add('ina.INa', concentration, 147.9)
        add('ical.ICaL', concentration, 26.7)
        add('iks.IKs', concentration, 415.8)
        
    elif drug == 'quinidine':
        add('ikr.IKr', concentration, 0.3)
        add('ina.INa', concentration, 16.6)
        add('ical.ICaL', concentration, 15.6)
    
    elif drug == 'sotalol':
        add('ikr.IKr', concentration, 111.4)
        add('ina.INa', concentration, 7013.9)
        add('ical.ICaL', concentration, 193.3)

    elif drug == 'sotalol-okada':
        add('ikr.IKr', concentration, 356.4)
        #add('ina.INa', concentration, 7013.9)
        #add('ical.ICaL', concentration, 193.3)
    
    elif drug == 'verapamil':
        add('ikr.IKr', concentration, 0.25)
        add('ina.INa', concentration, 32.5)
        add('ical.ICaL', concentration, 0.2)
        
    else:
        raise ValueError('Unknown drug: ' + str(drug))

    # Apply scalings
    for var, scale in scales.iteritems():
        v = model.get(var)
        v.set_rhs(myokit.Multiply(myokit.Number(scale), v.rhs()))
    
    # Simulate with modified model
    sim = myokit.Simulation(model)
    
    # Add stimulus
    if stimulate:
        protocol = myokit.pacing.blocktrain(period=cl, duration=5e-3,
            offset=5e-2)
        sim.set_protocol(protocol)
    
    # Pre-pace for some beats
    sim.pre(100*cl)
    
    # Log some beats and return
    log = ['engine.time', 'membrane.V']
    if store_currents:
        log += ['ina.INa','ical.ICaL','ikr.IKr','iks.IKs','inaca.INaCa']
    return sim.run(beats * cl, log=log, log_interval=DT).npview()
    
def simulate_tailored_aps(table, cells, drug=None, concentration=0, beats=2,
        cl=1, stimulate=True, store_currents=False, force=False):
    """
    Simulates 2 APs for all fitted cells or returns cached results
    """
    concentration = float(concentration)
    if concentration == 0 or drug is None or drug == 'paracetamol':
        filename = os.path.join(methods.RESULTS, 'sim-tailored-aps-base-'
            + str(cl) + '-' + str(beats)
            + '-' + ('paced' if stimulate else 'unpaced')
            + ('-with-currents' if store_currents else '')
            + '.csv')
    else:
        drug = str(drug)
        filename = os.path.join(methods.RESULTS, 'sim-tailored-aps-' + drug
            + '-' + str(concentration) + '-' + str(cl) + '-' + str(beats)
            + '-' + ('paced' if stimulate else 'unpaced')
            + ('-with-currents' if store_currents else '')
            + '.csv')
        
    if not force and os.path.isfile(filename):
        print('Using cached tailored model APs!')
        return myokit.DataLog.load_csv(filename).npview()

    print('Simulating tailored model APs!')
    aps = myokit.DataLog()
    aps.set_time_key('engine.time')
    
    for i, long_id in enumerate(ORDER):
        print('Simulating ' + long_id)

        scale = {}
        
        # Add inward current scales
        scale['ina.INa'] = 0.69     # From sodium experiment
        scale['ical.ICaL'] = 0.80   # From calcium experiment
        
        # Add outward current scales
        #scale['ina.INa'] = table['scale.ina.INa'][i]
        #scale['ical.ICaL'] = table['scale.ik1.ICaL'][i]
        #scale['ik1.IK1'] = table['scale.ik1.IK1'][i]
        #scale['ikr.IKr'] = table['scale.ikr.IKr'][i]        
        scale['iks.IKs'] = table['scale.iks.IKs'][i]
        #scale['ito.Ito'] = table['scale.ito.Ito'][i]
        #scale['if.If'] = table['scale.if.If'][i]
        scale['inaca.INaCa'] = table['scale.inaca.INaCa'][i]
        #scale['inak.INaK'] = table['scale.inak.INaK'][i]
        #scale['ibna.IbNa'] = table['scale.ibna.IbNa'][i]

        log = simulate_aps(scale, drug, concentration, beats, cl, stimulate,
            store_currents)
        aps['engine.time'] = log.time()
        aps['membrane.V', i] = log['membrane.V']
        if store_currents:
            for c in ['ina.INa','ical.ICaL','ikr.IKr','iks.IKs','inaca.INaCa']:
                aps[c,i] = log[c]
    del(log)
    
    aps = aps.npview()
    aps.save_csv(filename)
    
    return aps

def simulate_original_aps(drug=None, concentration=0, force=False, beats=2,
        cl=1):
    """
    Simulates APs with the unaltered model.
    """
    concentration = float(concentration)
    if concentration == 0 or drug is None or drug == 'paracetamol':
        filename = os.path.join(methods.RESULTS, 'sim-original-aps-base-'
            + str(cl) + '-' + str(beats) + '.csv')
    else:
        drug = str(drug)        
        filename = os.path.join(methods.RESULTS, 'sim-original-aps-' + drug
            + '-' + str(concentration) + '-' + str(cl) + '-' + str(beats)
            + '.csv')
    
    if not force and os.path.isfile(filename):
        print('Using cached original APs!')
        return myokit.DataLog.load_csv(filename).npview()

    print('Simulating original APs!')
    aps = simulate_aps({}, drug, concentration, beats, cl)
    aps.save_csv(filename)
    return aps

def simulate_ohara_aps(drug=None, concentration=0, beats=2, cl=1, force=False):
    """
    Simulates APs with the O'Hara model.
    """
    concentration = float(concentration)
    if concentration == 0 or drug is None or drug == 'paracetamol':
        filename = os.path.join(methods.RESULTS, 'sim-ohara-aps-base'
            + '-' + str(cl) + '-' + str(beats) + '.csv')
    else:
        drug = str(drug)        
        filename = os.path.join(methods.RESULTS, 'sim-ohara-aps-' + drug
            + '-' + str(concentration) + '-' + str(cl) + '-' + str(beats)
            +'.csv')
    
    if not force and os.path.isfile(filename):
        print('Using cached O\'Hara APs!')
        return myokit.DataLog.load_csv(filename).npview()

    print('Simulating O\'Hara APs!')
    aps = simulate_ohara_ap(drug, concentration, beats, cl)
    aps.save_csv(filename)
    return aps

def simulate_ohara_ap(drug=None, concentration=0, beats=2, cl=1):
    """
    Generates APs with the O'Hara model.
    """
    # Convert CL to milliseconds
    cl *= 1000
    
    # Load model
    model = myokit.load_model('ohara-2011.mmt')

    # Load protocol
    protocol = myokit.pacing.blocktrain(period=cl, duration=0.5, offset=50)
    
    # Add drug block scalings
    scales = {}
    def hill(x, ic50, h=1):
        return 1 / (1 + (x / ic50)**h)
    def add(current, x, ic50):
        f = hill(x, ic50)
        if current in scales:
            scales[current] *= f
        else:
            scales[current] = f
    
    # Add drug scales
    if drug is None or drug == 'paracetamol':
        pass
    
    elif drug == 'dofetilide':
        add('ikr.IKr', concentration, 5.2e-3)
        add('ina.INa', concentration, 147.9)
        add('ical.ICaL', concentration, 26.7)
        add('iks.IKs', concentration, 415.8)
        
    elif drug == 'quinidine':
        add('ikr.IKr', concentration, 0.3)
        add('ina.INa', concentration, 16.6)
        add('ical.ICaL', concentration, 15.6)
    
    elif drug == 'sotalol':
        add('ikr.IKr', concentration, 111.4)
        add('ina.INa', concentration, 7013.9)
        add('ical.ICaL', concentration, 193.3)

    elif drug == 'sotalol-okada':
        add('ikr.IKr', concentration, 356.4)
        #add('ina.INa', concentration, 7013.9)
        #add('ical.ICaL', concentration, 193.3)
      
    elif drug == 'verapamil':
        add('ikr.IKr', concentration, 0.25)
        add('ina.INa', concentration, 32.5)
        add('ical.ICaL', concentration, 0.2)
        
    else:
        raise ValueError('Unknown drug: ' + str(drug))

    # Apply scalings
    for var, scale in scales.iteritems():
        v = model.get(var)
        v.set_rhs(myokit.Multiply(myokit.Number(scale), v.rhs()))
    
    # Simulate with modified model
    sim = myokit.Simulation(model, protocol)
    
    # Pre-pace for some beats
    p = myokit.ProgressPrinter()
    sim.pre(100*cl, progress=p)
    
    # Log some beats
    p = myokit.ProgressPrinter()
    log = sim.run(beats * cl, log=['engine.time', 'membrane.V'],
        log_interval=DT*1e3, progress=p).npview()
    
    # Convert to seconds and volts
    log['engine.time'] *= 1e-3
    log['membrane.V'] *= 1e-3
    
    # Return
    return log

def load_outward_peaks(force=False):
    """
    Loads the peak outward current per cell.
    """
    # Return cached version
    cachefile = os.path.join(methods.RESULTS, 'exp-peaks.csv')
    if not force and os.path.isfile(cachefile):
        print('Using cached peak outward current data!')
        log = myokit.DataLog.load_csv(cachefile).npview()
        return log['peaks']
    print('Gathering peak outward current data!')
    
    peaks = []
    for i, long_id in enumerate(ORDER):
        t, v = load(i)
        peaks.append(np.max(v))

    # Store
    log = myokit.DataLog()
    log['peaks'] = peaks
    log.save_csv(cachefile)
    
    # Return
    return peaks
    
        























