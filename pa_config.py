import sys, os
from math import sqrt

# Exact location of PAGEN scripts
PAGEN_HOME = os.environ['HOME']+'/src/simpimc/scripts/pagen'
sys.path.append(PAGEN_HOME)
from GenPairAction import run

# Constants
tau = 0.125 # Time step
L = 100.0 # Box size
D = 3 # physical dimension

# Species
e = {'type': 'e', 'lambda': 0.5, 'Z': -1.0}
p = {'type': 'p', 'lambda': 0.0002723089072243553, 'Z': 1.0}

# Potential
potential = {}
potential['function'] = lambda Z1,Z2,r: Z1*Z2/r
potential['r_min'] = 0.0001 # first grid point
potential['r_max'] = 100. # last grid point
potential['n_grid'] = 1000 # number grid points
potential['grid_type'] = "OPTIMIZED" # grid type (LOG, LINEAR, OPTIMIZED (Ilkka only!), LOGLIN (David only!))
potential['r_paste'] = L/4. # where the LOG and LINEAR grids would be pasted for LOGLIN

# Squarer
squarer = {}
squarer['type'] = "Ilkka" # Ilkka, David, or None
squarer['tau'] = tau # desired timestep of PIMC simulation
squarer['n_d'] = D # dimension
squarer['r_max'] = 100.0 # maximum distance on grid
squarer['n_grid'] = 100 # number of grid points
squarer['grid_type'] = "OPTIMIZED" # grid type (LOG, LINEAR, OPTIMIZED (Ilkka only!), LOGLIN (David only!))
squarer['r_paste'] = L/4. # where the LOG and LINEAR grids would be pasted for LOGLIN
squarer['n_square'] = 33 # total number of squarings to reach lowest temperature
squarer['n_order'] = -1 # order of off-diagonal PA fit: -1 = no fit (direct spline, Ilkka only!), 0 = only diagonal, 1-3 = fit off-diagonal to 1-3 order
squarer['n_temp'] = 1 # number of temperatures for which to calculate the pair action (David only!)

# Long-range breakup
breakup = {}
breakup['type'] = 'OptimizedEwald' # OptimizedEwald, StandardEwald, or None
breakup['n_d'] = D # dimension
breakup['L'] = L # length of box
breakup['tau'] = tau # desired timestep of PIMC simulation
breakup['r_min'] = 0.0001 # first grid point
breakup['r_max'] = sqrt(breakup['n_d'])*breakup['L']/2. # last grid point
breakup['r_cut'] = breakup['L']/2. # r cutoff for ewald
breakup['k_cut'] = 14./(L/2.) # k cutoff for ewald
breakup['n_grid'] = 1000 # number of grid points
breakup['grid_type'] = "OPTIMIZED" # grid type (LOG, LINEAR, OPTIMIZED (Ilkka only!), LOGLIN (David only!))
breakup['r_paste'] = L/4. # where the LOG and LINEAR grids would be pasted for LOGLIN
breakup['n_knots'] = 10 # number of knots in spline (probably fine)
breakup['n_images'] = 10 # Naive check

# Pair action objects
pa_objects = [
{'species_a': e, 'species_b': e, 'potential': potential, 'breakup': breakup, 'squarer': squarer},
{'species_a': e, 'species_b': p, 'potential': potential, 'breakup': breakup, 'squarer': squarer},
{'species_a': p, 'species_b': p, 'potential': potential, 'breakup': breakup, 'squarer': squarer},
]

# Run
run(pa_objects)
