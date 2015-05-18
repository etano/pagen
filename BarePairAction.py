import sys, os
from math import sqrt

# Units
units = {'energy':'H', 'distance':'A'}

# Constants
tau = 0.125 # Time step
L = 100.0 # Box size

# Particles
particles = [{'type': 'e', 'lambda': 0.5, 'Z': -1.0},
             {'type': 'p-importance', 'lambda': 0.0002723089072243553, 'Z': 1.0}]

# Potential
a1, a2 = 1., 1.
potential = {}
potential['function'] = lambda Z1,Z2,r: log(1.+a1*exp(-Z1*Z2*r/a2))
potential['rMin'] = 0.0001 # first grid point
potential['rMax'] = 100. # last grid point
potential['nGrid'] = 1000 # number grid points
potential['gridType'] = "OPTIMIZED" # grid type (LINEAR, LOG, OPTIMIZED (Ilkka only!))

# Long-range breakup
breakup = {}
breakup['D'] = 3 # dimension
breakup['L'] = L # length of box
breakup['rMin'] = 0.0001 # first grid point
breakup['rMax'] = sqrt(breakup['D'])*breakup['L']/2. # last grid point
breakup['rCut'] = breakup['L']/2. # r cutoff for ewald
breakup['nGrid'] = 1000 # number of grid points
breakup['gridType'] = "OPTIMIZED" # grid type (LINEAR, LOG, OPTIMIZED (Ilkka only!))
breakup['nKnots'] = 10 # number of knots in spline (probably fine)
breakup['nImages'] = 10 # Naive check
breakup['GenPIMCPairActionFile'] = 0 # Whether or not to generate the *.PairAction file

# Pair action objects (object type, kspace cutoff, breakup type)
# type : 2 - dU/dBeta, 1 - U, 0 - V
# breakup : 2 - Short-ranged only, 1 - Optimized breakup, 0 - Classical Ewald breakup
objects = [{'type': 0, 'breakup': 2, 'kCut': 0.}]

# Exact location of PAGEN scripts
PAGEN_HOME = os.environ['HOME']+'/src/pagen'

#-----------------------------
# DO NOT EDIT BELOW THIS LINE
#-----------------------------
sys.path.append(PAGEN_HOME)
from GenPairAction import *
run(units,particles,potential,squarer,breakup,objects)
