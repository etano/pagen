import sys
from math import *
from numpy import *
import subprocess

def GetUnique(a):
  seen = set()
  return [x for x in a if str(x) not in seen and not seen.add(str(x))]

def GenIlkkaSquarerInput(D,tau,L,particles,nGrid,nSquare):
    N = len(particles)
    f = open('INPUT','w')
    f.write('LABEL\n')
    f.write('%i %i\n' % (N, 2*N)) # of species, # of particles
    for [type,lam,Z] in particles:
        f.write('2 %f %f 1 64 5 0\n' % (Z, 1./(2.*lam))) # of particles, charge (Z), mass, spin, Trotter #, multilevel L, 0 for bolzmannons and 1 for fermions (i.e. exchange)
    f.write('%i\n' % (2*N)) # of particles to be moved (top to bottom)
    f.write('100\n') # frequency of calculating observables
    f.write('%f 0\n' % (tau)) # tau, 0 or temperature, 1
    f.write('5000 100000\n') # of blocks, # of MC steps in a block
    f.write('Squaring %i %i 1\n' % (nGrid,nSquare)) # text, grid points, # of squaring steps
    f.write('Box %f %f %f\n' % (L,L,L)) # text, box dimensions in atomic units (Bohr radius)
    f.write('NumOfImages 0\n') # # of images
    f.write('DisplaceMove 0 0.0\n') # text, 0 if not used (1 if used), 0.1 would mean 10 percent of moves are displace moves
    f.write('Estimator 1\n') # 0=IK_type (atoms and if all quantum particles), 1=thermal estimator
    f.write('PairCorr 0.0 10.0 200\n') # grid_start, grid_end, # of grid points (if 0 then not calculated)
    f.write('LongRange 0 0.0\n')
    f.write('GetConf 0\n')
    f.write('LevelUpdate 1 0.40 0.55\n')
    f.close()

# Input parameters
particles = [['e',0.5,1.0],['p',0.0002723089072243553,-1.0]]

## PIMC simulation
D = 3 # dimension
T = 0.025 # desired temperature of PIMC simulation
tau = 0.125 # desired timestep of PIMC simulation
M = int((1/T)/tau) # number of time slices
tau = (1/T)/M # used tau

## Squarer
L = 100.0 # length of box
nGrid = 200 # number of grid points
nSquare = 33 # total number of squarings to reach lowest temperature

# Create Ilkka Squarer input
print '**** Creating Squarer Inputs ****'
GenIlkkaSquarerInput(D,tau,L,particles,nGrid,nSquare)
subprocess.call(['ilkkaSquarer'])

## Long-range breakup
L = 10.0 # length of box
r0 = 0.01 # first grid point
rMax = L/2. # last grid point
nGrid = 400 # number of grid points
nMax = 10 # index of k cutoff for ewald
rCut = 5.0 # r cutoff for ewald
breakupType = 1 # 2 - Short-ranged only, 1 - Optimized breakup, 0 - Classical Ewald breakup
breakupObject = 1 # 2 - dU/dBeta, 1 - U, 0 - V
gridType = "LINEAR" # LOG/LINEAR
nKnots = 14
nImages = 10 # Naive Check

paIndex = 0
for i in xrange(0, len(particles)):
    for j in xrange(i, len(particles)):
        paIndex += 1
        [type1, lam1, Z1] = particles[i]
        [type2, lam2, Z2] = particles[j]
        print type1, lam1, Z1, type2, lam2, Z2

        # Write potential
        print '**** Writing potential to file ****'
        f = open('v.'+str(paIndex)+'.txt','w')
        if gridType=="LINEAR":
          gridIndex = 0
          rs = linspace(r0,20*rMax,num=20*nGrid,endpoint=True)
        elif gridType=="LOG":
          gridIndex = 1
          rs = logspace(log10(r0),log10(20*rMax),num=20*nGrid,endpoint=True)
        else:
          print 'Unrecognized grid'
        for r in rs:
          f.write('%.10E %.10E\n' % (r,Z1*Z2/r))
        f.close()

        # Perform breakup
        print '**** Performing breakup ****'
        if breakupType != 2:
          subprocess.call(['ewald',str(L),str(nMax),str(r0),str(rCut),str(nGrid),str(gridIndex),str(Z1*Z2),str(breakupType),str(breakupObject),str(paIndex),str(nKnots),str(tau),str(nImages)])
          rMax = 0.75*sqrt(3)*L
          kCut = 2*pi*nMax/L
        else:
          f = open('rData.txt','w')
          if gridType=="LINEAR":
            rs = linspace(r0,L/2.,num=nGrid,endpoint=True)
          elif gridType=="LOG":
            rs = logspace(log10(r0),log10(L/2.),num=nGrid,endpoint=True)
          else:
            print 'Unrecognized grid'
          for r in rs:
            f.write('%.10E %.10E %.10E\n' % (r,Z1*Z2/r,0.0))
          f.close()
          rMax = L/2.
