from math import pi, sqrt
from numpy import loadtxt, unique
import subprocess
import FitPA
from FixTail import FixTail
from GenGrid import GenGrid
import h5py as h5
import os, sys

PAGEN_HOME=os.path.dirname(os.path.realpath(__file__))

def GenIlkkaSquarerInput(particles,squarer):
    N = len(particles)
    f = open('INPUT','w')
    f.write('LABEL\n')
    f.write('%i %i\n' % (N, 2*N)) # of species, # of particles
    for p in particles:
        f.write('2 %f %f 1 64 5 0\n' % (p['Z'], 1./(2.*p['lambda']))) # of particles, charge (Z), mass, spin, Trotter number, multilevel L, 0 for bolzmannons and 1 for fermions (i.e. exchange)
    f.write('%i\n' % (2*N)) # of particles to be moved (top to bottom)
    f.write('100\n') # frequency of calculating observables
    f.write('%f 0\n' % (squarer['tau'])) # tau, 0 or temperature, 1
    f.write('5000 100000\n') # of blocks, # of MC steps in a block
    f.write('Squaring %i %i 1\n' % (squarer['nGrid'],squarer['nSquare'])) # text, grid points, # of squaring steps
    f.write('Box %f %f %f\n' % (squarer['rMax'],squarer['rMax'],squarer['rMax'])) # text, box dimensions in atomic units (Bohr radius)
    f.write('NumOfImages 0\n') # # of images
    f.write('DisplaceMove 0 0.0\n') # text, 0 if not used (1 if used), 0.1 would mean 10 percent of moves are displace moves
    f.write('Estimator 1\n') # 0=IK_type (atoms and if all quantum particles), 1=thermal estimator
    f.write('PairCorr 0.0 10.0 200\n') # grid_start, grid_end, number of grid points (if 0 then not calculated)
    f.write('LongRange 0 0.0\n')
    f.write('GetConf 0\n')
    f.write('LevelUpdate 1 0.40 0.55\n')
    f.close()

def GenPairActionInput(prefix,type1,lam1,type2,lam2,D,longRange):
    f = open(prefix+'.PairAction','w')
    f.write('Section (Fits)')
    f.write('\n{')
    f.write('\n  string Type="IlkkaFit";')
    f.write('\n  int NumOffDiagonalTerms = -1;')
    f.write('\n  Section (Particle1)')
    f.write('\n  {')
    f.write('\n    string Name = "'+type1+'";')
    f.write('\n    double lambda ='+str(lam1)+';')
    f.write('\n    int Ndim = '+str(D)+';')
    f.write('\n  }')
    f.write('\n  Section (Particle2)')
    f.write('\n  {')
    f.write('\n    string Name = "'+type2+'";')
    f.write('\n    double lambda ='+str(lam2)+';')
    f.write('\n    int Ndim = '+str(D)+';')
    f.write('\n  }')
    f.write('\n  string Ilkkadmfile = "'+prefix+'.h5";')
    if longRange:
        f.write('\n  bool longRange = true;')
    f.write('\n}')
    f.close()

def Square(particles,squarer):
    # Perform squaring
    GenIlkkaSquarerInput(particles,squarer)
    subprocess.call([PAGEN_HOME+'/ilkkaSquarer/ilkkaSquarer'])

    # Check for duplicates
    paIndex = 0
    combos = []
    for i in xrange(0, len(particles)):
        for j in xrange(i, len(particles)):
            paIndex += 1
            [Z1Z2,lam1,lam2] = [particles[i]['Z']*particles[j]['Z'],particles[i]['lambda'],particles[j]['lambda']]
            try:
                prevIndex = combos.index([Z1Z2,abs(lam1-lam2)]) + 1
                subprocess.call(['cp','ud.'+str(prevIndex)+'.txt','ud.'+str(paIndex)+'.txt'])
                subprocess.call(['cp','us.'+str(prevIndex)+'.txt','us.'+str(paIndex)+'.txt'])
                subprocess.call(['cp','dud.'+str(prevIndex)+'.txt','dud.'+str(paIndex)+'.txt'])
                subprocess.call(['cp','dus.'+str(prevIndex)+'.txt','dus.'+str(paIndex)+'.txt'])
            except:
                combos.append([Z1Z2,abs(lam1-lam2)])

def Breakup(particles,potential,squarer,breakup,objects):
    paIndex = 0
    for i in xrange(0, len(particles)):
        for j in xrange(i, len(particles)):

            # Start from 1
            paIndex += 1

            # Assign particle attributes
            [type1, lam1, Z1] = particles[i]['type'], particles[i]['lambda'], particles[i]['Z']
            [type2, lam2, Z2] = particles[j]['type'], particles[j]['lambda'], particles[j]['Z']
            tau = squarer['tau']
            print '****************************************'
            print '****', type1, ', lam1 =', lam1, ', Z1 =', Z1
            print '****', type2, ', lam2 =', lam2, ', Z2 =', Z2
            print '**** tau =', tau

            # Write potential
            f = open('v.'+str(paIndex)+'.txt','w')
            rs = GenGrid(potential)
            for r in rs:
                f.write('%.10E %.10E\n' % (r,potential['function'](Z1,Z2,r)))
            f.close()

            # Write grid
            f = open('grid.'+str(paIndex)+'.txt','w')
            rs = GenGrid(breakup)
            for r in rs:
                f.write('%.10E\n' % (r))
            f.close()

            # Perform breakup
            for o in objects:
                if o['type'] == 2:
                    paPrefix = 'du'
                    cofactor = 1.
                elif o['type'] == 1:
                    paPrefix = 'u'
                    cofactor = tau
                elif o['type'] == 0:
                    paPrefix = 'v'
                    cofactor = 1.
                print '****\n', '****', paPrefix, '\n****'
    
                # Fix tail
                subprocess.call(['cp','-n',paPrefix+'d.'+str(paIndex)+'.txt',paPrefix+'d.'+str(paIndex)+'.orig.txt']) # Backup original
                subprocess.call(['cp',paPrefix+'d.'+str(paIndex)+'.orig.txt',paPrefix+'d.'+str(paIndex)+'.txt']) # Replace with original
                if o['type'] != 0 and o['breakup'] != 0:
                    nPointsToFit = 10 # number of points to fit the tail with
                    asymptote = lambda r: cofactor*potential['function'](Z1,Z2,r)
                    try:
                        FixTail(paPrefix+'d.'+str(paIndex)+'.txt',nPointsToFit,asymptote)
                    except:
                        break
                elif o['type'] == 0:
                    subprocess.call(['cp',paPrefix+'.'+str(paIndex)+'.txt',paPrefix+'d.'+str(paIndex)+'.txt'])

                # Do breakup
                if o['breakup'] != 2:
                    if breakup['gridType']=="LINEAR":
                        gridIndex = 0
                    elif breakup['gridType']=="LOG":
                        gridIndex = 1
                    elif breakup['gridType']=="OPTIMIZED":
                        gridIndex = 2
                    else:
                        print 'Unrecognized grid:', breakup['gridType']
                    subprocess.call([PAGEN_HOME+'/ewald/bin/ewald',str(breakup['L']),str(o['kCut']),str(breakup['rMin']),
                                             str(breakup['rCut']),str(breakup['nGrid']),str(gridIndex),
                                             str(Z1*Z2),str(o['breakup']),str(o['type']),str(paIndex),
                                             str(breakup['nKnots']),str(tau),str(breakup['nImages'])])
                else:
                    subprocess.call(['cp',paPrefix+'d.'+str(paIndex)+'.txt',paPrefix+'d.'+str(paIndex)+'.r.txt'])
    
                # Fit off-diagonal
                if o['type'] != 0 and squarer['nOrder'] > 0:
                    FitPA.main(['',squarer['nOrder'],paPrefix+'s.'+str(paIndex)+'.txt',0])
                elif o['type'] == 0 and o['breakup'] != 2:
                    subprocess.call(['cp',paPrefix+'.'+str(paIndex)+'.r.txt',paPrefix+'d.'+str(paIndex)+'.r.txt'])
                    subprocess.call(['cp',paPrefix+'.'+str(paIndex)+'.k.txt',paPrefix+'d.'+str(paIndex)+'.k.txt'])
    
            # Write to h5 file
            prefix = type1+'-'+type2
            f = h5.File(prefix+'.h5','w')
            info = f.create_group('Info')
            info.attrs.create('type1',type1)
            info.attrs.create('type2',type2)
            info.attrs.create('lam1',lam1)
            info.attrs.create('lam2',lam2)
            info.attrs.create('Z1',Z1)
            info.attrs.create('Z2',Z2)
            info.create_dataset('Z1',data=[Z1])
            info.create_dataset('Z2',data=[Z2])
            info.create_dataset('tau',data=[tau])
            for o in objects:
                if o['type'] == 2:
                    paPrefix = 'du'
                elif o['type'] == 1:
                    paPrefix = 'u'
                elif o['type'] == 0:
                    paPrefix = 'v'
                paGroup = f.create_group(paPrefix)

                # Write out diagonal PA
                paSubgroup = paGroup.create_group('diag')
                paArray = loadtxt(paPrefix+'d.'+str(paIndex)+'.txt', comments='#')
                paSubgroup.create_dataset('r',data=paArray[:,0])
                paSubgroup.create_dataset('nr',data=len(paArray[:,0]))
                paSubgroup.create_dataset(paPrefix+'_r',data=paArray[:,1])
                if o['breakup'] != 2:
                    paArray = loadtxt(paPrefix+'d.'+str(paIndex)+'.r.txt', comments='#')
                    paSubgroup.create_dataset(paPrefix+'Long_r0',data=[paArray[0,1]])
                    paSubgroup.create_dataset('rLong',data=paArray[1:,0])
                    paSubgroup.create_dataset('nrLong',data=len(paArray[1:,0]))
                    paSubgroup.create_dataset(paPrefix+'Long_r',data=paArray[1:,1])
                    paArray = loadtxt(paPrefix+'d.'+str(paIndex)+'.k.txt', comments='#')
                    paSubgroup.create_dataset(paPrefix+'Long_k0',data=[paArray[0,1]])
                    paSubgroup.create_dataset('k',data=paArray[:,0])
                    paSubgroup.create_dataset('nk',data=len(paArray[:,0]))
                    paSubgroup.create_dataset(paPrefix+'Long_k',data=paArray[:,1])
    
                # Write out off diagonal PA
                if o['type'] != 0:
                    paSubgroup = paGroup.create_group('offDiag')
                    paArray = loadtxt(paPrefix+'s.'+str(paIndex)+'.txt', comments='#')
                    xs = unique(paArray[:,0])
                    ys = unique(paArray[:,1])
                    zs = paArray[:,2].reshape((len(xs),len(ys)))
                    paSubgroup.create_dataset('x',data=xs)
                    paSubgroup.create_dataset('nx',data=len(xs))
                    paSubgroup.create_dataset('y',data=ys)
                    paSubgroup.create_dataset('ny',data=len(ys))
                    paSubgroup.create_dataset(paPrefix+'_xy',data=zs)

                    # Write out fit to off diagonal PA if desired
                    for iOrder in range(1,squarer['nOrder']+1):
                        paSubSubgroup = paSubgroup.create_group('A.'+str(iOrder))
                        paArray = loadtxt(paPrefix+'s.'+str(paIndex)+'.A.'+str(iOrder)+'.txt', comments='#')
                        paSubSubgroup.create_dataset('r',data=paArray[:,0])
                        paSubSubgroup.create_dataset('nr',data=len(paArray[:,0]))
                        paSubSubgroup.create_dataset('A',data=paArray[:,1])
    
            f.flush()
            f.close()
    
            # Build PairAction File
            if breakup['GenPIMCPairActionFile']:
                print 'Creating pairAction file...'
                longRange = 0
                if objects[0]['breakup'] != 2:
                    longRange = 1
                GenPairActionInput(prefix,type1,lam1,type2,lam2,squarer['D'],longRange)

            print '****************************************'

