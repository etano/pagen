from math import pi, sqrt
from numpy import loadtxt, unique
import subprocess
import FitPA
from FixTail import FixTail
from GenGrid import GenGrid
import h5py as h5
import os, sys

PAGEN_HOME=os.path.dirname(os.path.realpath(__file__))

def Square(particles,squarer):
    # Do nothing
    print 'No squaring required for bare action'

def Breakup(particles,potential,squarer,breakup,objects):
    paIndex = 0
    for i in xrange(0, len(particles)):
        for j in xrange(i+1, len(particles)):

            # Start from 1
            paIndex += 10

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
                if o['type'] == 0:
                    paPrefix = 'v'
                    cofactor = 1.
                    print '****\n', '****', paPrefix, '\n****'
    
                    # Copy file over
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
                if o['type'] == 0:
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
            f.flush()
            f.close()

            print '****************************************'

