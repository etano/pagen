import sys, os
import subprocess
from math import sqrt, pi
from numpy import loadtxt
import DavidParse
from GenGrid import GenGrid

PAGEN_HOME=os.path.dirname(os.path.realpath(__file__))

def GetUnique(a):
  seen = set()
  return [x for x in a if str(x) not in seen and not seen.add(str(x))]

def GenPotgenInput(prefix,type1,type2,lam1,lam2,Z1Z2,L,D,tau,grid_type,n_grid,r_min,r_max,r_cut,r_paste,k_cut,n_temp,n_square,breakup):
    # Determine temperatures
    if n_temp > 8:
        print 'WARNING: Max n_temp is 8!'
        n_temp = 8
    min_tau = tau
    max_tau = min_tau*(2**(n_temp-1))
    if n_square < n_temp:
        n_square = n_temp

    print 'Creating '+prefix+'.in'
    f = open(prefix+'.in','w')
    f.write(' UNITS H A')
    f.write('\n TYPE '+type1+' %f' % (lam1))
    f.write('\n TYPE '+type2+' %f' % (lam2))
    f.write('\n GRID %i %s %f %f %f' % (n_grid,grid_type,r_min,r_max,r_paste))
    f.write('\n SQUARER %f %i %i 3 30 %i' % (1./max_tau,n_temp,D,n_square))
    box_string = ' '.join([str(L) for i in range(D)])
    if breakup == 2:
        f.write('\n POT COUL %f %f 0.D0' % (r_cut,Z1Z2))
    if breakup == 1:
        f.write('\n POT COULOPT %f %f %f %i 0.D0 1.0 %s' % (r_cut,k_cut,Z1Z2,D,box_string))
    elif breakup == 0:
        f.write('\n POT COULLR %f %f %f %i 0.D0 1.0 %s' % (r_cut,k_cut,Z1Z2,D,box_string))
    f.close()

def GenPairActionInput(prefix,type1,lam1,type2,lam2,D,long_range):
    f = open(prefix+'.PairAction','w')
    f.write('Section (Fits)')
    f.write('\n{')
    f.write('\n  string Type="DavidFit";')
    f.write('\n  int NumOffDiagonalTerms = 3;')
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
    f.write('\n  string Daviddmfile = "'+prefix+'.h5";')
    if long_range:
        f.write('\n  bool long_range = true;')
    f.write('\n}')
    f.close()

def Square(pa_object):
    species_a = pa_object['species_a']
    species_b = pa_object['species_b']
    squarer = pa_object['squarer']
    breakup = pa_object['breakup']
    potential = pa_object['potential']

    # Prefix
    prefix = species_a['type']+'_'+species_b['type']

    # Squarer
    print 'Performing squaring procedure...'
    subprocess.call([PAGEN_HOME+'/davidSquarer/sqdir/squarer',prefix+'_sq'])

    # Density Matrix Parser
    print 'Parsing density matrix...'
    DavidParse.main(['',prefix+'_sq.dm'])

    ## Build PairAction File
    #print 'Creating pairAction file...'
    #long_range = 0
    #if breakup['type'] != 'None':
    #    long_range = 1
    #GenPairActionInput(prefix,type1,lam1,type2,lam2,squarer['n_d'],long_range)

def Breakup(pa_object):
    species_a = pa_object['species_a']
    species_b = pa_object['species_b']
    squarer = pa_object['squarer']
    breakup = pa_object['breakup']
    potential = pa_object['potential']

    # Assign breakup index
    if breakup['type'] == 'StandardEwald':
        breakup_index = 0
    elif breakup['type'] == 'OptimizedEwald':
        breakup_index = 1
    elif breakup['type'] == 'None':
        breakup_index = 2

    # Assign particle attributes
    [type1, lam1, Z1] = species_a['type'], species_a['lambda'], species_a['Z']
    [type2, lam2, Z2] = species_b['type'], species_b['lambda'], species_b['Z']
    tau = squarer['tau']
    print '****************************************'
    print '****', type1, ', lam1 =', lam1, ', Z1 =', Z1
    print '****', type2, ', lam2 =', lam2, ', Z2 =', Z2
    print '**** tau =', tau

    # Generate potgen input
    prefix = species_a['type']+'_'+species_b['type']
    GenPotgenInput(prefix+'_sq',type1,type2,lam1,lam2,Z1*Z2,
                   breakup['L'],breakup['n_d'],squarer['tau'],breakup['grid_type'],
                   breakup['n_grid'],breakup['r_min'],breakup['r_max'],breakup['r_cut'],
                   breakup['r_paste'],breakup['k_cut'],
                   squarer['n_temp'],squarer['n_square'],breakup_index)

    # Write potential
    f = open(prefix+'_sq_v_diag.dat','w')
    rs = GenGrid(potential)
    for r in rs:
        f.write('%.10E %.10E\n' % (r,potential['function'](Z1,Z2,r)))
    f.close()

    # Write grid
    f = open(prefix+'_sq_grid.dat','w')
    rs = GenGrid(breakup)
    for r in rs:
        f.write('%.10E\n' % (r))
    f.close()

    # Object is only potential
    [object_type,object_index,cofactor] = ['v',0,1.]

    # Do breakup
    if breakup['type'] != 'None':
        subprocess.call([PAGEN_HOME+'/davidSquarer/sqdir/potgen_lr',prefix+'_sq'])
    else:
        subprocess.call([PAGEN_HOME+'/davidSquarer/sqdir/potgen_sr',prefix+'_sq'])
