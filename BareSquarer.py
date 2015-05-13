from math import pi, sqrt
from numpy import loadtxt, unique
import subprocess
import FitPA
from FixTail import FixTail
from GenGrid import GenGrid
import h5py as h5
import os, sys

PAGEN_HOME=os.path.dirname(os.path.realpath(__file__))

def Square(pa_object):
    # Do nothing
    print 'No squaring required for bare action'

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

    # Write potential
    prefix = species_a['type']+'_'+species_b['type']
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

    # Perform breakup
    [object_type,object_index,cofactor] = ['v',0,1.]
    print '****\n', '****', object_type, '\n****'
    if breakup['type'] != 'None':
        if breakup['grid_type']=="LINEAR":
            grid_index = 0
        elif breakup['grid_type']=="LOG":
            grid_index = 1
        elif breakup['grid_type']=="OPTIMIZED":
            grid_index = 2
        else:
            print 'Unrecognized grid:', breakup['grid_type']
        subprocess.call([PAGEN_HOME+'/ewald/bin/ewald',str(breakup['L']),str(breakup['k_cut']),str(breakup['r_min']),str(breakup['r_cut']),str(breakup['n_grid']),str(grid_index),str(Z1*Z2),str(breakup_index),str(object_index),prefix,str(breakup['n_knots']),str(tau),str(breakup['n_images'])])
    else:
        subprocess.call(['cp',prefix+'_sq_'+object_type+'_diag.dat',prefix+'_sq_'+object_type+'_diag_r.dat'])

    # Write to h5 file
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

    # Write out diagonal PA
    pa_group = f.create_group(object_type)
    pa_subgroup = pa_group.create_group('diag')
    pa_array = loadtxt(prefix+'_sq_'+object_type+'_diag.dat', comments='#')
    pa_subgroup.create_dataset('r',data=pa_array[:,0])
    pa_subgroup.create_dataset('n_r',data=len(pa_array[:,0]))
    pa_subgroup.create_dataset(object_type+'_r',data=pa_array[:,1])
    if breakup['type'] != 'None':
        pa_array = loadtxt(prefix+'_sq_'+object_type+'_diag_r.dat', comments='#')
        pa_subgroup.create_dataset(object_type+'_long_r_0',data=[pa_array[0,1]])
        pa_subgroup.create_dataset('r_long',data=pa_array[1:,0])
        pa_subgroup.create_dataset('n_r_long',data=len(pa_array[1:,0]))
        pa_subgroup.create_dataset(object_type+'_long_r',data=pa_array[1:,1])
        pa_array = loadtxt(prefix+'_sq_'+object_type+'_diag_k.dat', comments='#')
        pa_subgroup.create_dataset(object_type+'_long_k_0',data=[pa_array[0,1]])
        pa_subgroup.create_dataset('k',data=pa_array[:,0])
        pa_subgroup.create_dataset('n_k',data=len(pa_array[:,0]))
        pa_subgroup.create_dataset(object_type+'_long_k',data=pa_array[:,1])
    f.flush()
    f.close()

    print '****************************************'

