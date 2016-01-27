from math import pi, sqrt
from numpy import loadtxt, unique
import subprocess
import FitPA
from FixTail import FixTail
from GenGrid import GenGrid
import Ewald
import h5py as h5
import os, sys

PAGEN_HOME=os.path.dirname(os.path.realpath(__file__))

def GenIlkkaSquarerInput(species_a, species_b, squarer):
    f = open('INPUT','w')
    f.write('LABEL\n')
    if (species_a['type'] == species_b['type']):
        f.write('%i %i\n' % (1, 2)) # of species, # of particles
        f.write('2 %f %f 1 64 5 0\n' % (species_a['Z'], 1./(2.*species_a['lambda']))) # of particles, charge (Z), mass, spin, Trotter number, multilevel L, 0 for bolzmannons and 1 for fermions (i.e. exchange)
    else:
        f.write('%i %i\n' % (2, 2)) # of species, # of particles
        f.write('1 %f %f 1 64 5 0\n' % (species_a['Z'], 1./(2.*species_a['lambda']))) # of particles, charge (Z), mass, spin, Trotter number, multilevel L, 0 for bolzmannons and 1 for fermions (i.e. exchange)
        f.write('1 %f %f 1 64 5 0\n' % (species_b['Z'], 1./(2.*species_b['lambda']))) # of particles, charge (Z), mass, spin, Trotter number, multilevel L, 0 for bolzmannons and 1 for fermions (i.e. exchange)
    f.write('%i\n' % (2)) # of particles to be moved (top to bottom)
    f.write('100\n') # frequency of calculating observables
    f.write('%f 0\n' % (squarer['tau'])) # tau, 0 or temperature, 1
    f.write('5000 100000\n') # of blocks, # of MC steps in a block
    f.write('Squaring %i %i 1\n' % (squarer['n_grid'],squarer['n_square'])) # text, grid points, # of squaring steps
    f.write('Box %f %f %f\n' % (squarer['r_max'],squarer['r_max'],squarer['r_max'])) # text, box dimensions in atomic units (Bohr radius)
    f.write('NumOfImages 0\n') # # of images
    f.write('DisplaceMove 0 0.0\n') # text, 0 if not used (1 if used), 0.1 would mean 10 percent of moves are displace moves
    f.write('Estimator 1\n') # 0=IK_type (atoms and if all quantum particles), 1=thermal estimator
    f.write('PairCorr 0.0 10.0 200\n') # grid_start, grid_end, number of grid points (if 0 then not calculated)
    f.write('LongRange 0 0.0\n')
    f.write('GetConf 0\n')
    f.write('LevelUpdate 1 0.40 0.55\n')
    f.close()

def Square(pa_object):
    # Perform squaring
    GenIlkkaSquarerInput(pa_object['species_a'],pa_object['species_b'],pa_object['squarer'])
    subprocess.call([PAGEN_HOME+'/ilkkaSquarer/ilkkaSquarer'])

    # Move files
    prefix = pa_object['species_a']['type']+'_'+pa_object['species_b']['type']
    subprocess.call(['cp','ud.1.txt',prefix+'_sq_u_diag.dat'])
    subprocess.call(['cp','us.1.txt',prefix+'_sq_u_offdiag.dat'])
    subprocess.call(['cp','dud.1.txt',prefix+'_sq_du_diag.dat'])
    subprocess.call(['cp','dus.1.txt',prefix+'_sq_du_offdiag.dat'])
    subprocess.call(['rm','ud.1.txt'])
    subprocess.call(['rm','us.1.txt'])
    subprocess.call(['rm','dud.1.txt'])
    subprocess.call(['rm','dus.1.txt'])

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
    for [object_type,object_index,cofactor] in [['du',2,1.],['u',1,tau],['v',0,1.]]:
        print '****\n', '****', object_type, '\n****'

        # Fix tail
        subprocess.call(['cp','-n',prefix+'_sq_'+object_type+'_diag.dat',prefix+'_sq_'+object_type+'_diag_orig.dat']) # Backup original
        subprocess.call(['cp',prefix+'_sq_'+object_type+'_diag_orig.dat',prefix+'_sq_'+object_type+'_diag.dat']) # Replace with original
        if object_type != 'v' and breakup['type'] != 'StandardEwald':
            n_points_to_fit = 10 # number of points to fit the tail with
            asymptote = lambda r: cofactor*potential['function'](Z1,Z2,r)
            try:
                FixTail(prefix+'_sq_'+object_type+'_diag.dat',n_points_to_fit,asymptote)
            except:
                break

        # Do breakup
        if breakup['type'] != 'None':
            if object_type == 'u':
                Ewald.run(breakup,object_type,prefix,Z1*Z2,tau)
            else:
                Ewald.run(breakup,object_type,prefix,Z1*Z2,1.)
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
    for [object_type,cofactor] in [['du',1.],['u',tau],['v',1]]:
        pa_group = f.create_group(object_type)

        # Write out diagonal PA
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

        # Write out off diagonal PA
        if object_type != 'v':
            pa_subgroup = pa_group.create_group('off_diag')
            pa_array = loadtxt(prefix+'_sq_'+object_type+'_offdiag.dat', comments='#')
            xs = unique(pa_array[:,0])
            ys = unique(pa_array[:,1])
            zs = pa_array[:,2].reshape((len(xs),len(ys)))
            pa_subgroup.create_dataset('x',data=xs)
            pa_subgroup.create_dataset('n_x',data=len(xs))
            pa_subgroup.create_dataset('y',data=ys)
            pa_subgroup.create_dataset('n_y',data=len(ys))
            pa_subgroup.create_dataset(object_type+'_xy',data=zs)

    f.flush()
    f.close()

    print '****************************************'

