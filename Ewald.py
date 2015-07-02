import sys, os
from math import pi, ceil, floor, erf, cos, sin
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline
from scipy.integrate import romberg as integrate
from GenGrid import *

class LPQHIBasis(object):
    """Docstring for Basis """

    def __init__(self,L,r_c,n_d,n_knots):
        self.L = L
        self.r_c = r_c
        self.n_d = n_d
        self.n_knots = n_knots
        self.S = np.array([np.array([   1.,   0.,   0., -10.,  15.,  -6.]),
                           np.array([   0.,   1.,   0.,  -6.,   8.,  -3.]),
                           np.array([   0.,   0.,  0.5, -1.5,  1.5, -0.5])])

        # Set up box
        self.SetUpBox()

        # Set delta
        self.SetDelta()

    def SetUpBox(self):
        self.box,self.i_box = [],[]
        self.vol = 1.
        for d_i in range(self.n_d):
            self.box.append(self.L) # TODO: Currently only cubic box
            self.i_box.append(1./self.box[d_i])
            self.vol *= self.box[d_i]
        self.box = np.array(self.box)
        self.prefactor = 4.*pi/self.vol

    def SetNKnots(self,n_knots):
        self.n_knots = n_knots
        if (self.r_c != 0.):
            self.SetDelta()

    def SetRC(self,r_c):
        self.r_c = r_c
        if (self.n_knots != 0):
            self.SetDelta()

    def SetDelta(self):
        self.delta = self.r_c/(self.n_knots-1.)
        self.i_delta = 1./self.delta

    def GetNElements(self):
        return 3*self.n_knots # TODO: Why 3?

    def EPlus(self,i,k,n):
        eye = complex(0.,1.)
        if (n == 0):
            e1 = complex(cos(k*self.delta)-1., sin(k*self.delta))
            e2 = complex(cos(k*self.delta*i), sin(k*self.delta*i))
            return eye*(-(e1*e2/k))
        else:
            t1 = complex(cos(k*self.delta*(i+1)), sin(k*self.delta*(i+1)))
            t2 = -(n/self.delta)*self.EPlus(i,k,n-1)
            return eye*((-(t1+t2)/k))

    def EMinus(self,i,k,n):
        eye = complex(0.,1.)
        if (n == 0):
            e1 = complex(cos(k*self.delta)-1., -sin(k*self.delta))
            e2 = complex(cos(k*self.delta*i), sin(k*self.delta*i))
            return eye*(-(e1*e2/k))
        else:
            sign = -1. if (n & 1) else 1.
            t1 = sign*complex(cos(k*self.delta*(i-1)), sin(k*self.delta*(i-1)))
            t2 = -(n/self.delta)*self.EMinus(i,k,n-1)
            return eye*((-(t1+t2)/k))

    def DPlus(self,i,k,n):
        z1 = self.EPlus(i,k,n+1)
        z2 = self.EPlus(i,k,n)
        return (self.prefactor/k)*(self.delta*z1.imag + i*self.delta*z2.imag)

    def DMinus(self,i,k,n):
        z1 = self.EMinus(i,k,n+1)
        z2 = self.EMinus(i,k,n)
        return (-self.prefactor/k)*(self.delta*z1.imag + i*self.delta*z2.imag)

    def h(self,n,r):
        i = n/3
        alpha = n - 3*i
        ra = self.delta*(i-1)
        rb = self.delta*i
        rc = self.delta*(i+1)
        rc = min(self.r_c, rc)
        if (r > ra) and (r <= rb):
            sum = 0.
            prod = 1.
            for j in range(0,6):
                sum += self.S[alpha,j]*prod
                prod *= (rb-r)*self.i_delta
            for j in range(0,alpha):
                sum *= -1.
            return sum
        elif (r > rb) and (r <= rc):
            sum = 0.
            prod = 1.
            for j in range(0,6):
                sum += self.S[alpha,j]*prod
                prod *= (r-rb)*self.i_delta
            return sum
        return 0.

    def c(self,m,k):
        i = m/3
        alpha = m - 3*i
        sum = 0.
        if (i == 0):
            for n in range(0,6):
                sum += self.S[alpha,n]*self.DPlus(i,k,n)
        elif (i == self.n_knots-1):
            for n in range(0,6):
                sign = -1. if (alpha+n)&1 else 1.
                sum += self.S[alpha,n]*self.DMinus(i,k,n)*sign
        else:
            for n in range(0,6):
                sign = -1. if (alpha+n)&1 else 1.
                sum += self.S[alpha,n]*(self.DPlus(i,k,n) + self.DMinus(i,k,n)*sign)
        return sum

class EwaldBreakup(object):
    """Docstring for EwaldBreakup """

    def __init__(self,n_d,L,breakup_type,object_string,prefix,grid_type,r_min,r_max,n_points,r_cut,k_cut,z_1_z_2,cofactor,n_knots):
        """@todo: to be defined1 """
        self.n_d = n_d
        self.L = L
        self.breakup_type = breakup_type
        self.object_string = object_string
        self.prefix = prefix
        self.grid_type = grid_type
        self.r_min = r_min
        self.r_max = r_max
        self.n_points = n_points
        self.r_cut = r_cut
        self.k_cut = k_cut
        self.z_1_z_2 = z_1_z_2
        self.cofactor = cofactor
        self.n_knots = n_knots

        # Set up box
        self.SetUpBox()

        # Set up k vectors
        self.SetUpKs()

        # Set up grid
        self.SetUpGrid()

    def SetUpBox(self):
        self.box,self.i_box = [],[]
        self.vol = 1.
        for d_i in range(self.n_d):
            self.box.append(self.L)
            self.i_box.append(1./self.L)
            self.vol *= self.L # TODO: Currently only cubic box
        self.box = np.array(self.box)

    def Include(self,k):
        mag_k = np.sqrt(np.dot(k,k))
        if (mag_k <= self.k_cut):
            if (abs(k[0]) > 0.):
                return 1
            elif ((self.n_d > 1) and (k[0] == 0.) and (abs(k[1]) > 0.)):
                return 1
            elif ((self.n_d > 2) and (k[0] == 0.) and (k[1] == 0.) and (abs(k[2]) > 0.)):
                return 1
            else:
                return 0
        else:
            return 0

    def GenKis(self,max_k_i):
        k_is = []
        for n_x in range(-max_k_i[0],max_k_i[0]+1):
            if self.n_d == 1:
                k_is.append(np.array([n_x]))
            else:
                for n_y in range(-max_k_i[1],max_k_i[1]+1):
                    if self.n_d == 2:
                        k_is.append(np.array([n_x,n_y]))
                    else:
                        for n_z in range(-max_k_i[2],max_k_i[2]+1):
                            k_is.append(np.array([n_x,n_y,n_z]))
        return k_is

    def SetUpKs(self):
        # Set up k box/vol
        self.k_box = []
        self.k_vol = 1.
        for d_i in range(self.n_d):
            self.k_box.append(2.*pi/self.box[d_i])
            self.k_vol *= self.k_box[d_i]
        self.k_avg = self.k_vol**(1./self.n_d)

        # Set up max k indices
        max_k_i = []
        for d_i in range(self.n_d):
            max_k_i.append(int(ceil(1.1*self.k_cut/self.k_box[d_i])))

        # Generate k indices
        k_is = self.GenKis(max_k_i)

        # Set up ks
        self.ks = []
        self.mag_ks = []
        for k_i in k_is:
            k = k_i*self.k_box
            if self.Include(k):
                self.ks.append(k)
                self.mag_ks.append(np.sqrt(np.dot(k,k)))

    def Addk(self,mag_ks,mag_k,degeneracy=1):
        try:
            mag_ks[str(mag_k)] += degeneracy
        except:
            mag_ks[str(mag_k)] = degeneracy

    def ExtendKs(self,k_cont,k_max):
        # Set up max k indices
        max_k_i = []
        for d_i in range(self.n_d):
            max_k_i.append(int(ceil(1.1*k_cont/self.k_box[d_i])))

        # Generate k indices
        print '......generating ks...'
        k_is = self.GenKis(max_k_i)

        # Set up discrete ks
        print '......adding discrete ks...'
        mag_ks_dict = {}
        for k_i in k_is:
            k = k_i*self.k_box
            mag_k = np.sqrt(np.dot(k,k))
            if (mag_k > self.k_cut) and (mag_k < k_cont):
                self.Addk(mag_ks_dict,mag_k)

        # Set up continuous ks
        print '......adding continuous ks...'
        N = 4000 # TODO: This is just fixed
        delta_k = (k_max-k_cont)/N
        for i in range(N):
            k1 = k_cont + delta_k*i
            k2 = k1 + delta_k
            k = 0.5*(k1+k2)
            vol = (4.*pi/3.)*(k2*k2*k2-k1*k1*k1) # FIXME: Only 3D!
            degeneracy = vol/self.k_vol
            self.Addk(mag_ks_dict,k,degeneracy)

        # Create opt_mag_ks
        print '......creating ks object...'
        self.opt_mag_ks = []
        for mag_k in mag_ks_dict.iterkeys():
            self.opt_mag_ks.append([float(mag_k),mag_ks_dict[mag_k]])

    def SetUpGrid(self):
        self.rs = GenGrid({'grid_type':self.grid_type,
                           'r_min':self.r_min,
                           'r_max':self.r_max,
                           'n_grid':self.n_points})
        self.r_min = self.rs[0]
        self.r_max = self.rs[-1]
        self.n_points = len(self.rs)

    def DoBreakup(self):
        if self.breakup_type == 'OptimizedEwald':
            self.OptimizedBreakup()
        elif self.breakup_type == 'StandardEwald':
            self.StandardBreakup()
        else:
            print 'ERROR: Unrecognized breakup type.'
            sys.exit(2)

    def StandardBreakup(self):
        print 'Performing standard Ewald breakup...'

        # Set alpha TODO: Fixed for now
        alpha = np.sqrt(self.k_cut/(2.*self.r_cut))

        # Long range r space potential part
        f = open(self.prefix+'_sq_'+self.object_string+'_diag_r.dat','w')
        v_l_0 = 2.*self.cofactor*self.z_1_z_2*alpha/np.sqrt(pi)
        f.write('%.10E %.10E\n'%(0.,v_l_0))
        for i in range(self.n_points):
            r = self.rs[i]
            f.write('%.10E %.10E\n'%(r,self.cofactor*self.z_1_z_2*erf(alpha*r)/r))
        f.close()

        # Long range k space potential part
        f = open(self.prefix+'_sq_'+self.object_string+'_diag_k.dat','w')
        if (self.n_d == 2): # FIXME: Check this is correct
            f_v_s_0 = -2.*self.cofactor*np.sqrt(pi)*self.z_1_z_2/(alpha*self.vol)
        elif (self.n_d == 3):
            f_v_s_0 = -4.*self.cofactor*pi*self.z_1_z_2/(4.*alpha*alpha*self.vol)
        f.write('%.10E %.10E\n'%(0.,f_v_s_0))
        f_v_ls = []
        for k in self.ks:
            k_2 = np.dot(k,k)
            mag_k = np.sqrt(k_2)
            f_v_l = 0.
            if (self.n_d == 2): #FIXME: Check this is correct
                f_v_l = (2.*self.cofactor*np.sqrt(pi)*self.z_1_z_2/(mag_k*self.vol))*exp(-k_2/(4.*alpha*alpha))
            elif (self.n_d == 3):
                f_v_l = (4.*pi*self.cofactor*self.z_1_z_2/(k_2*self.vol))*exp(-k_2/(4.*alpha*alpha))
            f_v_ls.append([mag_k, f_v_l])
        mag_k_prev = -1
        for [mag_k,f_v_l] in sorted(f_v_ls):
            if abs(mag_k-mag_k_prev) > 1.e-8:
                f.write('%.10E %.10E\n'%(mag_k,f_v_l))
            mag_k_prev = mag_k
        f.close()

    def CalcXkCoul(self,k,r):
        x_k = 0.
        if (self.n_d == 2): # FIXME: This probably isn't right
            return -self.cofactor*(2.*pi*self.z_1_z_2/k)*cos(k*r)
        elif (self.n_d == 3):
            return -self.cofactor*(4.*pi*self.z_1_z_2/(k*k))*cos(k*r)

    def CalcXk(self, v_r_spline, r_cut, k, r_max):
        # Tolerances
        abs_tol = 1.e-11
        rel_tol = 1.e-11

        # Integrand
        v_integrand = lambda r: r*sin(k*r)*v_r_spline(r)

        # Calculate x_k
        r_first = r_cut + ((pi/k)-(r_cut % (pi/k)))
        if (self.n_d == 2): # FIXME: This probably isn't right
            prefactor = -2.*pi
        elif (self.n_d == 3):
            prefactor = -4.*pi/k
        x_k = 0.
        if (r_max >= r_first):
            # First segment
            x_k += prefactor * integrate(v_integrand, r_cut, r_first, divmax=20)

            # Other segments
            if (int(k) != 0):
                n_pi_k = int(k)*pi/k # TODO: Manually fixed number currently
            else:
                n_pi_k = pi/k
            n_seg = max(int(floor((r_max-r_first)/n_pi_k)),1)
            for i in range(n_seg):
                x_k += prefactor * integrate(v_integrand, r_first+i*n_pi_k, r_first+(i+1)*n_pi_k, divmax=20)
            r_end = r_first + n_seg*n_pi_k
        elif (r_max >= r_cut):
            x_k += prefactor * integrate(v_integrand, r_cut, r_max, divmax=20)
            r_end = r_max
        else:
            r_end = r_cut

        # Add in analytic part after r_end
        x_k += self.CalcXkCoul(k,r_end)

        return x_k

    def OptimizedBreakup(self):
        print 'Performing optimized Ewald breakup...'

        # Read potential and spline it
        print '...reading potential...'
        data = np.loadtxt(self.prefix+'_sq_'+self.object_string+'_diag.dat')
        rs,r_min,r_max = data[:,0],data[0,0],data[-1,0]
        vs = data[:,1]
        v_r_spline = spline(rs, vs)

        # Set up basis
        print '...forming LPQHI basis...'
        basis = LPQHIBasis(self.L,self.r_cut,self.n_d,self.n_knots)

        # Extend k space to continuum
        print '...extending k space...'
        k_cont = 50.*self.k_avg
        k_max = 50.*pi/basis.delta
        self.ExtendKs(k_cont,k_max)
        n_k = len(self.opt_mag_ks)
        n_r = basis.GetNElements()

        # Determine r max from tolerances
        v_tol = 1.e-5 # TODO: This is fixed
        v_c = self.cofactor*self.z_1_z_2/r_max
        if (abs(v_r_spline(r_max) - v_c) > v_tol):
            print 'WARNING: |v(r_max) - v_{c}(r_max)| = ',abs(v_r_spline(r_max) - v_c),'>',v_tol,'with r_max =',r_max,' v(r_max) = ',v_r_spline(r_max),' v_{c}(r_max) = ',v_c
        else:
            i = 1
            r_max = rs[i]
            while (abs(v_r_spline(r_max) - self.cofactor*self.z_1_z_2/r_max) > v_tol) and (i+1 < len(rs)-1):
                i += 1
                r_max = rs[i]
        if (r_max < self.r_cut):
            r_max = self.r_cut
        v_c = self.cofactor*self.z_1_z_2/r_max
        print '...setting r_max = ',r_max,'with |v(r_max) - v_{c}(r_max)| <',v_tol,'...'

        # Calculate x_k
        print '...calculating Xk...'
        x_k = []
        tot_x_k = 0.
        percent_i = 1
        for k_i in range(n_k):
            k = self.opt_mag_ks[k_i][0]
            x_k.append(self.CalcXk(v_r_spline, self.r_cut, k, r_max)/self.vol)
            if (float(k_i)/float(n_k)) > percent_i*0.1:
                print '......', percent_i*10, '% complete...'
                percent_i += 1
            tot_x_k += x_k[k_i]

        # Fill in c_n_k
        print '...filling in c_n_k...'
        c_n_k = np.zeros((n_r,n_k))
        for n in range(n_r):
            for k_i in range(n_k):
                c_n_k[n,k_i] = basis.c(n,self.opt_mag_ks[k_i][0])

        # Fill in A and b
        print '...filling in A and b...'
        A = np.zeros((n_r,n_r))
        b = np.zeros((n_r))
        for l in range(n_r):
            for k_i in range(n_k):
                b[l] += self.opt_mag_ks[k_i][1] * x_k[k_i] * c_n_k[l,k_i]
                for n in range(n_r):
                    A[l,n] += self.opt_mag_ks[k_i][1] * c_n_k[l,k_i] * c_n_k[n,k_i]

        # Add constraints
        t = np.zeros((n_r))
        adjust = np.ones((n_r)) # TODO: Currently no constraints

        # Reduce for constraints
        n_r_c = n_r
        for i in range(n_r):
            if not adjust[i]:
                n_r_c -= 1

        # Build constrained A_c and b_c
        A_c = np.zeros((n_r_c,n_r_c))
        b_c = np.zeros((n_r_c))
        j = 0
        for col in range(n_r):
            if adjust[col]:
                i = 0
                for row in range(n_r):
                    if adjust[row]:
                        A_c[i,j] = A[row,col]
                        i += 1
                j += 1
            else:
                for row in range(n_r):
                    b[row] -= A[row,col]*t[col]
        j = 0
        for row in range(n_r):
            if adjust[row]:
                b_c[j] = b[row]
            j += 1

        # Do SVD
        print '...performing SVD...'
        U, S, V = np.linalg.svd(A_c, full_matrices=True)

        # Get maximum value in S
        s_max = S[0]
        for i in range(1,n_r_c):
            s_max = max(S[i],s_max)

        # Check for negative singular values
        for i in range(n_r_c):
            if S[i] < 0.:
                print 'WARNING: Negative singular value.'

        # Assign inverse S
        breakup_tol = 1.e-16
        i_S = np.zeros((n_r_c))
        n_singular = 0
        for i in range(n_r_c):
            if (S[i] < breakup_tol*s_max):
                i_S[i] = 0.
            else:
                i_S[i] = 1./S[i]
            if (i_S[i] == 0.):
                n_singular += 1
        if (n_singular > 0):
            print 'WARNING: There were',n_singular,'singular values.'

        # Compute t_n, removing singular values
        t_c = np.zeros((n_r_c))
        for i in range(n_r_c):
            coef = 0.
            for j in range(n_r_c):
                coef += U[j,i]*b_c[j]
            coef *= i_S[i]
            for k in range(n_r_c):
                t_c[k] += coef*V[i,k]

        # Copy t_c values into t
        j = 0
        for i in range(n_r):
            if adjust[i]:
                t[i] = t_c[j]
                j += 1

        # Calculate chi-squared
        chi_2 = 0.
        for k_i in range(n_k):
            y_k = x_k[k_i]
            for n in range(n_r):
                y_k -= c_n_k[n,k_i]*t[n]
            chi_2 += self.opt_mag_ks[k_i][1]*y_k*y_k
        print '...chi^2 = ', chi_2,'...'

        # Compose real space part
        print '...composing real space part...'
        v_l_0 = 0.
        for n in range(n_r):
            v_l_0 += t[n]*basis.h(n,0.)
        v_l = np.zeros((self.n_points))
        for i in range(self.n_points):
            r = self.rs[i]
            if (r <= self.r_cut):
                for n in range(n_r):
                    v_l[i] += t[n]*basis.h(n,r)
            else:
                v_l[i] = v_r_spline(r)
        v_l_spline = spline(self.rs, v_l)

        # Get k=0 components (short)
        print '...computing k=0 components...'
        def v_short_integrand(r):
            if r < self.r_min:
                r = self.r_min
            return r*r*(v_r_spline(r) - v_l_spline(r))
        f_v_s_0 = -integrate(v_short_integrand, 1.e-100, self.r_cut, divmax=100)
        if (self.n_d == 2):
            f_v_s_0 *= 2.*pi/self.vol # FIXME: Probably wrong for 2D
        elif (self.n_d == 3):
            f_v_s_0 *= 4.*pi/self.vol

        # Get k=0 components (long)
        def v_long_integrand(r):
            if r < self.r_min:
                r = self.r_min
            return r*r*v_l_spline(r)
        f_v_l_0 = -integrate(v_long_integrand, 1.e-100, self.r_cut, divmax=100)
        if (self.n_d == 2):
            f_v_s_0 *= 2.*pi/self.vol # FIXME: Probably wrong for 2D
        elif (self.n_d == 3):
            f_v_l_0 *= 4.*pi/self.vol

        # Write r space part to file
        f = open(self.prefix+'_sq_'+self.object_string+'_diag_r.dat','w')
        f.write('%.10E %.10E\n'%(0.,v_l_0))
        for i in range(self.n_points):
            r = self.rs[i]
            f.write('%.10E %.10E\n'%(r,v_l[i]))
        f.close()

        # Compose k space part
        f = open(self.prefix+'_sq_'+self.object_string+'_diag_k.dat','w')
        f.write('%.10E %.10E\n'%(0.,f_v_s_0))
        f_v_ls = []
        for k in self.ks:
            k_2 = np.dot(k,k)
            mag_k = np.sqrt(k_2)
            f_v_l = 0.
            for n in range(n_r):
                f_v_l += t[n]*basis.c(n,mag_k)
            f_v_l -= self.CalcXk(v_r_spline, self.r_cut, mag_k, self.r_max)/self.vol
            f_v_ls.append([mag_k, f_v_l])
        mag_k_prev = -1
        for [mag_k,f_v_l] in sorted(f_v_ls):
            if abs(mag_k-mag_k_prev) > 1.e-8:
                f.write('%.10E %.10E\n'%(mag_k,f_v_l))
            mag_k_prev = mag_k
        f.close()

    def ComputeMadelung(self):
        print 'Computing madelung constant from breakup...'

        # K space part
        data = np.loadtxt(self.prefix+'_sq_'+self.object_string+'_diag_k.dat')
        mag_ks = data[:,0]
        f_v_ls = data[:,1]/self.z_1_z_2
        f_v_l_0 = f_v_ls[0]

        # R space part
        data = np.loadtxt(self.prefix+'_sq_'+self.object_string+'_diag_r.dat')
        rs = data[1:,0]
        v_ls = data[1:,1]/self.z_1_z_2
        v_l_0 = v_ls[0]

        # Spline v short
        v_s_spline = spline(rs, (self.cofactor/rs) - v_ls)

        # Set up test system
        half_L = self.L/2.
        N = 8
        xs = []
        xs.append(np.array([0,0,0]))
        xs.append(np.array([half_L,half_L,0]))
        xs.append(np.array([half_L,0,half_L]))
        xs.append(np.array([0,half_L,half_L]))
        xs.append(np.array([half_L,0,0]))
        xs.append(np.array([0,half_L,0]))
        xs.append(np.array([0,0,half_L]))
        xs.append(np.array([half_L,half_L,half_L]))
        qs = np.array([1.,1.,1.,1.,-1.,-1.,-1.,-1.])

        # Compute short ranged part
        v_s = 0.
        for i in range(N-1):
            for j in range(i+1,N):
                r = xs[i] - xs[j]
                for d_i in range(self.n_d):
                    r[d_i] -= floor(r[d_i]*self.i_box[d_i] + 0.5)*self.L
                mag_r = np.sqrt(np.dot(r,r))
                if (mag_r <= self.r_cut):
                    v_s -= qs[i]*qs[i]*v_s_spline(mag_r)

        # Match k space pieces
        f_v_l = np.zeros((len(self.ks)))
        for i in range(len(self.ks)):
            mag_k_i = np.sqrt(np.dot(self.ks[i],self.ks[i]))
            found_me = 0
            for j in range(len(mag_ks)):
                if (abs(mag_k_i-mag_ks[j])<1.e-4):
                    if not found_me:
                        found_me = 1
                        f_v_l[i] = f_v_ls[j]

        # Compute long ranged part
        v_l = 0.
        for i in range(len(self.ks)):
            mag_k_i = np.sqrt(np.dot(self.ks[i],self.ks[i]))
            if (mag_k_i < self.k_cut):
                re,im = 0.,0.
                for j in range(N):
                    h = np.dot(self.ks[i],xs[j])
                    re += qs[j]*cos(h)
                    im -= qs[j]*sin(h)
                for j in range(N):
                    h = np.dot(self.ks[i],xs[j])
                    v_l += 0.5*qs[j]*(re*cos(h) - im*sin(h))*f_v_l[i]

        # Compute self interacting terms
        v_self = 0.
        for i in range(N):
            v_self -= 0.5*qs[i]*qs[i]*v_l_0

        # Compute neutralizing background
        v_b = 0.
        for i in range(N):
            v_b += 0.5*N*N*f_v_l_0

        # Compute Madelung constant
        print self.L*(v_s + v_l + v_self)/N

    def ComputeMadelungNaive(self,n_images):
        print 'Computing Madelung constant from images...'

        # Read potential and spline it
        data = np.loadtxt(self.prefix+'_sq_'+self.object_string+'_diag.dat')
        v_r_spline = spline(data[:,0], data[:,1]/self.z_1_z_2)
        r_max = data[-1,0]

        # Set up test system
        half_L = self.L/2.
        N = 8
        xs = []
        xs.append(np.array([0,0,0]))
        xs.append(np.array([half_L,half_L,0]))
        xs.append(np.array([half_L,0,half_L]))
        xs.append(np.array([0,half_L,half_L]))
        xs.append(np.array([half_L,0,0]))
        xs.append(np.array([0,half_L,0]))
        xs.append(np.array([0,0,half_L]))
        xs.append(np.array([half_L,half_L,half_L]))
        qs = np.array([1.,1.,1.,1.,-1.,-1.,-1.,-1.])

        # Compose images
        n_is = []
        for n_x in range(-n_images,n_images+1):
            if (self.n_d > 1):
                for n_y in range(-n_images,n_images+1):
                    if (self.n_d > 2):
                        for n_z in range(-n_images,n_images+1):
                            n_is.append(np.array([n_x,n_y,n_z]))
                    else:
                        n_is.append(np.array([n_x,n_y]))
            else:
                n_is.append(np.array([n_x]))

        # Compute sum over images FIXME: assumes Coulomb tail
        v_s = 0.
        for i in range(N-1):
            for j in range(i+1,N):
                r_0 = xs[i] - xs[j]
                for n_i in n_is:
                    r = r_0 + n_i*self.box
                    mag_r = np.sqrt(np.dot(r,r))
                    if (mag_r > r_max):
                        v_s += self.cofactor*qs[i]*qs[j]/mag_r
                    else:
                        v_s += qs[i]*qs[j]*v_r_spline(mag_r)

        # Compute self energy
        v_self = 0.
        for i in range(N):
            for n_i in n_is:
                r = n_i*self.box
                mag_r = np.sqrt(np.dot(r,r))
                if (mag_r == 0.):
                    v_self += 0.
                elif (mag_r > r_max):
                    v_self += self.cofactor*qs[i]*qs[i]/mag_r
                else:
                    v_self += qs[i]*qs[i]*v_r_spline(mag_r)

        # Compute madelung constant
        print self.box[0]*(v_s + 0.5*v_self)/N

    def ComputeMadelungExact(self):
        print 'Computing Madelung constant exactly (bare 3D Coulomb only)...'
        print '-1.747564594633182190636212035544397403481'

def run(settings,object_type,prefix,z_1_z_2,cofactor):
    # Set up system
    e = EwaldBreakup(settings['n_d'],settings['L'],settings['type'],object_type,prefix,settings['grid_type'],settings['r_min'],settings['r_max'],settings['n_grid'],settings['r_cut'],settings['k_cut'],z_1_z_2,cofactor,settings['n_knots'])

    # Do the breakup
    e.DoBreakup()
    e.ComputeMadelung()
    e.ComputeMadelungNaive(settings['n_images'])
    if (object_type == 'v'):
        e.ComputeMadelungExact()

def usage():
    print "Usage: %s n_d L breakup_type object_string prefix grid_type r_min r_max n_points r_cut k_cut z_1_z_2 cofactor n_knots" % os.path.basename(sys.argv[0])
    sys.exit(2)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    if "-h" in argv or "--help" in argv:
        usage()

    try:
        n_d = int(argv[1])
        L = float(argv[2])
        breakup_type = argv[3]
        object_string = argv[4]
        prefix = argv[5]
        grid_type = argv[6]
        r_min = float(argv[7])
        r_max = float(argv[8])
        n_points = int(argv[9])
        r_cut = float(argv[10])
        k_cut = float(argv[11])
        z_1_z_2 = float(argv[12])
        cofactor = float(argv[13])
        n_knots = int(argv[14])
        n_images = int(argv[15])
    except:
        usage()

    # Set up system
    e = EwaldBreakup(n_d,L,breakup_type,object_string,prefix,grid_type,r_min,r_max,n_points,r_cut,k_cut,z_1_z_2,cofactor,n_knots)

    # Do the breakup
    e.DoBreakup()
    e.ComputeMadelung()
    e.ComputeMadelungNaive(n_images)
    e.ComputeMadelungExact()

if __name__ == "__main__":
    sys.exit(main())

