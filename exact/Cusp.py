# Pollock, E. L., Comp. Phys. Comm. 52 , 49 (1988).
# http://www.sciencedirect.com/science/article/pii/0010465588901713

import sys
import os
from scipy.misc import factorial, comb
from scipy.special import zeta, gamma
from math import pi, sqrt

def kappa(n, D, kappas, mus):
    if kappas[n-1] != 0.:
        return kappas[n-1]
    else:
        tot = -mu(n, D, mus)
        for m in range(1,n):
            tot -= (float(m)/float(n))*kappa(m, D, kappas, mus)*mu(n-m, D, mus)
        kappas[n-1] = tot
        return tot

def mu(n, D, mus):
    if mus[n-1] != 0.:
        return mus[n-1]
    else:
        u = 0.
        if D == 3:
            if n == 1:
                u = -sqrt(pi)
            elif(n % 2 == 0):
                u = float((-1)**(n-2))/float(factorial(n-2,exact=1)) * factorial((n-2)/2,exact=1) * zeta(n,1)
            else:
                u = float((-1)**(n-2))/float(factorial(n-2,exact=1)) * gamma(1.+(n-2.)/2.) * zeta(n,1)
        elif D ==2:
            if(n % 2 == 0):
                u = float((-1)**(n-1))/float(factorial(n-1,exact=1)) * gamma(1.+(n-1.)/2.) * (float(1-(2**(n+1)))/sqrt(pi)) * zeta(n+1,1)
            else:
                u = float((-1)**(n-1))/float(factorial(n-1,exact=1)) * factorial((n-1)/2,exact=1) * (float(1-(2**(n+1)))/sqrt(pi)) * zeta(n+1,1)
        mus[n-1] = u
        return u

def GetPs(nOrder, D):
    kappas = [0.]*nOrder
    mus = [0.]*nOrder
    kappa(nOrder, D, kappas, mus)
    return kappas

def gam(tau, lam, Z1Z2):
    return tau*(Z1Z2*Z1Z2)/lam

def u00(tau, lam, Z1Z2, D, tol):
    nOrder = 10
    g = gam(tau, lam, Z1Z2)
    P = GetPs(nOrder, D)
    #for p in P:
    #    print p
    u = 0.
    oldu = 1.e100
    for j in range(1,len(P)+1):
        if (Z1Z2 > 0.):
            u += P[j-1] * (g**(j/2.))
        else:
            u += ((-1)**j) * P[j-1] * (g**(j/2.))
        if abs(u-oldu) < tol:
            return u
        else:
            oldu = u
    print 'Warning: answer not within given tolerance!'
    return u

def du00dBeta(tau, lam, Z1Z2, D, tol):
    dtau = 1e-8
    nOrder = 50
    return (u00(tau+dtau,lam,Z1Z2,D,tol) - u00(tau,lam,Z1Z2,D,tol))/dtau

def usage():
    print "Usage: %s tau lam1 lam2 Z1Z2 D" % os.path.basename(sys.argv[0])
    sys.exit(2)

def main(argv=None):
    if argv is None:
        argv = sys.argv
    if "-h" in argv or "--help" in argv:
        usage()

    try:
        tau = float(sys.argv[1])
        lam1 = float(sys.argv[2])
        lam2 = float(sys.argv[3])
        Z1Z2 = float(sys.argv[4])
        D = int(sys.argv[5])
        lam = lam1*lam2/(lam1 + lam2)
        m1 = 1./(2.*lam1)
        m2 = 1./(2.*lam2)
        m12 = m1*m2/(m1+m2)
        lam12 = 1./(2.*m12)
        xkappa = lam12
        z = Z1Z2/xkappa
        print u00(tau,1/xkappa,z,D,1e-6)
        print du00dBeta(tau,1./xkappa,z,D,1e-4)
    except:
        usage()

if __name__ == "__main__":
    sys.exit(main())
