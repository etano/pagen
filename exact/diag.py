import rhocoul
import sys
import os
from math import sqrt, exp, cos, sin, pi, log, log10, acos
import numpy as np

def run(tau,xkappa,z,D):
    import Cusp
    #du00 = Cusp.du00dBeta(tau,1./xkappa,z,D,1e-9)
    #print 'du00', du00

    u00 = Cusp.u00(tau,1./xkappa,z,D,1e-9)
    print 0.0, u00
    
    kmax = 10
    lmax = 80
    nmax = 100
    rMin = 0.1
    L = 6.0
    rMax = L/2.
    rN = 30
    gridType = "LINEAR"
    if gridType == "LOG":
        rGrid = np.logspace(log10(rMin), log10(rMax), num=rN, endpoint=True)
    elif gridType == "LINEAR":
        rGrid = np.linspace(rMin, rMax, num=rN, endpoint=True)
    elif gridType == "OPT":
        rGrid = [rMin]
        dr = 2.*(rMax-rMin)/rN
        a = 10.
        f = [-a]
        for i in range(1,rN):
            f.append(f[i-1] + 2.*a/rN)
            rGrid.append(rGrid[i-1] + (1. - (1./(exp(f[i])+1.)))*dr)
        rGrid.append(rMax)
    iFourXkappaTau = 1./(4.*xkappa*tau)
    iFourPiXkappaTauToHalfD = 1./(4.*pi*xkappa*tau)**(D/2.)
    urrs = []
    f = open('Vr.dat','w')
    g = open('Vr0.dat','w')
    for r in rGrid:
        ra = np.array((r,0.,0.))
        rb = np.array((r,0.,0.))
        r1 = sqrt(np.dot(ra,ra))
        r2 = sqrt(np.dot(rb,rb))
        if (r1 != 0 and r2 != 0):
            theta = acos(np.dot(ra,rb)/(r1*r2))
        else:
            theta = 0
        s2 = r1*r1 + r2*r2 - 2*r1*r2*cos(theta)
        coul = rhocoul.RhoCoul(r1,r2,theta,tau,xkappa,z,kmax,lmax,nmax,D)
        free = iFourPiXkappaTauToHalfD * exp(-s2*iFourXkappaTau)
        #print coul, free
        urr = -log(coul/free)
        print r, urr
        f.write('%f %f\n'%(r,urr/(xkappa*z*tau)))
        g.write('%f %f\n'%(r,u00/(xkappa*z*tau)))
    f.close()
    g.close()

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
    except:
        usage()

    run(tau,xkappa,z,D)

if __name__ == "__main__":
    sys.exit(main())
