import rhocoul
from math import sqrt, exp, cos, sin, pi, log, log10, acos
import numpy as np

tau = 0.125
lam1 = 0.5
lam2 = 0.5
lam2 = 0.00002285749742
Z1Z2 = -1.0
D = 3

lam = lam1*lam2/(lam1 + lam2)
m1 = 1./(2.*lam1)
m2 = 1./(2.*lam2)
m12 = m1*m2/(m1+m2)
lam12 = 1./(2.*m12)
xkappa = lam12
z = Z1Z2/xkappa

import Cusp
u00 = Cusp.u00(tau,1./xkappa,z,D,1e-9)
print 'u00', u00

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
    coul = rhocoul.RhoCoul(r1,r2,theta,tau,xkappa,Z1Z2/xkappa,kmax,lmax,nmax,D)
    free = iFourPiXkappaTauToHalfD * exp(-s2*iFourXkappaTau)
    #print coul, free
    urr = -log(coul/free)
    print r, urr
    f.write('%f %f\n'%(r,urr/(Z1Z2*tau)))
    g.write('%f %f\n'%(r,u00/(Z1Z2*tau)))
f.close()
g.close()
