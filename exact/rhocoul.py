import numpy as np
import coulcc
from math import sqrt, exp, cos, sin, pi, log, log10, acos
from numpy import unique
from scipy import integrate
import sys

def RhoCoul(r1,r2,theta,tau,xkappa,z,kmax,lmax,nmax,D):
    uklim = sqrt(float(kmax)/(tau*xkappa))
    rhocoul = 0.
    if D == 2:
        rhocoul = integrate.romberg(CGrand2D, 0., 2.*uklim, args=(r1,r2,theta,tau,xkappa,z,lmax),divmax=100)
    elif D == 3:
        rhocoul = integrate.romberg(CGrand3D, 0., 2.*uklim, args=(r1,r2,theta,tau,xkappa,z,lmax),divmax=100)
    if (z > 0.):
        return rhocoul
    else:
        if D == 2:
            rhocoul += Bound2D(r1,r2,theta,tau,xkappa,z,nmax)
        elif D == 3:
            rhocoul += Bound3D(r1,r2,theta,tau,xkappa,z,nmax)
    return rhocoul

def Bound2D(r1,r2,theta,tau,xkappa,z,nmax):
    epsbs = 1.e-9
    sumonn = 0.
    bound = 0.
    for n in range(1,nmax+1):
        twonm1in = 1./float(2*n-1)
        foo = 2.*abs(z)*twonm1in
        gibbs = exp(tau*xkappa*z*z*twonm1in*twonm1in)
        x1 = 2.*r1*abs(z)*twonm1in
        x2 = 2.*r2*abs(z)*twonm1in
        sumonl = 0.
        for l in range(0,n):
            nmlm1 = n-l-1
            a1 = 1.0
            a2 = 1.0
            for j in range(nmlm1,0,-1):
                reco = float(nmlm1-j+1)/float(j*(2*1+j))
                a1 = 1. - reco*x1*a1
                a2 = 1. - reco*x2*a2
            if (l == 0):
                fnorm = 1.
                combo = 1.
                power = foo
                r1tol = 1.
                r2tol = 1.
            else:
                fnorm /= float((n-l)*(n+l-1))
                combo *= float((n-l)*(n+l-1))/float(2*l*(2*l-1))
                power *= foo
                r1tol *= r1
                r2tol *= r2
            polyr1 = combo*a1
            polyr2 = combo*a2
            radfac = power*sqrt(fnorm*twonm1in)
            rad1 = radfac*exp(-.5*x1)*r1tol*polyr1
            rad2 = radfac*exp(-.5*x2)*r2tol*polyr2
            if (l == 0):
                sumonl = rad1*rad2/(2.*pi)
            else:
                sumonl += 2.*rad1*rad2*cos(l*theta)/(2.*pi)
        sumonn += gibbs*sumonl
        oldbound = bound
        tail = (n-1)*(n-1)*gibbs*sumonl/(2.*n-1.)
        bound = sumonn + tail
        if (abs(bound-oldbound) < epsbs*abs(bound)):
            return bound
    return bound

def Bound3D(r1,r2,theta,tau,xkappa,z,nmax):
    epsbd = 1.e-9
    sumonn = 0.
    bound = 0.
    r12 = sqrt(r1*r1 + r2*r2 - 2.*r1*r2*cos(theta))
    x = (r1 + r2 + r12)/2.
    y = (r1 + r2 - r12)/2.
    zabs = abs(z)
    for n in range(1,nmax+1):
        fxp = 1.
        zxn = zabs*x/float(n)
        fx = 1. - zxn
        fyp = 1.
        zyn = zabs*y/float(n)
        fy = 1. - zyn
        if (n != 1):
            for i in range(1,n):
                fxn = ((2.*i + 1. - zxn)*fx - i*fxp)/(i + 1.)
                fxp = fx
                fx = fxn
                if (x != y):
                    fyn = ((2.*i + 1. - zyn)*fy - i*fyp)/(i + 1.)
                    fyp = fy
                    fy = fyn
        prefac = (1./(4.*pi))*exp(.25*tau*xkappa*z*z/(n*n))*(.5*(zabs**3)/(n**5))*exp(-.5*zabs*(r1+r2)/n)
        if (x != y):
            term = ((n**3)/zabs)*prefac*(fxp*fy-fx*fyp)/r12
        elif (x != 0.):
            term = (n**2)*prefac*((n**2/(zabs*r1))*((fx-fxp)**2) + fx*fxp)
        else:
            term = (n**2)*prefac
        sumonn += term
        oldbound = bound
        tail = (n-1)*(n-1)*term/float(2*n-1)
        bound = sumonn + tail
        if (abs(bound-oldbound) < epsbd*abs(bound)):
            return bound
    return bound

def CGrand2D(wavek,r1,r2,theta,tau,xkappa,z,lmax):
    cgrand = 0.
    if (wavek == 0.):
        return cgrand
    ifail = 0
    r1k = complex(wavek*r1, 0.)
    r2k = complex(wavek*r2, 0.)
    fo1 = complex(z/(2.*wavek), 0.)
    gibbs = exp(-tau*xkappa*wavek*wavek)

    fc = np.asfortranarray(np.zeros(lmax), dtype='F')
    fcp = np.asfortranarray(np.zeros(lmax), dtype='F')
    gc = np.asfortranarray(np.zeros(lmax), dtype='F')
    gcp = np.asfortranarray(np.zeros(lmax), dtype='F')
    sig = np.asfortranarray(np.zeros(lmax), dtype='F')
    r1kl = np.asfortranarray(np.zeros(lmax), dtype='f')
    r2kl = np.asfortranarray(np.zeros(lmax), dtype='f')

    if (r1 != 0.):
        c1 = sqrt(2./(pi*r1))
        fc,gc,fcp,gcp,sig = coulcc.coulcc(r1k,fo1,-.5,lmax,1,0,ifail)
        for l in range(0,lmax):
            r1kl[l] = c1*fc[l].real
    else:
        r1kl[0] = sqrt(2.*wavek/(exp(pi*z/wavek)+1.))
        for l in range(1,lmax):
            r1kl[l] = 0.

    if (r1 == r2):
        for l in range(0,lmax):
            r2kl[l] = r1kl[l]
    elif (r2 != 0.):
        c2 = sqrt(2./(pi*r2))
        fc,gc,fcp,gcp,sig = coulcc.coulcc(r2k,fo1,-.5,lmax,1,0,ifail)
        for l in range(0,lmax):
            r2kl[l] = c2*fc[l].real
    else:
        r2kl[0] = sqrt(2.*wavek/(exp(pi*z/wavek)+1.))
        for l in range(1,lmax):
            r2kl[l] = 0.

    sum = r1kl[0]*r2kl[0]
    for l in range(1,lmax):
        sum += 2.*r1kl[l]*r2kl[l]*cos(l*theta)
    sum *= gibbs/(2.*pi)

    return sum

def CGrand3D(wavek,r1,r2,theta,tau,xkappa,z,lmax):
    cgrand = 0.
    if (wavek == 0.):
        return cgrand
    ifail = 0
    r12 = sqrt(r1*r1 + r2*r2 - 2.*r1*r2*cos(theta))
    x = (r1+r2+r12)/2.
    y = (r1+r2-r12)/2.
    xk = wavek*x
    yk = wavek*y
    xkz = complex(xk, 0.)
    ykz = complex(yk, 0.)
    fo1 = z/(2.*wavek)
    gibbs = exp(-tau*xkappa*wavek*wavek)
    eta1 = complex(fo1, 0.)

    fc = np.asfortranarray(np.zeros(1), dtype='F')
    fcp = np.asfortranarray(np.zeros(1), dtype='F')
    gc = np.asfortranarray(np.zeros(1), dtype='F')
    gcp = np.asfortranarray(np.zeros(1), dtype='F')
    sig = np.asfortranarray(np.zeros(1), dtype='F')

    if (x == y):
        if (x == 0.):
            cgrand = (z/(2.*pi))*wavek*gibbs/(exp(z*pi/wavek)-1.)
        else:
            fc,gc,fcp,gcp,sig = coulcc.coulcc(xkz,eta1,0.,1,3,0,ifail)
            fx = fc[0].real
            fxd = fcp[0].real
            fo2 = wavek*wavek*gibbs*(fxd*fxd + (1. - z/(wavek*wavek*x))*fx*fx)/(2.*pi*pi)
            cgrand = fo2
    else:
        fc,gc,fcp,gcp,sig = coulcc.coulcc(xkz,eta1,0,1,3,0,ifail)
        fx = fc[0].real
        fxd = fcp[0].real
        if (y < 1.e-8):
            fy = 0.
            fyd = sqrt(2.*pi*fo1/(exp(2.*pi*fo1)-1.))
        else:
            fc,gc,fcp,gcp,sig = coulcc.coulcc(ykz,eta1,0,1,3,0,ifail)
            fy = fc[0].real
            fyd = fcp[0].real
        cgrand = wavek*gibbs*(fx*fyd - fy*fxd)/(2.*pi*pi*r12)
    return cgrand

def Addk(k, ks, degeneracy=1):
    ki = 0
    while (ki < len(ks) and abs(k-ks[ki][0]) > 1.0e-12):
      ki += 1
    if ki == len(ks):
        ks.append([k,degeneracy])
    else:
        ks[ki][1] += degeneracy

def Genks(kMax,D,L):
    b = 2.*pi/L
    nMax = int((kMax/b) + 1.)
    print nMax

    ks = []
    for n1 in range(-nMax,nMax):
        if D == 1:
            k = 2*pi*np.array([n1])/L
            Addk(sqrt(np.dot(k,k)), ks)
        else:
            for n2 in range(-nMax,nMax):
                if D == 2:
                    k = 2*pi*np.array([n1,n2])/L
                    Addk(sqrt(np.dot(k,k)), ks)
                else:
                    for n3 in range(-nMax,nMax):
                        k = 2.*pi*np.array([n1,n2,n3])/L
                        Addk(sqrt(np.dot(k,k)), ks)
    return ks

def AddContinuumks(ks, kMin, kMax):
    kVol = 1.
    for i in range(0,D):
        kVol *= 2*pi/L
    kAvg = kvol**(1./D)

    N = 4000
    deltak = (kMax-kMin)/N
    nK = len(ks)
    for i in range(0,N):
      k1 = kMin + deltak*i
      k2 = k1 + deltak
      k = 0.5*(k1+k2)
      vol = (4.*pi/3.)*(k2*k2*k2-k1*k1*k1)
      degeneracy = vol/kVol
      Addk(k, ks, degeneracy)
      nK += degeneracy

#tau = 0.1
#lam1 = 0.5
#lam2 = 0.000272
#lam2 = 0.5
#Z1Z2 = 1.0
#D = 3
#
#lam = lam1*lam2/(lam1 + lam2)
#m1 = 1./(2.*lam1)
#m2 = 1./(2.*lam2)
#m12 = m1*m2/(m1+m2)
#lam12 = 1./(2.*m12)
#xkappa = lam12
#z = Z1Z2/xkappa
#
#import Cusp
#u00 = Cusp.u00(tau,1./xkappa,z,D,1e-9)
#print 'u00', u00
#
#kmax = 5
#lmax = 80
#nmax = 100
#rMin = 1e-4
#L = 5.0
#rMax = L/2.
#rN = 100
#gridType = "OPT"
#if gridType == "LOG":
#    rGrid = np.logspace(log10(rMin), log10(rMax), num=rN, endpoint=True)
#elif gridType == "OPT":
#    rGrid = [rMin]
#    dr = 2.*(rMax-rMin)/rN
#    a = 10.
#    f = [-a]
#    for i in range(1,rN):
#        f.append(f[i-1] + 2.*a/rN)
#        rGrid.append(rGrid[i-1] + (1. - (1./(exp(f[i])+1.)))*dr)
#    rGrid.append(rMax)
#f = open('Vr.dat','w')
#thetaGrid = np.linspace(0., 2.*pi, num=10)
#thetaGrid = 0.
#iFourXkappaTau = 1./(4.*xkappa*tau)
#iFourPiXkappaTauToHalfD = 1./(4.*pi*xkappa*tau)**(D/2.)
#urrs = []
#f = open('Vr.dat','w')
#for r in rGrid:
#    print r
#    for rp in rGrid:
#        for theta in thetaGrid:
#            #urr = 0.
#            #urr_primitive = 0.
#            #ra = np.array((r,0.,0.))
#            #rb = np.array((r,0.,0.))
#            #r1 = sqrt(np.dot(ra,ra))
#            #r2 = sqrt(np.dot(rb,rb))
#            #if (r1 != 0 and r2 != 0):
#            #    print r1, r2, (-s*s + 2.*q*q + 0.5*z*z)/(2.*q*q - 0.5*z*z)
#            #    theta = acos(np.dot(ra,rb)/(r1*r2))
#            #    theta = acos((-s*s + 2.*q*q + 0.5*z*z)/(2.*q*q - 0.5*z*z))
#            #else:
#            #    theta = 0
#            #q = 0.5*(r + rp)
#            s2 = r*r + rp*rp - 2*r*rp*cos(theta)
#            #z = r - rp
#            rhocoul = RhoCoul(r,rp,theta,tau,xkappa,Z1Z2/xkappa,kmax,lmax,nmax,D)
#            rhofree = iFourPiXkappaTauToHalfD * exp(-s2*iFourXkappaTau)
#            urr = -log(rhocoul/rhofree)
#            #urr_primitive += (tau/2.)*Z1Z2*((1./r1) + (1./r2))
#            #durr = urr - urr_primitive
#            #urrs.append(urr)
#            #vc = tau*(Z1Z2/r) + durr
#            #print r, urr, urr_primitive, durr, vc
#            #print q, s, z, urr
#            f.write('%f %f %f %f\n'%(r,rp,theta,urr))
#f.close()
##from scipy.interpolate import interp1d
##from cmath import exp as cexp
##u = interp1d(rGrid,np.array(urrs), kind='cubic')
##print 'Generating ks'
###kAvg = 2.*pi/L
###ks = Genks(4*kAvg,D,L)
###kCont = 50.*kAvg
###delta = 2.5/(20-1)
###kMax = 20.*pi/delta
##print 'Adding continuum ks'
###AddContinuumks(ks, kCont, kMax)
##f = open('Vk.dat','w')
##g = open('Vr.0.dat','r')
##for line in g:
##    k = float(line.split()[0])
##    q = k
##    if q != 0.:
##        uk_integrand = lambda r: 4.*pi*r*sin(q*r)*u(r)/q
##        uk_short = integrate.quad(uk_integrand,rMin,rMax)[0]
##        uk_long = 4.*pi*Z1Z2*tau*cos(q*rMax)/(q*q)
##        uk_coul = 4.*pi*Z1Z2*tau/(q*q)
##        f.write('%f %f\n'%(q,q*q*(uk_short+uk_long)/(4.*pi*Z1Z2*tau)))
##        #f.write('%f %f\n'%(q,q*q*(uk_coul)/(4.*pi*Z1Z2*tau)))
##g.close()
##f.close()
