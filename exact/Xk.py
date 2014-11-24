from math import sqrt, pi, cos, sin
import numpy as np
from scipy import integrate
import mpmath as mp
mp.dps = 15; mp.pretty = True
from scipy.interpolate import InterpolatedUnivariateSpline

def Addk(k, degeneracy, kpoints):
  try:
      kpoints[k] += degeneracy
  except:
      kpoints[k] = degeneracy

def SetkVecs(kc, kCont, kMax, kpoints, L):
  numk = 0
  b = [1./L,1./L,1./L]
  maxIndex = [0,0,0]
  for i in range(len(maxIndex)):
    b[i] = 2*pi*b[i]
    maxIndex[i] = int(np.ceil(kCont/b[i]))
  print maxIndex

  k = np.zeros((3))
  for ix in range(-maxIndex[0],maxIndex[0]+1):
    k[0] = ix*b[0]
    for iy in range(-maxIndex[1],maxIndex[1]+1):
      k[1] = iy*b[1]
      for iz in range(-maxIndex[2],maxIndex[2]+1):
        k[2] = iz*b[2]
        k2 = np.dot(k,k)
        if (k2 > kc*kc and k2 < kCont*kCont):
          Addk(sqrt(k2),1,kpoints)
          numk += 1

  # Now, add kpoints to the list with approximate degeneracy.
  kvol = b[0]*b[1]*b[2]
  N = 4000
  deltak = (kMax-kCont)/N
  for i in range(N):
    k1 = kCont + deltak*i
    k2 = k1 + deltak
    k = 0.5*(k1+k2)
    vol = 4.*pi/3.*(k2*k2*k2-k1*k1*k1)
    degeneracy = vol/kvol
    Addk(k, degeneracy, kpoints)
    numk += degeneracy

  print "Total k vecs = ", numk
  print "non-degenerate k vecs = ", len(kpoints)

def XkCoul(k, r, Z1Z2):
  return -(4.*pi*Z1Z2/(k*k))*cos(k*r)

# This calculates the quantity 
# \f$ X_k \equiv -\frac{4 \pi}{\Omega k} \int_{r_c}^\infty dr \, r \sin(kr) V(r).\f$
def CalcXk(V, k, rMax, rEnd, Z1Z2):
  def integrand(r):
      if r > rEnd:
          return sin(k*r)
      else:
          return r*sin(k*r)*V(r)
  # DIVIDE INTEGRATION INTO SEGMENTS AND ADD THEM TOGETHER

  ## First segment
  rFirst = rMax + ((pi/k)-(rMax % (pi/k)))
  Xk = -(4.0*pi/k) * integrate.quad(integrand, rMax, rFirst)[0]

  ## Subsequent segments
  rMaxMax = 20*rMax
  nPi = 256*pi
  nSeg = int((rMaxMax-rFirst)/(nPi/k))
  for i in range(nSeg):
      Xk += -(4.0*pi/k) * integrate.quad(integrand, rFirst + i*(nPi/k), rFirst + (i+1)*(nPi/k))[0]

  # Add in the analytic part that I ignored
  # Multiply analytic term by tau only for U -- do not multiply
  # for dU or V.
  Xk += XkCoul(k, rFirst + (nSeg)*(nPi/k), Z1Z2)
  return Xk


kpoints = {}
L = 10.
rMax = L/2.
Z1Z2 = 1.0
kvol = (2.*pi/L)**3
kavg = kvol**(1./3.)
kCont = 50.0 * kavg
delta = 0.384615
kMax = 20.0*pi/delta
kc = 1.88496
print kavg, kc, kCont, delta, kMax

# Get potential
print 'Retrieving potential'
rs,vs = [],[]
f = open('v.1.txt')
for line in f:
  line = line.split()
  rs.append(float(line[0]))
  vs.append(float(line[1]))
rEnd = rs[-1]
V = InterpolatedUnivariateSpline(np.array(rs), np.array(vs), k='3')

# Setup k vecs
print 'Setup k vecs'
SetkVecs(kc, kCont, kMax, kpoints, L)

# Calc Xk
print 'Calculate Xk'
tot = 0.
Xk = {}
for k in kpoints:
  tmpXk = XkCoul(k, rMax, Z1Z2)
  Xk[k] = CalcXk(V, k, rMax, rEnd, Z1Z2)
  tot += abs(tmpXk - Xk[k])
print tot
