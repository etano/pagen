import sys
from numpy import linspace, logspace, array, meshgrid, ravel, unique, log10, zeros
from scipy import interpolate, optimize
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# Inputs
nOrder = int(sys.argv[1])
fileName = sys.argv[2]

# Get data
print 'Reading data...'
f = open(fileName,'r')
xs,ys,zs,xds,yds,zds = [],[],[],[],[],[]
for line in f:
    if len(line.split()) == 3:
        [x,y,z] = [float(a) for a in line.split()]
        xs.append(x)
        ys.append(y)
        zs.append(z)
        if x == y:
            xds.append(x)
            yds.append(z)
            zds.append(z)

# Spline z
print 'Interpolating...'
u_xs = unique(array(xs))
u_ys = unique(array(ys))
u_zs = zeros((len(u_xs),len(u_ys)))
for i in range(len(u_xs)):
    for j in range(len(u_ys)):
        u_zs[i][j] = zs[i*len(u_ys) + j]
z = interpolate.RectBivariateSpline(u_xs,u_ys,u_zs)

# Plot original data and spline
print 'Plotting...'
X, Y = meshgrid(u_xs, u_ys)
Z = array([z(x,y) for x,y in zip(ravel(X), ravel(Y))])
Z = Z.reshape(X.shape)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xs, ys, zs)
ax.plot_surface(X, Y, Z)
plt.show()

# Subtract out diagonal
print 'Subtracting diagonal...'
zd = interpolate.InterpolatedUnivariateSpline(xds,zds)
zmzds = []
for i in range(len(xs)):
    if xs[i] == ys[i]:
        zmzds.append(0.)
    else:
        q = (xs[i] + ys[i])/2.
        zmzds.append(zs[i] - zd(q))
        #zmzds.append(zs[i] - (xs[i]*zd(xs[i]) + ys[i]*zd(ys[i]))/(xs[i] + ys[i]))
        #zmzds.append(zs[i] - (zd(xs[i]) + zd(ys[i]))/2.)

# Spline z - zd
print 'Interpolating...'
u_xs = unique(array(xs))
u_ys = unique(array(ys))
u_zmzds = zeros((len(u_xs),len(u_ys)))
for i in range(len(u_xs)):
    for j in range(len(u_ys)):
        u_zmzds[i][j] = zmzds[i*len(u_ys) + j]
zmzd = interpolate.RectBivariateSpline(u_xs,u_ys,u_zmzds)

# Plot original data and spline
print 'Plotting...'
X, Y = meshgrid(u_xs, u_ys)
Z = array([zmzd(x,y) for x,y in zip(ravel(X), ravel(Y))])
Z = Z.reshape(X.shape)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xs, ys, zmzds)
ax.plot_surface(X, Y, Z)
plt.show()

# Fitting function to off-diagonal terms
def f_zod(s, *A):
    tot = 0.
    for i in range(nOrder):
        tot += A[i] * (s**(2.*i))
    return tot

# Get all unique q values
qs = []
for x,y in zip(xs,ys):
    qs.append((x+y)/2.)
qs = unique(qs)

# Fit A's for each unique q
print 'Fitting A\'s...'
nPts = 100
As = []
for q in qs:
    s_space = linspace(0.,2.*q,num=nPts)
    zmzd_space = zeros((nPts))
    for i in range(nPts):
        zmzd_space[i] = zmzd(q + s_space[i]/2., q - s_space[i]/2.)
    popt,pcov = optimize.curve_fit(f_zod,s_space,zmzd_space,p0=zeros((nOrder)),xtol=1.e-10,ftol=1.e-10)
    As.append(popt)
As = array(As).T

# Spline each A
Ais = []
for i in range(nOrder):
    Ais.append(interpolate.InterpolatedUnivariateSpline(qs,As[i]))

# Plot As
print 'Plotting...'
for i in range(nOrder):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(qs,As[i],'o')
    t_qs = logspace(log10(1e-8),log10(max(qs)),num=100)
    t_Ais = Ais[i](t_qs)
    ax.set_xscale('log')
    ax.plot(t_qs,t_Ais,'-')
    plt.show()

# Spline z - zd
print 'Interpolating...'
u_zs_new = zeros((len(u_xs),len(u_ys)))
for i in range(len(u_xs)):
    for j in range(len(u_ys)):
        q = (u_xs[i] + u_ys[j])/2.
        s = u_xs[i] - u_ys[j]
        params = []
        for i in range(nOrder):
            params.append(Ais[i](q))
        u_zs_new[i][j] = zd(q) + f_zod(s,*params)
z_new = interpolate.RectBivariateSpline(u_xs,u_ys,u_zs_new)

# Plot original data and spline
print 'Plotting...'
X, Y = meshgrid(u_xs, u_ys)
Z = array([z_new(x,y) for x,y in zip(ravel(X), ravel(Y))])
Z = Z.reshape(X.shape)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xs, ys, zs)
ax.plot_surface(X, Y, Z)
plt.show()
