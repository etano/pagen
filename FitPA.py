import sys
from math import pi, sqrt
from numpy import linspace, logspace, array, meshgrid, ravel
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# Inputs
nOrder = int(sys.argv[1])
fileName = sys.argv[2]

# Get data
f = open(fileName,'r')
xs,ys,zs = [],[],[]
for line in f:
    if len(line.split()) == 3:
        [x,y,z] = [float(a) for a in line.split()]
        xs.append(x)
        ys.append(y)
        zs.append(z)

# Make spline
spline = interpolate.SmoothBivariateSpline(array(xs),array(ys),array(zs),kx=3,ky=3)

# Plot original data and spline
X, Y = meshgrid(xs, ys)
zs = array([spline(x,y) for x,y in zip(ravel(X), ravel(Y))])
Z = zs.reshape(X.shape)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z)
plt.show()

