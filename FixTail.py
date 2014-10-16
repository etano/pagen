import sys
from scipy.optimize import curve_fit
from numpy import log10, array, exp

def FindMinMax(xs,ys,nPointsToFit):
    xys = list(zip(xs,ys))
    x_last = xys[-1][0]
    y_last = xys[-1][1]
    nPointsToFit = 10
    nTrue = 0
    for i in range(1,len(xys)):
        x = xys[-i-1][0]
        y = xys[-i-1][1]
        if y > y_last:
            nTrue += 1
            if nTrue == 1:
                xMax = x
        else:
            nTrue = 0
        if nTrue == nPointsToFit:
            xMin = x
            break
        x_last = x
        y_last = y
    xMin = pow(10.,xMin)
    xMax = pow(10.,xMax)
    return xMin, xMax

def FixTail(fileName,nPointsToFit,asymptote):
    # Read file
    f = open(fileName,'r')
    xs,ys = [],[]
    for line in f:
        [x,y] = [float(i) for i in line.split()]
        xs.append(x)
        ys.append(y)
    f.close()
    
    # Find xmin and xmax
    logxs,logys = [],[]
    for (x,y) in zip(xs,ys):
        if (x > 0):
            logxs.append(log10(x))
            logys.append(log10(abs(y-asymptote(x))))
    xMin, xMax = FindMinMax(logxs,logys,nPointsToFit)
    print 'Setting xMin=', xMin, ', xMax=', xMax, '...'

    # Select data and take difference with bare Coulomb
    logxs,logys = [],[]
    for (x,y) in zip(xs,ys):
        if x > xMin and x < xMax:
            logxs.append(log10(x))
            logys.append(log10(abs(y-asymptote(x))))
    
    # Fit data
    def fn(x,m,b):
        return m*x + b
    popt,pcov = curve_fit(fn,logxs,logys)
    
    # Write fit to file
    f = open(fileName,'w')
    for (x,y) in zip(xs,ys):
        if x > xMin:
            f.write('%.10e %.10e\n' % (x, 10**fn(log10(x),*popt) + asymptote(x)))
        else:
            f.write('%.10e %.10e\n' % (x,y))
    f.close()
