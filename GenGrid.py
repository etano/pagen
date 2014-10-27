from numpy import linspace, logspace, log10, exp

def GenGrid(o):
    rs = []
    if o['gridType']=="LINEAR":
        rs = linspace(o['rMin'],o['rMax'],num=o['nGrid'],endpoint=True)
    elif o['gridType']=="LOG":
        rs = logspace(log10(o['rMin']),log10(o['rMax']),num=o['nGrid'],endpoint=True)
    elif o['gridType']=="OPTIMIZED":
        rs = [o['rMin']]
        f0 = 0.
        a = exp(-0.58015)*pow(1.*o['nGrid'],0.506494) # Empirical factor
        dr = o['rMax']/(o['nGrid']-a)
        for iGrid in range(o['nGrid']):
            fi = f0 + 2.*(iGrid+1)*a/o['nGrid']
            rs.append(rs[iGrid] + (1. - (1./(exp(fi)+1.)))*dr)
    else:
        print 'Unrecognized grid'
    return rs
