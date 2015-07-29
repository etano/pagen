from numpy import linspace, logspace, log10, exp

def GenGrid(o):
    rs = []
    if o['grid_type']=="LINEAR":
        rs = linspace(o['r_min'],o['r_max'],num=o['n_grid'],endpoint=True)
    elif o['grid_type']=="LOG":
        rs = logspace(log10(o['r_min']),log10(o['r_max']),num=o['n_grid'],endpoint=True)
    elif o['grid_type']=="LOGLIN":
        dr=(o['r_max']-o['r_paste'])/(((o['n_grid']/2.)+1.)-1.)
        gdr=(o['r_paste']/o['r_min'])**(1./((o['n_grid']/2.)-1.))
        for grid_i in range(1,o['n_grid']+1):
            if(grid_i>o['n_grid']/2.):
                rs.append(o['r_paste']+dr*(grid_i-(o['n_grid']/2.)))
            else:
                rs.append(o['r_min']*gdr**(grid_i-1.))
    elif o['grid_type']=="OPTIMIZED":
        rs = [o['r_min']]
        f0 = 0.
        a = exp(-0.58015)*pow(1.*(o['n_grid']-1),0.506494) # Empirical factor
        dr = o['r_max']/((o['n_grid']-1)-a)
        for grid_i in range(o['n_grid']-1):
            fi = f0 + 2.*(grid_i+1)*a/(o['n_grid']-1)
            rs.append(rs[grid_i] + (1. - (1./(exp(fi)+1.)))*dr)
    else:
        print 'Unrecognized grid'
    return rs
