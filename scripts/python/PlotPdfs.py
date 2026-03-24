import netCDF4 as nc            # Reading NetCDF4 files.
import numpy   as np            # Array operations.
import matplotlib.pyplot as plt # Plot data
import sys

from matplotlib import rc       # Globally setup figure options
rc('text',       usetex=True)
rc('text.latex', preamble=r"\usepackage{fourier}")
rc('font',       family='serif', size=12)
rc('grid',       linestyle='dotted')
rc('axes',       grid=True)

# At home, screen 27''
rc('figure',     dpi=200)
rc('savefig',    dpi=100)

# getting data from stdin
if ( len(sys.argv) <= 2 ):
    print("Usage: python $0 variable list-of-files.")
    quit()

var = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

for file in setoffiles:
    plt.figure(figsize=(5,4))
    axes = plt.gca()

    data = nc.Dataset(file, 'r')
    data.set_auto_mask(False)

    t=data.variables['t'][:]
    z=data.variables['z'][:]
    f=data.variables[var][:,:,:]            # the 1. index is time, the 2. is vertical node, the 3. bins
    r=data.variables[var+'_range'][:,:,:]   # the 1. index is time, the 2. is vertical node
    
    nb = np.size(f,2)                       # number of bins, for readability below
    
    # create reference mesh
    bins = np.linspace(0., 1., nb)
    bins2d, z2d = np.meshgrid(bins, z)

    # consider only the first time in each file
    its = [0]
    for it in its:  
        for k in range(np.size(z)):             # update to the local range
            bins2d[k,:] = np.linspace(r[it, k, 0], r[it, k, 1], nb)
        
        plt.contourf(bins2d, z2d, f[it,:,:])
        plt.colorbar(label='count',format="%.2g")
        print('Total sample size {}'.format(np.sum(f[it,:,:])))
        
    plt.xlabel(r'{}'.format(var))
    plt.ylabel(r'height $z$')
    plt.title(r'It = {} in {}'.format(it,file))
    axes.spines['right'].set_visible(False)
    axes.spines['left'].set_position(('axes',-0.01))
    axes.get_yaxis().tick_left()
    axes.spines['top'].set_visible(False)
    axes.spines['bottom'].set_position(('axes',-0.01))
    axes.get_xaxis().tick_bottom()
    plt.tight_layout(pad=0.1)

    plt.savefig('figure1.pdf',bbox_inches='tight')
    
plt.show()
