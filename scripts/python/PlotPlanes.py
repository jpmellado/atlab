#!/usr/bin/python3

import numpy as np
import struct
import sys
import matplotlib.pyplot as plt

sizeofdata = 4 # in bytes
# sizeofdata = 1 # for gate files

# etype = ">" # big-endian
etype = "<" # little-endian

dtype = "f" # floating number
# dtype = 'B' # unsigned character, for gate files

# do not edit below this line
import atlab
nx, ny, nz = atlab.getGridSize()                        # getting grid size from tlab.ini, if necessary

# getting data from stdin
if ( len(sys.argv) <= 2 ):
    print("Usage: python $0 [xy,xz,yz] list-of-files.")
    quit()

planetype  = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

# processing data
fin = open('grid.'+planetype, 'rb')
if   ( planetype == 'xy' ):
    raw = fin.read(nx*4)
    x1 = np.array(struct.unpack(etype+'{}f'.format(nx), raw))
    raw = fin.read(ny*4)
    x2 = np.array(struct.unpack(etype+'{}f'.format(ny), raw))
elif ( planetype == 'xz' ):
    raw = fin.read(nx*4)
    x1 = np.array(struct.unpack(etype+'{}f'.format(nx), raw))
    raw = fin.read(nz*4)
    x2 = np.array(struct.unpack(etype+'{}f'.format(nz), raw))
elif ( planetype == 'yz' ):
    raw = fin.read(ny*4)
    x1 = np.array(struct.unpack(etype+'{}f'.format(ny), raw))
    raw = fin.read(nz*4)
    x2 = np.array(struct.unpack(etype+'{}f'.format(nz), raw))
fin.close()

nx1 = len(x1)
nx2 = len(x2)

for file in setoffiles:
    print("Processing file %s ..." % file)
    fin = open(file, 'rb')
    raw = fin.read()
    a = np.array(struct.unpack((etype+'{}'+dtype).format(int(fin.tell()/sizeofdata)), raw))
    a = a.reshape((nx2,nx1))
    fin.close()

    mean, std = np.mean( a ), np.std( a )
    plt.figure( figsize=(10,8) )
    plt.pcolormesh(x1,x2,a,shading='auto',vmin=mean-std,vmax=mean+std)
    # plt.contourf(x1,x2,a)
    plt.axis('equal')
    # plt.gca().set_xlim([x1[0],x1[-1]])
    # plt.gca().set_ylim([x2[0],x2[-1]])
    # plt.axis([ x1[0], x1[-1], x2[0], x2[-1]])
    plt.colorbar()
    plt.tight_layout(pad=2)
    plt.title(file)
    # plt.savefig("{}.jpg".format(file),dpi=150,bbox_inches='tight')
    plt.show()
    plt.close()
