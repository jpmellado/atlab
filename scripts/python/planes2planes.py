import numpy as np
import struct
import sys

selectedplanes = () # list of planes to be extracted. Leave empty to extract all of them
        
sizeofdata = 4      # in bytes

sizeofheader = -1   # read it from header in each file

# etype = ">"         # big-endian
etype = "<"         # little-endian

import atlab

# getting data from stdin
if ( len(sys.argv) <= 2 ):
    print("Usage: python $0 [xy,xz,yz] list-of-files.")
    quit()

planetype  = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

nx, ny, nz = atlab.getGridSize()                        # getting grid size from tlab.ini, if necessary

match planetype:
    case "xy":
        sizeofplane = nx*ny
        sizeofmask = len(str(nz))
    case "xz":
        sizeofplane = nx*nz
        sizeofmask = len(str(ny))
    case "yz":
        sizeofplane = ny*nz
        sizeofmask = len(str(nx))
    case _:
        print("Usage: python $0 [xy,xz,yz] list-of-files.")
        quit()

# further initialization
def tag(sizeofmask,number):
    a = str(number)
    for i in range(sizeofmask-len(a)):
        a = '0' + a
    return a

# extracting planes
for file in setoffiles:
    print("Processing file {:s} ...".format(file))
    if file == 'grid':
        x, y, z = atlab.getGrid(nx, ny, nz, 'grid')     # get grid as n-tuples
        rawx = struct.pack(etype+'{}f'.format(nx),*x)
        rawy = struct.pack(etype+'{}f'.format(ny),*y)
        rawz = struct.pack(etype+'{}f'.format(nz),*z)

        fout = open('grid.'+planetype,'wb')
        if   ( planetype == 'xy' ):
            fout.write(rawx)
            fout.write(rawy)
        elif ( planetype == 'xz' ):
            fout.write(rawx)
            fout.write(rawz)
        elif ( planetype == 'yz' ):
            fout.write(rawy)
            fout.write(rawz)
        fout.close()

    else:
        if sizeofheader == -1:
            fin = open(file, 'rb')
            raw = fin.read(4)
            sizeofheader = struct.unpack(etype+'{}i'.format(1), raw)[0]
            print(sizeofheader)
        
        fin = open(file, 'rb')        
        fin.seek(4)

        raw = fin.read(4)
        iteration = struct.unpack(etype+'{}i'.format(1), raw)[0]

        raw = fin.read(8)
        time = struct.unpack(etype+'{}d'.format(1), raw)[0]
        print('Processing iteration {}, time {:e}.'.format(iteration,time))

        raw = fin.read(int(sizeofheader-4-4-8))
        setofplanes = struct.unpack(etype+'{}i'.format(int((sizeofheader-4-4-8)/4)), raw)
        print('Planes {}.'.format(setofplanes))
        
        if not selectedplanes:              # plot all planes
            selectedplanes = setofplanes

        sizeofvars = 5

        for var in range(sizeofvars):
            for index, plane in enumerate(setofplanes):
                if plane in selectedplanes:
                    fout = open(file  +'.'+str(var) +'.'+tag(sizeofmask,plane),'wb')
                    fin.seek(sizeofheader +(index+var*len(setofplanes))*sizeofplane*sizeofdata,0)
                    raw = fin.read(sizeofplane*sizeofdata)
                    fout.write(raw)

                    fout.close()

        fin.close()
