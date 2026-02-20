import numpy as np
import struct
import sys

selectedplanes = ()     # list of planes to be extracted. Leave empty to extract all of them
sizeofdata = 4          # in bytes (4 for single precision, 8 for double)
# etype = ">"           # big-endian
etype = "<"             # little-endian
dtype = "f"             # floating-point number, single precision

# fileformat = 'raw'
fileformat = "NetCDF"

# do not edit below this line
import atlab

# getting data from stdin
if ( len(sys.argv) <= 2 ):
    print("Usage: python $0 [xy,xz,yz] list-of-files.")
    quit()

planetype  = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

nx, ny, nz = atlab.getGridSize()                        # getting grid size from tlab.ini, if necessary
x, y, z = atlab.getGrid(nx, ny, nz, 'grid')             # get grid as n-tuples

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
        fin = open(file, 'rb')

        raw = fin.read(4)
        sizeofheader = struct.unpack(etype+'{}i'.format(1), raw)[0]
        # print(sizeofheader)
        
        raw = fin.read(4)
        iteration = struct.unpack(etype+'{}i'.format(1), raw)[0]

        raw = fin.read(8)
        time = struct.unpack(etype+'{}d'.format(1), raw)[0]
        # print('Processing iteration {}, time {:e}.'.format(iteration,time))

        raw = fin.read(int(sizeofheader-4-4-8))
        setofplanes = struct.unpack(etype+'{}i'.format(int((sizeofheader-4-4-8)/4)), raw)
        # print('Planes {}.'.format(setofplanes))
        
        fin.seek(0, 2)                      # move the cursor to the end of the file
        sizeoffile = fin.tell()
        sizeofvars = int((sizeoffile-sizeofheader)/(sizeofplane*len(setofplanes)*sizeofdata))
        # print('Number of variables is {}'.format(sizeofvars))

        if not selectedplanes:              # plot all planes
            selectedplanes = setofplanes

        xa = np.array(x)
        ya = np.array(y)
        za = np.array(z)

        for var in range(sizeofvars):
            for index, plane in enumerate(setofplanes):
                if plane in selectedplanes:
                    fin.seek(sizeofheader +(index+var*len(setofplanes))*sizeofplane*sizeofdata,0)
                    raw = fin.read(sizeofplane*sizeofdata)

                    dst_file = file  +'.v'+str(var) +'.p'+tag(sizeofmask,plane)
                    match fileformat:
                        case "NetCDF":
                            match planetype:
                                case "xy":
                                    za = np.array([z[plane-1]])    # Note that the plane starts at 1 and not 0
                                    nz = len(za)
                                case "xz":
                                    ya = np.array([y[plane-1]])
                                    ny = len(ya)
                                case "yz":
                                    xa = np.array([x[plane-1]])
                                    nx = len(xa)
                            a = np.array(struct.unpack((etype+'{}'+dtype).format(int(nx*ny*nz)), raw))
                            atlab.writeNetCfd(dst_file, 'Var'+str(var), time, xa, ya, za, a)

                        case "raw":
                            fout = open(dst_file,'wb')
                            fout.write(raw)
                            fout.close()

        fin.close()
