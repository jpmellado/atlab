#!/usr/bin/python3

import numpy as np
import struct
import sys

setofplanes = [ 1,2,3,4,5 ]

sizeofdata = 4 # in bytes
# sizeofdata = 1 # for gate files

sizeofheader = -1   # read it from header in each file
# sizeofheader = 0
# sizeofheader = 36 # for gate files
# sizeofheader = 52 # for flow files
# sizeofheader = 44 # for scal files

# etype = ">" # big-endian
etype = "<" # little-endian

nx = 0 # number of points in Ox; if 0, then search tlab.ini
ny = 0 # number of points in Oy; if 0, then search tlab.ini
nz = 0 # number of points in Oz; if 0, then search tlab.ini

# do not edit below this line

# getting grid size from tlab.ini, if necessary
if ( nx == 0 ):
    for line in open('tlab.ini'):
        if line.lower().replace(" ","").startswith("imax="):
            nx = int(line.split("=",1)[1])
            break

if ( ny == 0 ):
    for line in open('tlab.ini'):
        if line.lower().replace(" ","").startswith("jmax="):
            ny = int(line.split("=",1)[1])
            break

if ( nz == 0 ):
    for line in open('tlab.ini'):
        if line.lower().replace(" ","").startswith("kmax="):
            nz = int(line.split("=",1)[1])
            break

print("Grid size is {}x{}x{}.".format(nx,ny,nz))

# getting data from stdin
if ( len(sys.argv) <= 2 ):
    print("Usage: python $0 [xy,xz,yz] list-of-files.")
    quit()

planetype  = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

if   ( planetype == 'xy' ):
    sizeofmask = len(str(nz))
elif ( planetype == 'xz' ):
    sizeofmask = len(str(ny))
elif ( planetype == 'yz' ):
    sizeofmask = len(str(nx))
else:
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
    print("Processing file %s ..." % file)
    if file == 'grid':
        fin = open('grid', 'rb')
        #
        fin.seek(56,0)
        raw = fin.read(nx*8)
        a = struct.unpack(etype+'{}d'.format(nx), raw)
        rawx = struct.pack(etype+'{}f'.format(nx),*a)
        #
        fin.seek(8,1)
        raw = fin.read(ny*8)
        a = struct.unpack(etype+'{}d'.format(ny), raw)
        rawy = struct.pack(etype+'{}f'.format(ny),*a)
        #
        fin.seek(8,1)
        raw = fin.read(nz*8)
        a = struct.unpack(etype+'{}d'.format(nz), raw)
        rawz = struct.pack(etype+'{}f'.format(nz),*a)
        #
        fin.close()

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
        if sizeofheader == -1:
            raw = fin.read(4)
            sizeofheader = struct.unpack(etype+'{}i'.format(1), raw)[0]
            # print(sizeofheader)
            # raw = fin.read(int(sizeofheader-4-8))
            # planeIds = np.array(struct.unpack(etype+'{}i'.format(int((sizeofheader-4-8)/4)), raw))
            # print(planeIds)
            # raw = fin.read(8)
            # time = struct.unpack(etype+'{}d'.format(1), raw)
            # print(time)
        for plane in setofplanes:
            fout = open(file+'.'+planetype+tag(sizeofmask,plane),'wb')
            if   ( planetype == 'xy' ):
                fin.seek(sizeofheader +(plane-1)*nx*ny*sizeofdata,0)
                raw = fin.read(nx*ny*sizeofdata)
                fout.write(raw)
            elif ( planetype == 'xz' ):
                fin.seek(sizeofheader +(plane-1)*nx*sizeofdata,0)
                for k in range(nz):
                    raw = fin.read(nx*sizeofdata)
                    fout.write(raw)
                    fin.seek((ny-1)*nx*sizeofdata,1)
            elif ( planetype == 'yz' ):
                fin.seek(sizeofheader +(plane-1)*sizeofdata,0)
                for jk in range(ny*nz):
                    raw = fin.read(sizeofdata)
                    fout.write(raw)
                    fin.seek((nx-1)*sizeofdata,1)

            fout.close()
        fin.close()
