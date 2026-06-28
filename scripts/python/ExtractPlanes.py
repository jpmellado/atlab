#!/usr/bin/python3

import numpy as np
import struct
import sys

setofplanes = [1, 2, 3, 4, 5]

sizeofdata = 4  # in bytes
# sizeofdata = 1 # for gate files

sizeofheader = -1  # read it from header in each file
# sizeofheader = 0
# sizeofheader = 36 # for gate files
# sizeofheader = 52 # for flow files
# sizeofheader = 44 # for scal files

# etype = ">" # big-endian
etype = "<"  # little-endian

# do not edit below this line
import atlab

# getting data from stdin
if len(sys.argv) <= 2:
    print("Usage: python $0 [xy,xz,yz] list-of-files.")
    quit()

planetype = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

nx, ny, nz = atlab.getGridSize()  # getting grid size from tlab.ini, if necessary

if planetype == "xy":
    sizeofmask = len(str(nz))
elif planetype == "xz":
    sizeofmask = len(str(ny))
elif planetype == "yz":
    sizeofmask = len(str(nx))
else:
    print("Usage: python $0 [xy,xz,yz] list-of-files.")
    quit()


# further initialization
def tag(sizeofmask, number):
    a = str(number)
    for i in range(sizeofmask - len(a)):
        a = "0" + a
    return a


# extracting planes
for file in setoffiles:
    print("Processing file %s ..." % file)
    if file == "grid":
        x, y, z = atlab.getGrid(nx, ny, nz, "grid")  # get grid as n-tuples

        rawx = struct.pack(etype + "{}f".format(nx), *x)
        rawy = struct.pack(etype + "{}f".format(ny), *y)
        rawz = struct.pack(etype + "{}f".format(nz), *z)

        fout = open("grid." + planetype, "wb")
        if planetype == "xy":
            fout.write(rawx)
            fout.write(rawy)
        elif planetype == "xz":
            fout.write(rawx)
            fout.write(rawz)
        elif planetype == "yz":
            fout.write(rawy)
            fout.write(rawz)
        fout.close()

    else:
        fin = open(file, "rb")
        if sizeofheader == -1:
            raw = fin.read(4)
            sizeofheader = struct.unpack(etype + "{}i".format(1), raw)[0]
            # print(sizeofheader)
            # raw = fin.read(int(sizeofheader-4-8))
            # planeIds = np.array(struct.unpack(etype+'{}i'.format(int((sizeofheader-4-8)/4)), raw))
            # print(planeIds)
            # raw = fin.read(8)
            # time = struct.unpack(etype+'{}d'.format(1), raw)
            # print(time)
        for plane in setofplanes:
            fout = open(file + "." + planetype + tag(sizeofmask, plane), "wb")
            if planetype == "xy":
                fin.seek(sizeofheader + (plane - 1) * nx * ny * sizeofdata, 0)
                raw = fin.read(nx * ny * sizeofdata)
                fout.write(raw)
            elif planetype == "xz":
                fin.seek(sizeofheader + (plane - 1) * nx * sizeofdata, 0)
                for k in range(nz):
                    raw = fin.read(nx * sizeofdata)
                    fout.write(raw)
                    fin.seek((ny - 1) * nx * sizeofdata, 1)
            elif planetype == "yz":
                fin.seek(sizeofheader + (plane - 1) * sizeofdata, 0)
                for jk in range(ny * nz):
                    raw = fin.read(sizeofdata)
                    fout.write(raw)
                    fin.seek((nx - 1) * sizeofdata, 1)

            fout.close()
        fin.close()
