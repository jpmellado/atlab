#!/usr/bin/python3

import sys
import numpy as np
import atlab

planes = [1, 2, 3, 4, 5]

dtype = "d"  # floating-point number, double precision
# dtype = "f"  # floating-point number, single precision
# dtype = "B"  # unsigned character, for gate files

# getting data from stdin
if len(sys.argv) <= 2:
    print("Usage: python $0 [xy,xz,yz] list-of-files.")
    quit()

planetype = sys.argv[1]
files = sorted(sys.argv[2:])

if planetype not in ["xy", "xz", "yz"]:
    print("Usage: python $0 [xy,xz,yz] list-of-files.")
    quit()

###########################################################
nx, ny, nz = atlab.getGridSize()  # getting grid size from tlab.ini

match planetype:
    case "xy":
        sizeofmask = len(str(nz))
    case "xz":
        sizeofmask = len(str(ny))
    case "yz":
        sizeofmask = len(str(nx))


# further initialization
def tag(sizeofmask, number):
    a = str(number)
    for i in range(sizeofmask - len(a)):
        a = "0" + a
    return a


###########################################################
# extracting planes
for file in files:
    print("Processing file %s ..." % file)

    field = atlab.Field(file, dtype=dtype)
    for plane in planes:
        fout = open(file + "." + planetype + tag(sizeofmask, plane), "wb")

        match planetype:
            case "xy":
                raw = field.readXYRaw(plane)

            case "xz":
                raw = field.readXZRaw(plane)

            case "yz":
                raw = field.readYZRaw(plane)

        fout.write(raw)
        fout.close()

    field.fin.close()
