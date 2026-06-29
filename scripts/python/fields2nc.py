#!/usr/bin/python3

import sys
import numpy as np
import atlab

dtype = "d"  # floating-point number, double precision
# dtype = "f"  # floating-point number, single precision
# dtype = "B"  # unsigned character, for gate files

# getting data from stdin
if len(sys.argv) <= 3:
    print("Usage: python $0 [3d,xy,xz,yz] varname list-of-files.")
    quit()

datatype = sys.argv[1]
varname = sys.argv[2]
files = sorted(sys.argv[3:])

if datatype not in ["3d", "xy", "xz", "yz"]:
    print("Usage: python $0 [3d,xy,xz,yz] list-of-files.")
    quit()

###########################################################
# process grid
grid = atlab.Grid()
grid.read()

x, y, z = np.copy(grid.x), np.copy(grid.y), np.copy(grid.z)
nx, ny, nz = np.size(x), np.size(y), np.size(z)

# handle grid data in case of planes
match datatype:
    case "xy":
        nz = 1
        z = np.array([0.0], dtype=np.float32)
    case "xz":
        ny = 1
        y = np.array([0.0], dtype=np.float32)
    case "yz":
        nx = 1
        x = np.array([0.0], dtype=np.float32)

###########################################################
# process files
for file in files:
    print("Processing file %s ..." % file)

    field = atlab.Field(file, dtype=dtype)
    # print(field.nx, field.ny, field.nz, field.nt, field.time)
    a = field.read()

    atlab.writeNetCDF(file, varname, field.time, x, y, z, a)

    field.fin.close()
