#!/usr/bin/python3

import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import atlab

# getting data from stdin
if len(sys.argv) <= 2:
    print("Usage: python $0 [xy,xz,yz] list-of-files.")
    quit()

planetype = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

if planetype not in ["xy", "xz", "yz"]:
    print("Usage: python $0 [xy,xz,yz] list-of-files.")
    quit()

###########################################################
# process grid
grid = atlab.Grid()
grid.read()

match planetype:
    case "xy":
        x1 = grid.x
        x2 = grid.y
    case "xz":
        x1 = grid.x
        x2 = grid.z
    case "yz":
        x1 = grid.y
        x2 = grid.z


nx1 = len(x1)
nx2 = len(x2)

###########################################################
# process files
for file in setoffiles:
    print("Processing file %s ..." % file)

    field = atlab.Field(file, nx=nx1, ny=nx2, nz=1, sizeofheader=0, dtype="d")
    a = field.read()
    a = a.reshape((nx2, nx1))
    field.fin.close()

    mean, std = np.mean(a), np.std(a)
    plt.figure(figsize=(10, 8))
    plt.pcolormesh(x1, x2, a, shading="auto", vmin=mean - std, vmax=mean + std)
    # plt.contourf(x1,x2,a)
    plt.axis("equal")
    # plt.gca().set_xlim([x1[0],x1[-1]])
    # plt.gca().set_ylim([x2[0],x2[-1]])
    # plt.axis([ x1[0], x1[-1], x2[0], x2[-1]])
    plt.colorbar()
    plt.tight_layout(pad=2)
    plt.title(file)
    # plt.savefig("{}.jpg".format(file),dpi=150,bbox_inches='tight')
    plt.show()
    plt.close()
