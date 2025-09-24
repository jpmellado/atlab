import numpy as np
import matplotlib.pyplot as plt
import struct
import sys

from matplotlib import rc

rc("text", usetex=True)
rc("text.latex", preamble=r"\usepackage{fourier}")
rc("font", family="serif", size=12)
rc("grid", linestyle="dotted")
rc("axes", grid=True)

# At home, screen 27''
rc("figure", dpi=200)
rc("savefig", dpi=100)

# etype = ">" # big-endian
etype = "<"  # little-endian

# getting data from stdin
if len(sys.argv) < 3:
    print("Usage: python $0 [x,y,z] list-of-grid-files")
    quit()

direction = sys.argv[1]
setoffiles = sorted(sys.argv[2:])

def main():
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
    for file in setoffiles:
        x, y, z = readGrid(file)

        if direction == 'x':
            grid = x
        elif direction == 'y':
            grid = y
        elif direction == 'z':
            grid = z

        n = np.size(grid)
        if n == 1:
            print('Only one grid point in that direction.')
            quit()

        s = np.linspace(1, n, num=n)        # computational, uniform grid

        spacing = np.gradient(grid, edge_order=2)
        axs[0].plot(grid, spacing, label=file)

        stretching = (
            np.gradient(spacing, edge_order=2) / spacing * 100.0
        )  # in percentage
        axs[1].plot(grid, stretching)

    for ax in axs:
        ax.set_ylim([0, None])
        ax.set_xlim([0, None])
        ax.set_xlabel(r"Node position")

        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_position(("axes", -0.01))
        ax.get_yaxis().tick_left()
        ax.spines["top"].set_visible(False)
        ax.spines["bottom"].set_position(("axes", -0.01))
        ax.get_xaxis().tick_bottom()

    axs[0].set_ylabel(r"spacing")
    axs[1].set_ylabel(r"stretching (\%)")
    axs[0].legend()
    plt.tight_layout(pad=0.1)

    plt.savefig("figure1.pdf", bbox_inches="tight")
    plt.show()


def readGrid(filename):
    fin = open(filename, "rb")

    fin.seek(4, 0)
    raw = fin.read(3 * 4)
    nx, ny, nz = struct.unpack(etype + "{}i".format(3), raw)
    print("Grid size is {}x{}x{}.".format(nx, ny, nz))

    # fin.seek(56, 0)
    fin.seek(4+4+8*3+4+4, 1)
    raw = fin.read(nx * 8)
    x = struct.unpack(etype + "{}d".format(nx), raw)

    fin.seek(8, 1)
    raw = fin.read(ny * 8)
    y = struct.unpack(etype + "{}d".format(ny), raw)

    fin.seek(8, 1)
    raw = fin.read(nz * 8)
    z = struct.unpack(etype + "{}d".format(nz), raw)

    fin.close()

    return x, y, z

if __name__ == "__main__":
    main()
