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
if len(sys.argv) < 2:
    print("Usage: python $0 list-of-grid-files")
    quit()

setoffiles = sorted(sys.argv[1:])


def main():
    global nx, ny, nz
    nx, ny, nz = readGridSize("tlab.ini")

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
    for file in setoffiles:
        s = np.linspace(1, nz, num=nz)
        x, y, z = readGrid(file)

        spacing = np.gradient(z, edge_order=2)
        axs[0].plot(z, spacing, label=file)

        stretching = (
            np.gradient(spacing, edge_order=2) / spacing * 100.0
        )  # in percentage
        axs[1].plot(z, stretching)

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

    fin.seek(56, 0)
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


def readGridSize(filename, nx=0, ny=0, nz=0):
    if nx == 0:
        for line in open(filename):
            if line.lower().replace(" ", "").startswith("imax="):
                nx = int(line.split("=", 1)[1])
                break

    if ny == 0:
        for line in open(filename):
            if line.lower().replace(" ", "").startswith("jmax="):
                ny = int(line.split("=", 1)[1])
                break

    if nz == 0:
        for line in open(filename):
            if line.lower().replace(" ", "").startswith("kmax="):
                nz = int(line.split("=", 1)[1])
                break

    print("Grid size is {}x{}x{}.".format(nx, ny, nz))

    return nx, ny, nz


if __name__ == "__main__":
    main()
