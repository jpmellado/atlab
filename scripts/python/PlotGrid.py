import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import atlab

# getting data from stdin
if len(sys.argv) <= 2:
    print("Usage: python $0 [x,y,z] list-of-grid-files")
    quit()

direction = sys.argv[1]
files = sorted(sys.argv[2:])

if direction not in ["x", "y", "z"]:
    print("Usage: python $0 [x,y,z] list-of-grid-files")
    quit()


###########################################################
def main():
    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))
    axs[0].set_ylabel(r"spacing")
    axs[1].set_ylabel(r"stretching (\%)")

    for file in files:
        grid = atlab.Grid()
        grid.read(filename=file)

        match direction:
            case "x":
                coors = grid.x
            case "y":
                coors = grid.y
            case "z":
                coors = grid.z

        n = np.size(coors)
        if n == 1:
            print("Only one grid point in that direction.")
            quit()

        s = np.linspace(1, n, num=n)  # computational grid, uniform axis

        spacing = np.gradient(coors, edge_order=2)
        sns.lineplot(x=coors, y=spacing, label=file, ax=axs[0])

        stretching = np.gradient(spacing, edge_order=2) / spacing * 100.0  # in percentage
        sns.lineplot(x=coors, y=stretching, ax=axs[1])

    for ax in axs:
        sns.despine(offset=5, ax=ax)
        ax.set_ylim([0, None])
        ax.set_xlim([0, None])
        ax.set_xlabel(r"Node position")
        ax.grid(ls="--")

    axs[0].legend()
    plt.tight_layout(pad=0.1)
    plt.savefig("figure-grid.pdf", bbox_inches="tight")
    plt.show()


###########################################################
if __name__ == "__main__":
    main()
