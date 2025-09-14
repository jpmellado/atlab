# Plotting eigenvalues of advection-diffusion equation
# Banded matrices given as lsh and rhs files following atlab structure
#   A du = B u    =>  du = A^{-1} B u 
# operators is then
#   c du = c  A^{-1} B u 
# where c is a constant
import numpy as np
import scipy.linalg
import matplotlib.pyplot as plt
import sys

from matplotlib import rc

rc("text", usetex=True)
rc("text.latex", preamble=r"\usepackage{fourier}")
rc("font", family="serif", size=12)
rc("grid", linestyle="dotted")
rc("axes", grid=True)

if len(sys.argv) not in [4, 7]:
    print("Usage: python $0 file-lhs file-rhs value [file-lhs file-rhs value].")
    quit()


def main():
    # construct matrix system
    A1, B1 = createMatrices(sys.argv[1], sys.argv[2])
    L = float(sys.argv[3]) * scipy.linalg.solve(A1, B1, assume_a="banded")

    if len(sys.argv) == 7:
        A2, B2 = createMatrices(sys.argv[4], sys.argv[5])
        L = L + float(sys.argv[6]) * scipy.linalg.solve(A2, B2, assume_a="banded")

    # obtain eigenvalues
    # lambdas = scipy.linalg.eigvals(L[1:, 1:])  # BC at the beginning of the interval
    # lambdas = scipy.linalg.eigvals(L[:-1, :-1])  # BC at the end of the interval
    lambdas = scipy.linalg.eigvals(L[1:-1, 1:-1])  # BCs at both ends

    # output information
    print("Maximum real part: ", np.max(np.real(lambdas)))

    plt.scatter(np.real(lambdas), np.imag(lambdas))
    plt.xlabel(r"real part")
    plt.ylabel(r"imaginary part")
    plt.savefig("figure1.pdf", bbox_inches="tight")
    plt.show()


def createMatrices(fileLhs, fileRhs):
    # read data
    lhs = np.loadtxt(fileLhs)
    rhs = np.loadtxt(fileRhs)

    nx = np.shape(lhs)[0]
    if np.shape(rhs)[0] != nx:
        print("Inconsistent array sizes.")
        quit()

    ndl = np.shape(lhs)[1]  # number of diagonals
    idl = int(ndl / 2)  # central diagonal
    ndr = np.shape(rhs)[1]
    idr = int(ndr / 2)

    # construct arrays
    A = np.diagflat(lhs[:, idl])
    for ic in range(1, idl + 1):
        A = A + np.diagflat(lhs[:-ic, idl + ic], ic)
        A = A + np.diagflat(lhs[ic:, idl - ic], -ic)
    # print(np.array_str(A, precision=2, suppress_small=True))

    B = np.diagflat(rhs[:, idr])
    for ic in range(1, idr + 1):
        B = B + np.diagflat(rhs[:-ic, idr + ic], ic)
        B = B + np.diagflat(rhs[ic:, idr - ic], -ic)
    # tlab uses extended domains at the boundary
    # by storing the possible additional node in the rhs
    # in the first element
    B[0, idr + 1] = rhs[0, 0]
    B[-1, -idr - 2] = rhs[-1, -1]
    # print(np.array_str(B, precision=2, suppress_small=True))

    return A, B


if __name__ == "__main__":
    main()
