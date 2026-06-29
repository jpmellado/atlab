import struct
import numpy as np


###########################################################
###########################################################
class Grid:
    def read(
        self,
        filename="grid",
        etype="<",  # little-endian
        # etype=">",  # big endian
    ):
        fin = open(filename, "rb")

        # grid file contains 4 bytes record lengths between records
        # total header size is 56 bytes
        fin.seek(4, 0)
        raw = fin.read(3 * 4)
        nx, ny, nz = struct.unpack(etype + "{}i".format(3), raw)
        print("Grid size is {}x{}x{}.".format(nx, ny, nz))

        fin.seek(4 + 4, 1)
        fin.seek(8 * 3, 1)  # advance over scalex, scaley, scalez

        fin.seek(4 + 4, 1)
        raw = fin.read(nx * 8)
        self.x = struct.unpack(etype + "{}d".format(nx), raw)

        fin.seek(4 + 4, 1)
        raw = fin.read(ny * 8)
        self.y = struct.unpack(etype + "{}d".format(ny), raw)

        fin.seek(4 + 4, 1)
        raw = fin.read(nz * 8)
        self.z = struct.unpack(etype + "{}d".format(nz), raw)

        fin.close()

    def writeNetCDF(filename):
        print("to be done.")


# getting grid size from tlab.ini in case grid file is not available
def getGridSize(nx=0, ny=0, nz=0, filename="tlab.ini"):
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


###########################################################
###########################################################
class Field:
    def __init__(
        self,
        filename,
        nx=0,
        ny=0,
        nz=0,
        sizeofheader=-1,  # read it from header in each file
        # sizeofheader=0,  # bytes
        # sizeofheader=36,  # bytes, for gate files
        # sizeofheader=52,  # bytes, for flow files
        # sizeofheader=44,  # bytes, for scal files
        etype="<",  # little-endian
        # etype=">",  # big endian
        dtype="f",  # floating-point number, single precision
        # dtype="f",  # floating-point number, double precision
        # dtype="B",  # unsigned character, for gate files
    ):

        self.fin = open(filename, "rb")

        self.sizeofheader = sizeofheader
        if sizeofheader == -1:  # read header
            raw = self.fin.read(4)
            self.sizeofheader = struct.unpack(etype + "{}i".format(1), raw)[0]
            raw = self.fin.read(4)
            self.nx = struct.unpack(etype + "{}i".format(1), raw)[0]
            raw = self.fin.read(4)
            self.ny = struct.unpack(etype + "{}i".format(1), raw)[0]
            raw = self.fin.read(4)
            self.nz = struct.unpack(etype + "{}i".format(1), raw)[0]
            raw = self.fin.read(4)
            self.nt = struct.unpack(etype + "{}i".format(1), raw)[0]
            raw = self.fin.read(8)
            self.time = struct.unpack(etype + "{}d".format(1), raw)[0]
        else:
            self.nx = nx
            self.ny = ny
            self.nz = nz

        self.dtype = dtype
        match dtype:
            case "f":  # signle precision
                self.sizeofdata = 4  # in bytes
            case "d":  # double precision
                self.sizeofdata = 8
            case "B":  # unsigned character, for gate files
                self.sizeofdata = 1  # in bytes

        self.etype = etype

        self.sizexy = self.nx * self.ny * self.sizeofdata
        self.sizex = self.nx * self.sizeofdata
        # self.sizexz = nx * nz * self.sizeofdata
        # self.sizeyz = ny * nz * self.sizeofdata

        return

    def readRaw(self):
        self.fin.seek(self.sizeofheader, 0)
        # self.fin.seek(-nx * ny * nz * sizeofdata, 2)  # read the field
        raw = self.fin.read()
        return raw

    def read(self):
        a = struct.unpack(
            (self.etype + "{}" + self.dtype).format(int(self.nx * self.ny * self.nz)),
            self.readRaw(),
        )
        return np.array(a).reshape((self.nz, self.ny, self.nx))

    def readXYRaw(self, plane):
        self.fin.seek(self.sizeofheader + (plane - 1) * self.sizexy, 0)
        raw = self.fin.read(self.sizexy)
        return raw

    def readXY(self, plane):
        a = struct.unpack(
            (self.etype + "{}" + self.dtype).format(int(self.nx * self.ny)), self.readXYRaw(plane)
        )

        return np.array(a).reshape((self.ny, self.nx))

    def readXZRaw(self, plane):
        self.fin.seek(self.sizeofheader + (plane - 1) * self.sizex, 0)
        raw = b""
        for k in range(self.nz):
            raw = raw + self.fin.read(self.sizex)
            self.fin.seek((self.ny - 1) * self.sizex, 1)
        return raw

    def readXZ(self, plane):
        a = struct.unpack(
            (self.etype + "{}" + self.dtype).format(int(self.nx * self.nz)), self.readXYRaw(plane)
        )

        return np.array(a).reshape((self.nz, self.nx))

    def readYZRaw(self, plane):
        self.fin.seek(self.sizeofheader + (plane - 1) * self.sizeofdata, 0)
        raw = b""
        for jk in range(self.ny * self.nz):
            raw = raw + self.fin.read(self.sizeofdata)
            self.fin.seek((self.nx - 1) * self.sizeofdata, 1)
        return raw

    def readYZ(self, plane):
        a = struct.unpack(
            (self.etype + "{}" + self.dtype).format(int(self.ny * self.nz)), self.readYZRaw(plane)
        )

        return np.array(a).reshape((self.nz, self.ny))


# writing minimal netcfd file, single precision
def writeNetCDF(filename, varname, time, x, y, z, field):
    import netCDF4 as nc

    # creating netcdf
    file_dst = nc.Dataset(filename + ".nc", "w")

    # create dimensions
    file_dst.createDimension("t", None)
    file_dst.createDimension("x", len(x))
    file_dst.createDimension("y", len(y))
    file_dst.createDimension("z", len(z))

    # create and write independent variables
    t_dst = file_dst.createVariable("t", "f4", ("t",))
    x_dst = file_dst.createVariable("x", "f4", ("x",))
    y_dst = file_dst.createVariable("y", "f4", ("y",))
    z_dst = file_dst.createVariable("z", "f4", ("z",))
    t_dst[:] = time
    x_dst[:] = x[:]
    y_dst[:] = y[:]
    z_dst[:] = z[:]

    # create and write field
    var_dst = file_dst.createVariable(
        varname,
        "f4",
        (
            "t",
            "z",
            "y",
            "x",
        ),
    )
    var_dst[0, :, :, :] = field.reshape((len(z), len(y), len(x)))

    file_dst.close()
