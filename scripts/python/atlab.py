import struct

if __name__ != "__main__":
    from __main__ import etype

# getting grid size from tlab.ini
def getGridSize(nx=0, ny=0, nz=0):
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

    return nx, ny, nz

# reading grid file and return n-tuples
def getGrid(nx, ny, nz, name='grid'):
    fin = open(name, 'rb')
    #
    fin.seek(56,0)
    raw = fin.read(nx*8)
    x = struct.unpack(etype+'{}d'.format(nx), raw)
    #
    fin.seek(8,1)
    raw = fin.read(ny*8)
    y = struct.unpack(etype+'{}d'.format(ny), raw)
    #
    fin.seek(8,1)
    raw = fin.read(nz*8)
    z = struct.unpack(etype+'{}d'.format(nz), raw)
    #
    fin.close()

    return x, y, z

# writing minimal netcfd file, single precision
def writeNetCfd(filename, varname, time, x, y, z, field):
    import netCDF4 as nc

    # creating netcdf
    file_dst = nc.Dataset(filename+'.nc', 'w')

    # create dimensions
    file_dst.createDimension('t',None)
    file_dst.createDimension('x',len(x))
    file_dst.createDimension('y',len(y))
    file_dst.createDimension('z',len(z))

    # create and write independent variables 
    t_dst = file_dst.createVariable('t', 'f4', ('t',))
    x_dst = file_dst.createVariable('x', 'f4', ('x',))
    y_dst = file_dst.createVariable('y', 'f4', ('y',))
    z_dst = file_dst.createVariable('z', 'f4', ('z',))
    t_dst[:] = time
    x_dst[:] = x[:]
    y_dst[:] = y[:]
    z_dst[:] = z[:]

    # create and write field
    var_dst = file_dst.createVariable(varname, 'f4', ('z','y','x',))
    var_dst[:] = field.reshape((len(z),len(y),len(x)))[:]

    file_dst.close()
 

