! debug code for levante where we are experiencing some problems possibly due to stripping of disks.
! mpifort -nbs -save-temps -heap-arrays -simd -vec-threshold50 -unroll-aggressive -march=core-avx2 -mtune=core-avx2 -qopt-prefetch -O3 -ipo vmpi_io.f90
program vmpi_io_levante
    use mpi_f08
    implicit none

    type(MPI_File) mpio_fh
    integer mpio_locsize
    type(MPI_Status) status
    integer(KIND=MPI_OFFSET_KIND) offset
    type(MPI_Datatype) subarray
    integer ims_err

    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.0d0)
    ! integer, parameter :: wp = dp               ! working precision is double
    integer, parameter :: wp = sp               ! working precision is single
    integer, parameter :: wi = selected_int_kind(9)

    integer, parameter :: sizeofreal = sizeof(1.0_wp)
    integer, parameter :: sizeofint = sizeof(1_wi)

    ! Decomposition along X and Y in ims_npro_i x ims_npro_j pencils
    integer ims_npro                            ! number of tasks
    ! integer, parameter :: ims_npro_i = 64       ! number of tasks in X
    ! integer, parameter :: ims_npro_j = 64       ! number of tasks in Y
    integer, parameter :: ims_npro_i = 128      ! number of tasks in X
    integer, parameter :: ims_npro_j = 128      ! number of tasks in Y
    integer ims_pro                             ! local task in global communicator
    integer ims_pro_i, ims_pro_j                ! local task in each directional communicator; here used only as offsets

    ! grid points
    integer nx, ny, nz

    ! number of variables
    integer, parameter :: nv = 3
    integer iv

    ! array
    real(wp), allocatable :: a(:, :)

    ! file name
    character(*), parameter :: name = 'test-io.'
    character(len=32) str

    ! ###############################################################
    call MPI_INIT(ims_err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ims_npro, ims_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, ims_pro, ims_err)

    ! initialize grid
    ! nx = 3072
    ! ny = 3072
    ! nz = 768
    nx = 6144
    ny = 6144
    nz = 1536
    nx = nx/ims_npro_i                      ! task-local number of grid points along X
    ny = ny/ims_npro_j                      ! task-local number of grid points along Z

    ims_pro_i = mod(ims_pro, ims_npro_i)    ! MPI offset
    ims_pro_j = ims_pro/ims_npro_i          ! MPI offset

    allocate (a(nx*ny*nz, nv))
    a = 0.0_wp

    select case (wp)
    case (sp)    ! single precision
        subarray = IO_Create_Subarray_XOY(nx, ny, nz, MPI_REAL4)
    case (dp)    ! double precision
        subarray = IO_Create_Subarray_XOY(nx, ny, nz, MPI_REAL8)
    end select
    offset = 0        ! no header
    !offset = 1*sizeofint
    mpio_locsize = nx*ny*nz
    do iv = 1, nv

        write (str, *) iv; str = trim(adjustl(name))//trim(adjustl(str))

        ! an old MPI_IO bug was that if PE 0 did this and not the others, then it hangs...
        if (ims_pro == 0) then
            !   open(55,file=str,status='unknown',form='unformatted',access='stream')
            !   write(55) nx*ny*nz
            !   close(55)
        end if

        call MPI_BARRIER(MPI_COMM_WORLD, ims_err)

        call MPI_FILE_OPEN(MPI_COMM_WORLD, str, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, mpio_fh, ims_err)
        !call MPI_FILE_OPEN(MPI_COMM_WORLD, str, MPI_MODE_WRONLY, MPI_INFO_NULL, mpio_fh, ims_err)

        select case (wp)
        case (sp)    ! single precision
            call MPI_File_set_view(mpio_fh, offset, MPI_REAL4, subarray, 'native', MPI_INFO_NULL, ims_err)
            call MPI_File_write_all(mpio_fh, a(:, iv), mpio_locsize, MPI_REAL4, status, ims_err)
        case (dp)    ! double precision
            call MPI_File_set_view(mpio_fh, offset, MPI_REAL8, subarray, 'native', MPI_INFO_NULL, ims_err)
            call MPI_File_write_all(mpio_fh, a(:, iv), mpio_locsize, MPI_REAL8, status, ims_err)
        end select

        call MPI_File_close(mpio_fh, ims_err)

    end do

    call MPI_FINALIZE(ims_err)

    stop

contains
    function IO_Create_Subarray_XOY(nx, ny, nz, locType) result(locSubarray)
        integer(wi), intent(in) :: nx, ny, nz
        type(MPI_Datatype), intent(in) :: locType

        type(MPI_Datatype) :: locSubarray
        integer, parameter :: ndims = 3
        integer(wi) :: sizes(ndims), locsize(ndims), offset(ndims)

        sizes = [nx*ims_npro_i, ny*ims_npro_j, nz]
        locsize = [nx, ny, nz]
        offset = [nx*ims_pro_i, ny*ims_pro_j, 0]

        call MPI_Type_create_subarray(ndims, sizes, locsize, offset, MPI_ORDER_FORTRAN, locType, locSubarray, ims_err)
        call MPI_Type_commit(locSubarray, ims_err)

    end function IO_Create_Subarray_XOY

end program
