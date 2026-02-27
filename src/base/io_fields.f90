#include "tlab_error.h"

#define USE_ACCESS_STREAM

!########################################################################
!#
!# Optional header/metadata
!# Unformatted records
!# No embedded record information
!#
!########################################################################

module IO_Fields
    use TLab_Constants, only: wp, wi, sp, dp, sizeofint, sizeofreal
    use TLab_Constants, only: lfile, wfile, efile
    use TLab_Constants, only: MAX_PARS, MAX_VARS
    use TLab_WorkFlow, only: TLab_Stop, TLab_Write_ASCII
    use TLab_Arrays, only: wrk3d
    use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_VARS
    use TLabMPI_PROCS, only: TLabMPI_Panic
#endif
    implicit none
    private

    public :: IO_Initialize
    public :: IO_Read_Fields, IO_Write_Fields
    public :: IO_Read_Field_INT1, IO_Write_Field_INT1

    public :: IO_TYPE_DOUBLE, IO_TYPE_SINGLE
    public :: io_header_q, io_header_s

    public :: io_subarray_dt
    public :: IO_Read_Subarray, IO_Write_Subarray
#ifdef USE_MPI
    public :: IO_Create_Subarray_XOY, IO_Create_Subarray_XOZ, IO_Create_Subarray_YOZ
#endif

    ! -------------------------------------------------------------------
    integer :: io_fileformat                        ! files format
    integer, parameter :: IO_MPIIO = 1
    integer, parameter :: IO_NETCDF = 2
    integer, parameter :: IO_NOFILE = 3

    integer :: io_datatype                          ! single or double precision
    integer, parameter :: IO_TYPE_SINGLE = 1
    integer, parameter :: IO_TYPE_DOUBLE = 2

    type :: io_header                               ! header information
        sequence
        integer size
        real(wp) params(MAX_PARS)
    end type
    type(io_header) :: io_header_q(1), io_header_s(MAX_VARS)

    type :: io_subarray_dt                          ! subarray types
        ! sequence
        integer :: precision = IO_TYPE_DOUBLE
#ifdef USE_MPI
        logical active, lpadding(3)
        type(MPI_Comm) communicator
        type(MPI_Datatype) subarray
        integer(KIND=MPI_OFFSET_KIND) offset
#else
        integer offset
#endif
    contains
        procedure :: read_double => io_subarray_read_double
        procedure :: read_single => io_subarray_read_single
        procedure :: read_int1 => io_subarray_read_int1
        generic :: read => read_double, read_single, read_int1
        !
        procedure :: write_double => io_subarray_write_double
        procedure :: write_single => io_subarray_write_single
        procedure :: write_int1 => io_subarray_write_int1
        generic :: write => write_double, write_single, write_int1
    end type io_subarray_dt

    ! -------------------------------------------------------------------
    type(io_subarray_dt) :: io_subarray_main
    character(len=64) name
    real(sp), pointer :: s_wrk(:) => null()

#ifdef USE_MPI
    type(MPI_File) mpio_fh
    integer mpio_locsize
    type(MPI_Status) status
#endif

contains
    !########################################################################
    !########################################################################
    subroutine IO_Initialize(locInifile)
        use TLab_Constants, only: ifile
        character(len=*), optional, intent(in) :: locInifile

        ! -------------------------------------------------------------------
        character(len=32) inifile
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        ! ###################################################################
        ! read
        if (present(locInifile)) then
            inifile = locInifile
        else
            inifile = ifile
        end if

        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'WorkFlow'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#FileFormat=<mpiio/NetCDF/None>')
        call TLab_Write_ASCII(bakfile, '#FileDatatype=<Double/Single>')

        call ScanFile_Char(bakfile, inifile, block, 'FileFormat', 'MpiIO', sRes)
        if (trim(adjustl(sRes)) == 'mpiio') then; io_fileformat = IO_MPIIO
        elseif (trim(adjustl(sRes)) == 'netcdf') then; io_fileformat = IO_NETCDF
        elseif (trim(adjustl(sRes)) == 'none') then; io_fileformat = IO_NOFILE
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong FileFormat.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'FileDatatype', 'Double', sRes)
        if (trim(adjustl(sRes)) == 'double') then; io_datatype = IO_TYPE_DOUBLE
        elseif (trim(adjustl(sRes)) == 'single') then; io_datatype = IO_TYPE_SINGLE
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong FileDatatype.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        ! ###################################################################
        ! initialize
        io_subarray_main%precision = io_datatype
#ifdef USE_MPI
        io_subarray_main%communicator = MPI_COMM_WORLD
#endif

        return
    end subroutine

#ifdef USE_MPI
    !########################################################################
    !########################################################################
    function IO_Create_Subarray_XOZ(nx, nz, locType) result(locSubarray)
        integer(wi), intent(in) :: nx, nz
        type(MPI_Datatype), intent(in) :: locType

        type(MPI_Datatype) :: locSubarray
        integer, parameter :: ndims = 2
        integer(wi) :: sizes(ndims), locsize(ndims), offset(ndims)

        sizes = [nx*ims_npro_i, nz]
        locsize = [nx, nz]
        offset = [nx*ims_pro_i, 0]

        call MPI_Type_create_subarray(ndims, sizes, locsize, offset, MPI_ORDER_FORTRAN, locType, locSubarray, ims_err)
        call MPI_Type_commit(locSubarray, ims_err)

    end function IO_Create_Subarray_XOZ

    !########################################################################
    !########################################################################
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

    !########################################################################
    !########################################################################
    function IO_Create_Subarray_YOZ(ny, nz, locType) result(locSubarray)
        integer(wi), intent(in) :: ny, nz
        type(MPI_Datatype), intent(in) :: locType

        type(MPI_Datatype) :: locSubarray
        integer, parameter :: ndims = 2
        integer(wi) :: sizes(ndims), locsize(ndims), offset(ndims)

        sizes = [ny*ims_npro_j, nz]
        locsize = [ny, nz]
        offset = [ny*ims_pro_j, 0]

        call MPI_Type_create_subarray(ndims, sizes, locsize, offset, MPI_ORDER_FORTRAN, locType, locSubarray, ims_err)
        call MPI_Type_commit(locSubarray, ims_err)

    end function IO_Create_Subarray_YOZ
#endif

    !########################################################################
    !########################################################################
    !# Read/write files of size (nx*ims_npro_i)x(ny*ims_npro_y)x(nz*ims_npro_z)
    subroutine IO_Read_Fields(fname, nx, ny, nz, nt, nfield, iread, a, params)
        character(LEN=*) fname
        integer, intent(in) :: nfield, iread   ! iread=0 reads all nfields, otherwise iread field
        integer(wi), intent(in) :: nx, ny, nz, nt
        real(wp), intent(out) :: a(nx*ny*nz, *)
        real(wp), intent(inout) :: params(:)

        ! -------------------------------------------------------------------
        integer(wi) header_offset, isize
        integer ifield, iz

        ! ###################################################################
        if (io_datatype == IO_TYPE_SINGLE) then
            call Write_Log('Reading single precision field '//trim(adjustl(fname)), nx, ny, nz)
            ! Pass memory address from double precision array to single precision array
            call c_f_pointer(c_loc(wrk3d), s_wrk, shape=[nx*ny*nz])
        else
            call Write_Log('Reading double precision field '//trim(adjustl(fname)), nx, ny, nz)
        end if

        ! ###################################################################
        select case (io_fileformat)

        case (IO_NOFILE)         ! Do nothing
            a(:, 1:nfield) = 0.0_wp

        case (IO_NETCDF)         ! To be implemented

        case DEFAULT              ! One file with header per field
#ifdef USE_MPI
            if (io_subarray_main%precision == IO_TYPE_SINGLE) then
                io_subarray_main%subarray = IO_Create_Subarray_XOY(nx, ny, nz, MPI_REAL4)
            else
                io_subarray_main%subarray = IO_Create_Subarray_XOY(nx, ny, nz, MPI_REAL8)
            end if
#endif

            ! -------------------------------------------------------------------
            ! read data
            iz = 0
            do ifield = 1, nfield
                if (iread == 0 .or. iread == ifield) then
                    iz = iz + 1
                    write (name, '(I2)') ifield
                    name = trim(adjustl(fname))//'.'//trim(adjustl(name))

                    ! -------------------------------------------------------------------
                    ! header
                    call Read_Header(header_offset, nx, ny, nz, nt, params)

                    ! -------------------------------------------------------------------
                    ! field
                    io_subarray_main%offset = header_offset
                    if (io_datatype == IO_TYPE_SINGLE) then
                        call io_subarray_main%read(s_wrk(:))
                        a(:, iz) = real(s_wrk(:), dp)
                    else
                        call io_subarray_main%read(a(:, iz))
                    end if

                end if
            end do

            ! -------------------------------------------------------------------
            ! process header info
            isize = (header_offset - 5*SIZEOFINT)/SIZEOFREAL ! Size of array params
#ifdef USE_MPI
            if (isize > 0) call MPI_BCAST(params, size(params), MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
#endif

        end select

        return
    end subroutine IO_Read_Fields

    !########################################################################
    !########################################################################
    subroutine IO_Read_Field_INT1(name, nx, ny, nz, nt, a, params)
        character(len=*) name
        integer(wi), intent(in) :: nx, ny, nz, nt
        integer(1), intent(out) :: a(nx*ny*nz)
        real(wp), intent(inout) :: params(:)

        ! -------------------------------------------------------------------
        integer(wi) header_offset, isize

        ! ###################################################################
        call Write_Log('Reading integer field '//trim(adjustl(name)), nx, ny, nz)

#ifdef USE_MPI
        io_subarray_main%subarray = IO_Create_Subarray_XOY(nx, ny, nz, MPI_INTEGER1)
#endif

        ! -------------------------------------------------------------------
        ! header
        call Read_Header(header_offset, nx, ny, nz, nt, params)

        ! -------------------------------------------------------------------
        ! field
        io_subarray_main%offset = header_offset
        call io_subarray_main%read(a(:))

        ! -------------------------------------------------------------------
        ! process header info
        isize = (header_offset - 5*SIZEOFINT)/SIZEOFREAL ! Size of array params
#ifdef USE_MPI
        if (isize > 0) call MPI_BCAST(params, size(params), MPI_REAL8, 0, MPI_COMM_WORLD, ims_err)
#endif

        return
    end subroutine IO_Read_Field_INT1

    !########################################################################
    !########################################################################
    ! Arbitrary list of fields, without header
    subroutine IO_Read_Subarray(locSubarrayPlan, fname, varname, data, sizes)
        type(io_subarray_dt), intent(in) :: locSubarrayPlan
        character(len=*), intent(in) :: fname
        integer(wi), intent(in) :: sizes(5) ! total size, lower bound, upper bound, stride, # variables
        character(len=*), intent(in) :: varname(sizes(5))
        real(wp), intent(out) :: data(sizes(1), sizes(5))

        ! -----------------------------------------------------------------------
        integer(wi) iv, isize
        character(len=128) line

        ! #######################################################################
        isize = (sizes(3) - sizes(2))/sizes(4) + 1

        if (locSubarrayPlan%precision == IO_TYPE_SINGLE) then
            line = 'Reading single precision field'
            call c_f_pointer(c_loc(wrk3d), s_wrk, shape=[isize])
        else
            line = 'Reading double precision field'
        end if

#ifdef USE_MPI
        if (locSubarrayPlan%active) then
#endif

            do iv = 1, sizes(5)
                name = trim(adjustl(fname))
                if (varname(iv) /= '') name = trim(adjustl(fname))//'.'//trim(adjustl(varname(iv)))
                call TLab_Write_ASCII(lfile, trim(adjustl(line))//' '//trim(adjustl(name))//'...')

                if (locSubarrayPlan%precision == IO_TYPE_SINGLE) then
                    call locSubarrayPlan%read(s_wrk(:))
                    data(sizes(2):sizes(3):sizes(4), iv) = real(s_wrk(1:isize), wp)
                else
                    call locSubarrayPlan%read(wrk3d(1:isize))
                    data(sizes(2):sizes(3):sizes(4), iv) = wrk3d(1:isize)
                end if

            end do

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine IO_Read_Subarray

    !########################################################################
    !########################################################################
    subroutine IO_Write_Fields(fname, nx, ny, nz, nt, nfield, a, locHeader)
        character(len=*), intent(in) :: fname
        integer, intent(in) :: nfield
        integer(wi), intent(in) :: nx, ny, nz, nt
        real(wp), intent(in) :: a(nx*ny*nz, nfield)
        type(io_header), intent(in), optional :: locHeader(:)

        ! -------------------------------------------------------------------
        integer(wi) header_offset
        integer ifield, ih

        ! ###################################################################
        if (io_datatype == IO_TYPE_SINGLE) then
            call Write_Log('Writing single precision field '//trim(adjustl(fname)), nx, ny, nz)
            ! Pass memory address from double precision array to single precision array
            call c_f_pointer(c_loc(wrk3d), s_wrk, shape=[nx*ny*nz])
        else
            call Write_Log('Writing double precision field '//trim(adjustl(fname)), nx, ny, nz)
        end if

        ! ###################################################################
        select case (io_fileformat)

        case (IO_NOFILE)         ! Do nothing

        case (IO_NETCDF)         ! To be implemented

        case DEFAULT              ! One file with header per field
#ifdef USE_MPI
            if (io_subarray_main%precision == IO_TYPE_SINGLE) then
                io_subarray_main%subarray = IO_Create_Subarray_XOY(nx, ny, nz, MPI_REAL4)
            else
                io_subarray_main%subarray = IO_Create_Subarray_XOY(nx, ny, nz, MPI_REAL8)
            end if
#endif

            do ifield = 1, nfield
                write (name, '(I2)') ifield
                name = trim(adjustl(fname))//'.'//trim(adjustl(name))

                ! -------------------------------------------------------------------
                ! header
                header_offset = 5*SIZEOFINT
                if (present(locHeader)) then
                    ih = min(size(locHeader), ifield)      ! use always the 1 if you enter with one, or the ifield
                    header_offset = header_offset + locHeader(ih)%size*SIZEOFREAL
                end if

                if (present(locHeader)) then
                    call Write_Header(nx, ny, nz, nt, locHeader(ih)%params(1:locHeader(ih)%size))
                else
                    call Write_Header(nx, ny, nz, nt)
                end if

                ! -------------------------------------------------------------------
                ! field
                io_subarray_main%offset = header_offset
                if (io_datatype == IO_TYPE_SINGLE) then
                    s_wrk(:) = real(a(:, ifield), sp)
                    call io_subarray_main%write(s_wrk(:))
                else
                    call io_subarray_main%write(a(:, ifield))
                end if

            end do

        end select

        return
    end subroutine IO_Write_Fields

    !########################################################################
    !########################################################################
    subroutine IO_Write_Field_INT1(name, nx, ny, nz, nt, a, params)
        character(len=*) name
        integer(wi), intent(in) :: nx, ny, nz, nt
        integer(1), intent(in) :: a(nx*ny*nz)
        real(wp), intent(in), optional :: params(:)

        integer(wi) header_offset

        ! ###################################################################
        call Write_Log('Writing integer field '//trim(adjustl(name)), nx, ny, nz)

        ! ###################################################################
#ifdef USE_MPI
        io_subarray_main%subarray = IO_Create_Subarray_XOY(nx, ny, nz, MPI_INTEGER1)
#endif

        ! -------------------------------------------------------------------
        ! header
        header_offset = 5*SIZEOFINT
        if (present(params)) then
            header_offset = header_offset + size(params)*SIZEOFREAL
        end if

        call Write_Header(nx, ny, nz, nt, params(:))

        ! -------------------------------------------------------------------
        ! field
        io_subarray_main%offset = header_offset
        call io_subarray_main%write(a(:))

        return
    end subroutine IO_Write_Field_INT1

    !########################################################################
    !########################################################################
    ! Arbitrary list of fields, without header
    subroutine IO_Write_Subarray(locSubarrayPlan, fname, varname, data, sizes)
        type(io_subarray_dt), intent(in) :: locSubarrayPlan
        character(len=*), intent(in) :: fname
        integer(wi), intent(in) :: sizes(5) ! total size, lower bound, upper bound, stride, # variables
        character(len=*), intent(in) :: varname(sizes(5))
        real(wp), intent(in) :: data(sizes(1), sizes(5))

        ! -----------------------------------------------------------------------
        integer(wi) iv, isize
        character(len=128) line

        ! #######################################################################
        isize = (sizes(3) - sizes(2))/sizes(4) + 1

        if (locSubarrayPlan%precision == IO_TYPE_SINGLE) then
            line = 'Writing single precision field'
            call c_f_pointer(c_loc(wrk3d), s_wrk, shape=[isize])
        else
            line = 'Writing double precision field'
        end if

#ifdef USE_MPI
        if (locSubarrayPlan%active) then
#endif

            do iv = 1, sizes(5)
                name = trim(adjustl(fname))
                if (varname(iv) /= '') name = trim(adjustl(fname))//'.'//trim(adjustl(varname(iv)))
                call TLab_Write_ASCII(lfile, trim(adjustl(line))//' '//trim(adjustl(name))//'...')

                if (locSubarrayPlan%precision == IO_TYPE_SINGLE) then
                    s_wrk(1:isize) = real(data(sizes(2):sizes(3):sizes(4), iv), sp)
                    call locSubarrayPlan%write(s_wrk(:))
                else
                    wrk3d(1:isize) = data(sizes(2):sizes(3):sizes(4), iv)
                    call locSubarrayPlan%write(wrk3d(1:isize))
                end if

            end do

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine IO_Write_Subarray

    !########################################################################
    !########################################################################
    subroutine Write_Log(tag, nx, ny, nz)
        character(len=*), intent(in) :: tag
        integer, intent(in) :: nx, ny, nz

        character(len=128) line
        character(len=32) str

        integer nx_total, ny_total, nz_total

#ifdef USE_MPI
        nx_total = nx*ims_npro_i
        ny_total = ny*ims_npro_j
        nz_total = nz
#else
        nx_total = nx
        ny_total = ny
        nz_total = nz
#endif

        line = trim(adjustl(tag))//' of size'
        write (str, *) nx_total; line = trim(adjustl(line))//' '//trim(adjustl(str))
        write (str, *) ny_total; line = trim(adjustl(line))//'x'//trim(adjustl(str))
        write (str, *) nz_total; line = trim(adjustl(line))//'x'//trim(adjustl(str))//'...'
        call TLab_Write_ASCII(lfile, line)

        return
    end subroutine

!########################################################################
!########################################################################
#define LOC_UNIT_ID 54
#define LOC_STATUS 'old'

    subroutine io_subarray_read_single(self, field)
        class(io_subarray_dt) self
        real(sp), intent(out) :: field(:)

#ifdef USE_MPI
        integer(wi) isize
        isize = size(field)
        call MPI_File_open(self%communicator, trim(adjustl(name)), &
                           MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_set_view(mpio_fh, self%offset, MPI_REAL4, self%subarray, 'native', MPI_INFO_NULL, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_read_all(mpio_fh, field, isize, MPI_REAL4, status, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_close(mpio_fh, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
#else
#include "tlab_open_file.h"
        read (LOC_UNIT_ID, POS=self%offset + 1) field(:)
        close (LOC_UNIT_ID)
#endif

        return
    end subroutine

!########################################################################
!########################################################################
    subroutine io_subarray_read_double(self, field)
        class(io_subarray_dt) self
        real(dp), intent(out) :: field(:)

#ifdef USE_MPI
        integer(wi) isize
        isize = size(field)
        call MPI_File_open(self%communicator, trim(adjustl(name)), &
                           MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_set_view(mpio_fh, self%offset, MPI_REAL8, self%subarray, 'native', MPI_INFO_NULL, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_read_all(mpio_fh, field, isize, MPI_REAL8, status, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_close(mpio_fh, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
#else
#include "tlab_open_file.h"
        read (LOC_UNIT_ID, POS=self%offset + 1) field(:)
        close (LOC_UNIT_ID)
#endif

        return
    end subroutine

!########################################################################
!########################################################################
    subroutine io_subarray_read_int1(self, field)
        class(io_subarray_dt) self
        integer(1), intent(out) :: field(:)

#ifdef USE_MPI
        integer(wi) isize
        isize = size(field)
        call MPI_File_open(self%communicator, trim(adjustl(name)), &
                           MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_set_view(mpio_fh, self%offset, MPI_INTEGER1, self%subarray, 'native', MPI_INFO_NULL, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_read_all(mpio_fh, field, isize, MPI_INTEGER1, status, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_close(mpio_fh, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
#else
#include "tlab_open_file.h"
        read (LOC_UNIT_ID, POS=self%offset + 1) field(:)
        close (LOC_UNIT_ID)
#endif
        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine Read_Header(offset, nx, ny, nz, nt, params)
        integer(wi), intent(out) :: offset
        integer(wi), intent(in) :: nx, ny, nz, nt
        real(wp), intent(inout) :: params(:)

        ! -------------------------------------------------------------------
        integer isize
        integer(wi) nx_loc, ny_loc, nz_loc, nt_loc
        integer nx_total, ny_total, nz_total

        !########################################################################
#ifdef USE_MPI
        if (ims_pro == 0) then
            nx_total = nx*ims_npro_i
            ny_total = ny*ims_npro_j
            nz_total = nz
#else
            nx_total = nx
            ny_total = ny
            nz_total = nz
#endif
#include "tlab_open_file.h"
            rewind (LOC_UNIT_ID)

            read (LOC_UNIT_ID) offset, nx_loc, ny_loc, nz_loc, nt_loc

            ! Check
            if (any([nx_total, ny_total, nz_total] /= [nx_loc, ny_loc, nz_loc])) then
                call TLab_Write_ASCII(wfile, __FILE__//'. Grid size mismatch.')
                ! call TLab_Write_ASCII(efile, __FILE__//'. Grid size mismatch.')
                ! call TLab_Stop(DNS_ERROR_DIMGRID)
            end if

            if (nt /= nt_loc) then
                call TLab_Write_ASCII(wfile, __FILE__//'. ItNumber mismatch. Filename value ignored.')
                ! nt = nt_loc
            end if

            isize = offset - 5*SIZEOFINT
            if (isize >= 0 .and. mod(isize, SIZEOFREAL) == 0) then
                ! isize = isize/SIZEOFREAL
                read (LOC_UNIT_ID) params(:)
                ! elseif (isize == 0) then
                !     continue ! no params to read; header format is correct
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Header format incorrect.')
                call TLab_Stop(DNS_ERROR_RECLEN)
            end if

#ifdef USE_MPI
        end if
        call MPI_BCAST(offset, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ims_err)
#endif

        return
    end subroutine Read_Header

#undef LOC_UNIT_ID
#undef LOC_STATUS

!########################################################################
!########################################################################
#define LOC_UNIT_ID 55
#define LOC_STATUS 'unknown'

    subroutine io_subarray_write_single(self, field)
        class(io_subarray_dt) self
        real(sp), intent(in) :: field(:)

#ifdef USE_MPI
        integer(wi) isize
        isize = size(field)
        call MPI_File_open(self%communicator, trim(adjustl(name)), &
                           ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, mpio_fh, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_set_view(mpio_fh, self%offset, MPI_REAL4, self%subarray, 'native', MPI_INFO_NULL, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_write_all(mpio_fh, field, isize, MPI_REAL4, status, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_close(mpio_fh, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
#else
#include "tlab_open_file.h"
        write (LOC_UNIT_ID, POS=self%offset + 1) field(:)
        close (LOC_UNIT_ID)
#endif

        return
    end subroutine

!########################################################################
!########################################################################
    subroutine io_subarray_write_double(self, field)
        class(io_subarray_dt) self
        real(dp), intent(in) :: field(:)

#ifdef USE_MPI
        integer(wi) isize
        isize = size(field)
        call MPI_File_open(self%communicator, trim(adjustl(name)), &
                           ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, mpio_fh, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_set_view(mpio_fh, self%offset, MPI_REAL8, self%subarray, 'native', MPI_INFO_NULL, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_write_all(mpio_fh, field, isize, MPI_REAL8, status, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_close(mpio_fh, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
#else
#include "tlab_open_file.h"
        write (LOC_UNIT_ID, POS=self%offset + 1) field(:)
        close (LOC_UNIT_ID)
#endif

        return
    end subroutine

!########################################################################
!########################################################################
    subroutine io_subarray_write_int1(self, field)
        class(io_subarray_dt) self
        integer(1), intent(in) :: field(:)

#ifdef USE_MPI
        integer(wi) isize
        isize = size(field)
        call MPI_File_open(self%communicator, trim(adjustl(name)), &
                           ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), MPI_INFO_NULL, mpio_fh, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_set_view(mpio_fh, self%offset, MPI_INTEGER1, self%subarray, 'native', MPI_INFO_NULL, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_write_all(mpio_fh, field, isize, MPI_INTEGER1, status, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
        call MPI_File_close(mpio_fh, ims_err)
        ! if (ims_err /= MPI_SUCCESS) call TLabMPI_Panic(__FILE__, ims_err)
#else
#include "tlab_open_file.h"
        write (LOC_UNIT_ID, POS=self%offset + 1) field(:)
        close (LOC_UNIT_ID)
#endif

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine Write_Header(nx, ny, nz, nt, params)
        integer(wi), intent(in) :: nx, ny, nz, nt
        real(wp), intent(in), optional :: params(:)

        ! -------------------------------------------------------------------
        integer(wi) offset
        integer nx_total, ny_total, nz_total

        !########################################################################
#ifdef USE_MPI
        if (ims_pro == 0) then
            nx_total = nx*ims_npro_i
            ny_total = ny*ims_npro_j
            nz_total = nz
#else
            nx_total = nx
            ny_total = ny
            nz_total = nz
#endif
#include "tlab_open_file.h"

            offset = 5*SIZEOFINT
            if (present(params)) then
                offset = offset + size(params)*SIZEOFREAL
            end if

            write (LOC_UNIT_ID) offset, nx_total, ny_total, nz_total, nt

            if (present(params)) then
                write (LOC_UNIT_ID) params(:)
            end if

            close (LOC_UNIT_ID)
#ifdef USE_MPI
        end if
        call MPI_BARRIER(MPI_COMM_WORLD, ims_err)
#endif
        return
    end subroutine Write_Header

#undef LOC_UNIT_ID
#undef LOC_STATUS

end module IO_Fields
