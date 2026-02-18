#include "tlab_error.h"

module Planes
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: efile, lfile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Memory, only: imax, jmax, kmax
    use TLab_Memory, only: inb_flow_array, inb_scal_array
    use TLab_Arrays, only: txc
    use TLab_Pointers_3D, only: pointers3d_dt
    use IO_Fields, only: io_subarray_dt
    implicit none
    private

    public :: planesX, planesY, planesZ
    public :: Planes_Save

    ! -------------------------------------------------------------------
    integer, parameter :: TYPE_NONE = 0
    integer, parameter :: TYPE_DEFAULT = 1        ! flow and scal arrays, and pressure.
    integer, parameter :: TYPE_DERIVATIVES = 2    ! Add log of enstrophy and scalar gradients

    type :: planes_dt
        integer :: type = TYPE_NONE
        integer :: size
        integer, allocatable :: nodes(:)
        real(wp), pointer :: data(:, :, :) => null()
        integer io(5)
        type(io_subarray_dt) io_subarray
    end type

    type, extends(planes_dt) :: planesX_dt
    contains
        procedure :: initialize => planesX_initialize_dt
        procedure :: save => planesX_save_dt
    end type

    type, extends(planes_dt) :: planesY_dt
    contains
        procedure :: initialize => planesY_initialize_dt
        procedure :: save => planesY_save_dt
    end type

    type, extends(planes_dt) :: planesZ_dt
    contains
        procedure :: initialize => planesZ_initialize_dt
        procedure :: save => planesZ_save_dt
    end type

    type(planesX_dt) :: planesX
    type(planesY_dt) :: planesY
    type(planesZ_dt) :: planesZ

contains
    ! ###################################################################
    ! ###################################################################
    subroutine planesX_initialize_dt(self, inifile)
        use TLab_Grid, only: x
        use IO_Fields, only: IO_TYPE_SINGLE
#ifdef USE_MPI
        use mpi_f08, only: MPI_REAL4
        use TLabMPI_VARS, only: ims_comm_y, ims_pro_i, ims_npro_i
        use IO_Fields, only: IO_Create_Subarray_YOZ, IO_TYPE_SINGLE
#endif
        class(planesX_dt), intent(out) :: self

        character(len=*), intent(in) :: inifile

        ! ###################################################################
        call planes_initialize(self, inifile, x)

        if (self%type == TYPE_NONE) return

        self%data(1:jmax, 1:self%size, 1:kmax) => txc(1:jmax*self%size*kmax, 1)
        self%io = [size(self%data), 1, size(self%data), 1, 1]

        self%io_subarray%offset = 0
        self%io_subarray%precision = IO_TYPE_SINGLE
#ifdef USE_MPI
        self%io_subarray%active = .false.  ! defaults
        if (ims_npro_i > 1) then
            if (ims_pro_i == (self%nodes(1)/imax)) self%io_subarray%active = .true.
            self%io_subarray%communicator = ims_comm_y
            self%io_subarray%subarray = IO_Create_Subarray_YOZ(jmax, self%size*kmax, MPI_REAL4)
        end if
#endif

        return
    end subroutine

    ! ###################################################################
    ! ###################################################################
    subroutine planesY_initialize_dt(self, inifile)
        use TLab_Grid, only: y
        use IO_Fields, only: IO_TYPE_SINGLE
#ifdef USE_MPI
        use mpi_f08, only: MPI_REAL4
        use TLabMPI_VARS, only: ims_comm_x, ims_pro_j, ims_npro_j
        use IO_Fields, only: IO_Create_Subarray_XOZ, IO_TYPE_SINGLE
#endif
        class(planesY_dt), intent(out) :: self

        character(len=*), intent(in) :: inifile

        ! ###################################################################
        call planes_initialize(self, inifile, y)

        if (self%type == TYPE_NONE) return

        self%data(1:imax, 1:self%size, 1:kmax) => txc(1:imax*self%size*kmax, 1)
        self%io = [size(self%data), 1, size(self%data), 1, 1]

        self%io_subarray%offset = 0
        self%io_subarray%precision = IO_TYPE_SINGLE
#ifdef USE_MPI
        self%io_subarray%active = .false.  ! defaults
        if (ims_npro_j > 1) then
            if (ims_pro_j == (self%nodes(1)/jmax)) self%io_subarray%active = .true.
            self%io_subarray%communicator = ims_comm_x
            self%io_subarray%subarray = IO_Create_Subarray_XOZ(imax, kmax*self%size, MPI_REAL4)
        end if
#endif

        return
    end subroutine

    ! ###################################################################
    ! ###################################################################
    subroutine planesZ_initialize_dt(self, inifile)
        use TLab_Grid, only: z
        use IO_Fields, only: IO_TYPE_SINGLE
#ifdef USE_MPI
        use mpi_f08, only: MPI_REAL4, MPI_COMM_WORLD
        use IO_Fields, only: IO_Create_Subarray_XOY
#endif
        class(planesZ_dt), intent(out) :: self

        character(len=*), intent(in) :: inifile

        ! ###################################################################
        call planes_initialize(self, inifile, z)

        if (self%type == TYPE_NONE) return

        self%data(1:imax, 1:jmax, 1:self%size) => txc(1:imax*jmax*self%size, 1)
        self%io = [size(self%data), 1, size(self%data), 1, 1]

        self%io_subarray%offset = 0
        self%io_subarray%precision = IO_TYPE_SINGLE
#ifdef USE_MPI
        self%io_subarray%active = .true.
        self%io_subarray%communicator = MPI_COMM_WORLD
        self%io_subarray%subarray = IO_Create_Subarray_XOY(imax, jmax, self%size, MPI_REAL4)
#endif

        return
    end subroutine

    ! ###################################################################
    ! ###################################################################
    subroutine planes_initialize(self, inifile, axis)
        use TLab_Grid, only: grid_dt
        class(planes_dt), intent(out) :: self
        type(grid_dt), intent(in) :: axis

        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block, tag
        character(len=128) eStr
        character(len=512) sRes

        integer, parameter :: nodes_size_max = 16
        integer :: nodes(nodes_size_max), nodes_size

        ! ###################################################################
        ! Read
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Planes'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        select case (trim(adjustl(axis%name)))
        case ('x')
            tag = 'PlanesI'
        case ('y')
            tag = 'PlanesJ'
        case ('z')
            tag = 'PlanesK'
        end select

        call ScanFile_Char(bakfile, inifile, block, tag, 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            nodes_size = nodes_size_max; call LIST_INTEGER(sRes, nodes_size, nodes)

            call ScanFile_Char(bakfile, inifile, block, trim(adjustl(tag))//'Type', 'default', sRes)
            if (trim(adjustl(sRes)) == 'none') then; self%type = TYPE_NONE
            elseif (trim(adjustl(sRes)) == 'default') then; self%type = TYPE_DEFAULT
            elseif (trim(adjustl(sRes)) == 'derivatives') then; self%type = TYPE_DERIVATIVES
            else
                call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong Type option.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

        end if

        if (self%type == TYPE_NONE) return

        !########################################################################
        ! Initialize data
        allocate (self%nodes, source=nodes(1:nodes_size))

        if (any(self%nodes(:) < 1) .or. any(self%nodes(:) > axis%size)) then
            call TLab_Write_ASCII(efile, __FILE__//'. Plane nodes out of bounds '//trim(adjustl(axis%name)))
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        self%size = (inb_flow_array + inb_scal_array + 1)*size(self%nodes)
        if (self%type == TYPE_DERIVATIVES) then
            self%size = self%size + (inb_scal_array + 1)*size(self%nodes)
        end if

        if (self%size > axis%size) then
            call TLab_Write_ASCII(efile, __FILE__//'. Array size is insufficient.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        return
    end subroutine planes_initialize

    ! ###################################################################
    ! ###################################################################
    subroutine Planes_Save()
        use TLab_Constants, only: small_wp
        use TLab_Constants, only: fmt_r
        use TLab_Arrays, only: q, s, txc
        use NSE_Pressure
        use TLab_Time
        use FI_VORTICITY_EQN
        use FI_GRADIENT_EQN

        ! -------------------------------------------------------------------
        character*32 str
        character*250 line1

        integer, parameter :: NMAX_VARS = 16
        integer iv, nvars
        type(pointers3d_dt) :: vars(NMAX_VARS)

        ! ###################################################################
        if (all([planesX%type, planesY%type, planesZ%type] == TYPE_NONE)) return

        write (line1, fmt_r) rtime
        write (str, *) itime; line1 = 'at It'//trim(adjustl(str))//' and time '//trim(adjustl(line1))//'.'

        ! ###################################################################
        ! define pointers to data
        ! general order of variables is [u, v, w, {scal1,...}, p, {log(entstrophy), log(grad(scal1)),...)}]
        nvars = 0

        do iv = 1, inb_flow_array
            nvars = nvars + 1; vars(nvars)%field(1:imax, 1:jmax, 1:kmax) => q(1:imax*jmax*kmax, iv)
        end do

        do iv = 1, inb_scal_array
            nvars = nvars + 1; vars(nvars)%field(1:imax, 1:jmax, 1:kmax) => s(1:imax*jmax*kmax, iv)
        end do

        call NSE_Pressure_Incompressible(q, s, txc(:, 2), txc(:, 3), txc(:, 4), txc(:, 5))
        nvars = nvars + 1; vars(nvars)%field(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 2)

        if (any([planesX%type, planesY%type, planesZ%type] == TYPE_DERIVATIVES)) then
            call FI_VORTICITY(imax, jmax, kmax, q(:, 1), q(:, 2), q(:, 3), txc(:, 3), txc(:, 4), txc(:, 5))
            txc(1:imax*jmax*kmax, 3) = log10(txc(1:imax*jmax*kmax, 3) + small_wp)
            nvars = nvars + 1; vars(nvars)%field(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 3)

            do iv = 1, inb_scal_array
                call FI_GRADIENT(imax, jmax, kmax, s(:, iv), txc(:, iv + 3), txc(:, iv + 4))
                txc(1:imax*jmax*kmax, iv + 3) = log10(txc(1:imax*jmax*kmax, iv + 3) + small_wp)
                nvars = nvars + 1; vars(nvars)%field(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, iv + 3)
            end do

        end if

        ! ###################################################################
        write (str, *) itime

        call TLab_Write_ASCII(lfile, 'Writing I-planes '//trim(adjustl(line1)))
        call planesX%save(vars(1:nvars), fname='planesI.'//trim(adjustl(str)))

        call TLab_Write_ASCII(lfile, 'Writing J-planes '//trim(adjustl(line1)))
        call planesY%save(vars(1:nvars), fname='planesJ.'//trim(adjustl(str)))

        call TLab_Write_ASCII(lfile, 'Writing K-planes '//trim(adjustl(line1)))
        call planesZ%save(vars(1:nvars), fname='planesK.'//trim(adjustl(str)))

        return
    end subroutine

    ! ###################################################################
    ! ###################################################################
    subroutine planesX_save_dt(self, vars, fname)
        use IO_Fields, only: IO_Write_Subarray
        class(planesX_dt) self
        type(pointers3d_dt), intent(in) :: vars(:)
        character(len=*), intent(in) :: fname

        integer offset, iv
        integer j, k

        ! ###################################################################
        if (self%type == TYPE_NONE) return

        offset = 0
        do iv = 1, size(vars)
            do k = 1, kmax
                do j = 1, jmax
                    self%data(j, 1 + offset:size(self%nodes) + offset, k) = vars(iv)%field(self%nodes(:), j, k)
                end do
            end do
            offset = offset + size(self%nodes)
        end do

        call IO_Write_Subarray(self%io_subarray, fname, [' '], self%data, self%io)

        return
    end subroutine planesX_save_dt

    ! ###################################################################
    ! ###################################################################
    subroutine planesY_save_dt(self, vars, fname)
        use IO_Fields, only: IO_Write_Subarray
        class(planesY_dt) self
        type(pointers3d_dt), intent(in) :: vars(:)
        character(len=*), intent(in) :: fname

        integer offset, iv

        ! ###################################################################
        if (self%type == TYPE_NONE) return

        offset = 0
        do iv = 1, size(vars)
            self%data(:, 1 + offset:size(self%nodes) + offset, :) = vars(iv)%field(:, self%nodes(:), :)
            offset = offset + size(self%nodes)
        end do

        call IO_Write_Subarray(self%io_subarray, fname, [' '], self%data, self%io)

        return
    end subroutine planesY_save_dt

    ! ###################################################################
    ! ###################################################################
    subroutine planesZ_save_dt(self, vars, fname)
        use IO_Fields, only: IO_Write_Subarray
        class(planesZ_dt) self
        type(pointers3d_dt), intent(in) :: vars(:)
        character(len=*), intent(in) :: fname

        integer offset, iv

        ! ###################################################################
        if (self%type == TYPE_NONE) return

        offset = 0
        do iv = 1, size(vars)
            self%data(:, :, 1 + offset:size(self%nodes) + offset) = vars(iv)%field(:, :, self%nodes(:))
            offset = offset + size(self%nodes)
        end do

        call IO_Write_Subarray(self%io_subarray, fname, [' '], self%data, self%io)

        return
    end subroutine planesZ_save_dt

end module Planes
