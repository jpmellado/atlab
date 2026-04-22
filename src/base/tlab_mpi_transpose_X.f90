#include "tlab_error.h"

! Circular transposition within directional communicators
module TLabMPI_Transpose_X
    use mpi_f08
    use TLab_Constants, only: wp, dp, sp, wi, sizeofreal
    use TLab_Constants, only: lfile, efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLabMPI_VARS, only: ims_err
    implicit none
    private

    public :: TLabMPI_Trp_Initialize_X
    public :: tmpi_transpose_x_dt
    public :: tmpi_transpose_y_dt
    public :: tmpi_trp_X, tmpi_trp_Y    ! general plans used in derivatives and other operators
    !                                   I wonder if here or somewhere else

    ! -----------------------------------------------------------------------
    type trp_mem_dt
        type(MPI_Datatype) :: type              ! derived types
        integer(wi), allocatable :: disp(:)     ! buffer displacements
        integer(wi), allocatable :: map(:)      ! processor mapping
    end type

    type :: tmpi_transpose_dt
        ! sequence
        integer :: mode                         ! asynchronous, sendrecv, alltoall
        type(trp_mem_dt) :: send                ! send information
        type(trp_mem_dt) :: recv                ! recv information
        type(MPI_Comm) :: comm                  ! communicator
        integer :: size_block_processes
        integer(wi) :: nlines
        integer(wi) :: size3d
    contains
        private
        procedure :: tmpi_trp_forward_real
        procedure :: tmpi_trp_forward_complex
        procedure :: tmpi_trp_backward_real
        procedure :: tmpi_trp_backward_complex
        generic, public :: forward => tmpi_trp_forward_complex, tmpi_trp_forward_real
        generic, public :: backward => tmpi_trp_backward_complex, tmpi_trp_backward_real
    end type tmpi_transpose_dt

    type, extends(tmpi_transpose_dt) :: tmpi_transpose_x_dt
    contains
        procedure :: initialize => tmpi_trp_initialize_x
    end type

    type, extends(tmpi_transpose_dt) :: tmpi_transpose_y_dt
    contains
        procedure :: initialize => tmpi_trp_initialize_y
    end type

    type(tmpi_transpose_x_dt) :: tmpi_trp_X
    type(tmpi_transpose_y_dt) :: tmpi_trp_Y

    ! -----------------------------------------------------------------------
    integer :: trp_mode_i, trp_mode_j                               ! Mode of transposition
    integer, parameter :: TLAB_MPI_TRP_NONE = 0
    integer, parameter :: TLAB_MPI_TRP_ASYNCHRONOUS = 1
    integer, parameter :: TLAB_MPI_TRP_SENDRECV = 2
    integer, parameter :: TLAB_MPI_TRP_ALLTOALL = 3

    type(MPI_Datatype) :: trp_datatype_i, trp_datatype_j            ! Transposition in double or single precision
    real(wp), allocatable, target :: wrk_mpi(:)                     ! 3D work array for datatype conversion
    real(sp), pointer :: a_wrk(:) => null(), b_wrk(:) => null()

    ! -----------------------------------------------------------------------
    integer(wi) :: trp_sizBlock_i, trp_sizBlock_j                   ! explicit sed/recv: group sizes of rend/recv messages

    type(MPI_Datatype), allocatable :: types_send(:), types_recv(:) ! alltoallw
    integer, allocatable :: counts(:)

    ! -----------------------------------------------------------------------
    type(MPI_Status), allocatable :: status(:)
    type(MPI_Request), allocatable :: request(:)
    integer ims_tag

contains
    ! ######################################################################
    ! ######################################################################
    subroutine TLabMPI_Trp_Initialize_X(inifile)
        use TLabMPI_VARS, only: xMpi, yMpi
        use TLab_Memory, only: isize_wrk3d, imax, jmax, kmax
        use TLab_Memory, only: TLab_Allocate_Real
        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        integer(wi) ip

        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes, line

        ! #######################################################################
        ! Read data
        bakfile = trim(adjustl(inifile))//'.bak'

        block = 'Parallel'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call ScanFile_Char(bakfile, inifile, block, 'TransposeModeI', 'asynchronous', sRes)
        if (trim(adjustl(sRes)) == 'none') then; trp_mode_i = TLAB_MPI_TRP_NONE
        elseif (trim(adjustl(sRes)) == 'asynchronous') then; trp_mode_i = TLAB_MPI_TRP_ASYNCHRONOUS
        elseif (trim(adjustl(sRes)) == 'sendrecv') then; trp_mode_i = TLAB_MPI_TRP_SENDRECV
        elseif (trim(adjustl(sRes)) == 'alltoall') then; trp_mode_i = TLAB_MPI_TRP_ALLTOALL
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong TransposeModeI option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'TransposeModeJ', 'asynchronous', sRes)
        if (trim(adjustl(sRes)) == 'none') then; trp_mode_j = TLAB_MPI_TRP_NONE
        elseif (trim(adjustl(sRes)) == 'asynchronous') then; trp_mode_j = TLAB_MPI_TRP_ASYNCHRONOUS
        elseif (trim(adjustl(sRes)) == 'sendrecv') then; trp_mode_j = TLAB_MPI_TRP_SENDRECV
        elseif (trim(adjustl(sRes)) == 'alltoall') then; trp_mode_j = TLAB_MPI_TRP_ALLTOALL
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong TransposeModeJ option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'TransposeTypeI', 'Double', sRes)
        if (trim(adjustl(sRes)) == 'double') then; trp_datatype_i = MPI_REAL8
        elseif (trim(adjustl(sRes)) == 'single') then; trp_datatype_i = MPI_REAL4
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong TransposeTypeI.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'TransposeTypeJ', 'Double', sRes)
        if (trim(adjustl(sRes)) == 'double') then; trp_datatype_j = MPI_REAL8
        elseif (trim(adjustl(sRes)) == 'single') then; trp_datatype_j = MPI_REAL4
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong TransposeTypeJ.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        ! #######################################################################
        ! Initialize

        ! Size of communication in explicit send/recv
#ifdef HLRS_HAWK
        ! On hawk, we tested that 192 yields optimum performance;
        ! Blocking will thus only take effect in very large cases
        trp_sizBlock_j = 192
        trp_sizBlock_i = 384
#else
        ! We assume that this will help to release some of the very heavy
        ! network load in transpositions on most systems
        trp_sizBlock_j = 64
        trp_sizBlock_i = 128
        ! trp_sizBlock_j=1e5   -- would essentially switch off the blocking
#endif

        if (xMpi%num_processors > trp_sizBlock_i) then
            write (line, *) trp_sizBlock_i
            line = 'Using blocking of '//trim(adjustl(line))//' in TLabMPI_TRP<F,B>_I'
            call TLab_Write_ASCII(lfile, line)
        end if

        if (yMpi%num_processors > trp_sizBlock_j) then
            write (line, *) trp_sizBlock_j
            line = 'Using blocking of '//trim(adjustl(line))//' in TLabMPI_TRP<F,B>_K'
            call TLab_Write_ASCII(lfile, line)
        end if

        allocate (status(2*max(trp_sizBlock_i, trp_sizBlock_j, xMpi%num_processors, yMpi%num_processors)))
        allocate (request(2*max(trp_sizBlock_i, trp_sizBlock_j, xMpi%num_processors, yMpi%num_processors)))

        ! -----------------------------------------------------------------------
        ! to use single transposition when running in double precision
        ! call TLab_Allocate_Real(__FILE__, wrk_mpi, [isize_wrk3d], 'wrk-mpi')
        ! isize_wrk3d is not yet defined; see if you need to move this somewhere else
        if (any([trp_datatype_j, trp_datatype_j] == MPI_REAL4)) then
            call TLab_Allocate_Real(__FILE__, wrk_mpi, [imax*jmax*kmax], 'wrk-mpi')
        end if

        ! -----------------------------------------------------------------------
        ! to use alltoallw
        allocate (counts(max(xMpi%num_processors, yMpi%num_processors)))!, zMpi%num_processors)))
        allocate (types_send(max(xMpi%num_processors, yMpi%num_processors)))!, zMpi%num_processors)))
        allocate (types_recv(max(xMpi%num_processors, yMpi%num_processors)))!, zMpi%num_processors)))
        counts(:) = 1

        ! -----------------------------------------------------------------------
        ! Create basic transposition plans used for partial X and partial Z; could be in another module...
        if (xMpi%num_processors > 1) then
            call tmpi_trp_X%initialize(imax, jmax*kmax, message='Ox derivatives.')
        end if

        if (yMpi%num_processors > 1) then
            call tmpi_trp_Y%initialize(jmax, imax*kmax, message='Oy derivatives.')
        end if

        return
    end subroutine

    ! ######################################################################
    ! ######################################################################
    subroutine tmpi_trp_initialize_x(self, nmax, npage, locStride, locType, message)
        use TLabMPI_VARS, only: xMpi
        class(tmpi_transpose_x_dt), intent(out) :: self
        integer(wi), intent(in) :: npage, nmax
        integer(wi), intent(in), optional :: locStride
        type(MPI_Datatype), intent(in), optional :: locType
        character(len=*), intent(in), optional :: message

        ! -----------------------------------------------------------------------
        integer(wi) i
        type(MPI_Datatype) :: datatype
        integer block_count, block_length, stride
        integer ims_ss, ims_rs
        character*64 str, line

        ! #######################################################################
        self%mode = trp_mode_i
        self%comm = xMpi%comm

        if (present(message)) &
            call TLab_Write_ASCII(lfile, 'Creating derived MPI types for '//trim(adjustl(message)))

        if (mod(npage, xMpi%num_processors) == 0) then
            self%nlines = npage/xMpi%num_processors
            allocate (self%send%disp(xMpi%num_processors), self%recv%disp(xMpi%num_processors))
            self%size3d = npage*nmax
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Ratio npage/npro not an integer.')
            call TLab_Stop(DNS_ERROR_PARPARTITION)
        end if

        block_count = self%nlines
        block_length = nmax

        ! Calculate array displacements in Forward Send/Receive
        self%send%disp(1) = 0
        self%recv%disp(1) = 0
        do i = 2, xMpi%num_processors
            self%send%disp(i) = self%send%disp(i - 1) + block_length*block_count
            self%recv%disp(i) = self%recv%disp(i - 1) + block_length
        end do

        ! #######################################################################
        if (present(locType)) then
            datatype = locType
        else
            datatype = trp_datatype_i
        end if

        stride = block_length                   ! stride = block_length because things are together
        call MPI_TYPE_VECTOR(block_count, block_length, stride, datatype, self%send%type, ims_err)
        call MPI_TYPE_COMMIT(self%send%type, ims_err)

        stride = nmax*xMpi%num_processors       ! stride is a multiple of nmax_total=nmax*xMpi%num_processors
        call MPI_TYPE_VECTOR(block_count, block_length, stride, datatype, self%recv%type, ims_err)
        call MPI_TYPE_COMMIT(self%recv%type, ims_err)

        ! -----------------------------------------------------------------------
        call MPI_TYPE_SIZE(self%send%type, ims_ss, ims_err)
        call MPI_TYPE_SIZE(self%recv%type, ims_rs, ims_err)

        if (ims_ss /= ims_rs) then
            write (str, *) ims_ss; write (line, *) ims_rs
            line = 'Send size '//trim(adjustl(str))//'differs from recv size '//trim(adjustl(line))
            call TLab_Write_ASCII(efile, line)
            call TLab_Stop(DNS_ERROR_MPITYPECHECK)
        end if

        ! -----------------------------------------------------------------------
        self%size_block_processes = trp_sizBlock_i

        ! -----------------------------------------------------------------------
        ! local PE mappings for explicit send/recv
        call explicit_mapping(self%send, self%recv, xMpi)

        return
    end subroutine tmpi_trp_initialize_x

    ! ######################################################################
    ! ######################################################################
    subroutine tmpi_trp_initialize_y(self, nmax, npage, locStride, locType, message)
        use TLabMPI_VARS, only: yMpi
        class(tmpi_transpose_y_dt), intent(out) :: self
        integer(wi), intent(in) :: npage, nmax
        integer(wi), intent(in), optional :: locStride
        type(MPI_Datatype), intent(in), optional :: locType
        character(len=*), intent(in), optional :: message

        ! -----------------------------------------------------------------------
        integer(wi) i
        type(MPI_Datatype) :: datatype
        integer block_count, block_length, stride
        integer ims_ss, ims_rs
        character*64 str, line

        ! #######################################################################
        self%mode = trp_mode_j
        self%comm = yMpi%comm

        if (present(message)) &
            call TLab_Write_ASCII(lfile, 'Creating derived MPI types for '//trim(adjustl(message)))

        if (mod(npage, yMpi%num_processors) == 0) then
            self%nlines = npage/yMpi%num_processors
            allocate (self%send%disp(yMpi%num_processors), self%recv%disp(yMpi%num_processors))
            self%size3d = npage*nmax
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Ratio npage/npro not an integer.')
            call TLab_Stop(DNS_ERROR_PARPARTITION)
        end if

        block_count = nmax
        block_length = self%nlines

        ! Calculate array displacements in Forward Send/Receive
        self%send%disp(1) = 0
        self%recv%disp(1) = 0
        do i = 2, yMpi%num_processors
            self%send%disp(i) = self%send%disp(i - 1) + block_length
            self%recv%disp(i) = self%recv%disp(i - 1) + block_length*block_count
        end do

        ! #######################################################################
        if (present(locType)) then
            datatype = locType
        else
            datatype = trp_datatype_i
        end if

        stride = npage
        call MPI_TYPE_VECTOR(block_count, block_length, stride, datatype, self%send%type, ims_err)
        call MPI_TYPE_COMMIT(self%send%type, ims_err)

        stride = block_length       ! stride = block_length to put things together
        call MPI_TYPE_VECTOR(block_count, block_length, stride, datatype, self%recv%type, ims_err)
        call MPI_TYPE_COMMIT(self%recv%type, ims_err)

        ! -----------------------------------------------------------------------
        call MPI_TYPE_SIZE(self%send%type, ims_ss, ims_err)
        call MPI_TYPE_SIZE(self%recv%type, ims_rs, ims_err)

        if (ims_ss /= ims_rs) then
            write (str, *) ims_ss; write (line, *) ims_rs
            line = 'Send size '//trim(adjustl(str))//'differs from recv size '//trim(adjustl(line))
            call TLab_Write_ASCII(efile, line)
            call TLab_Stop(DNS_ERROR_MPITYPECHECK)
        end if

        ! -----------------------------------------------------------------------
        self%size_block_processes = trp_sizBlock_j

        ! -----------------------------------------------------------------------
        ! local PE mappings for explicit send/recv
        call explicit_mapping(self%send, self%recv, yMpi)

        return
    end subroutine tmpi_trp_initialize_y

    subroutine explicit_mapping(send, recv, axis)
        use TLabMPI_VARS, only: mpi_axis_dt
        type(trp_mem_dt), intent(inout) :: send                ! send information
        type(trp_mem_dt), intent(inout) :: recv                ! recv information
        type(mpi_axis_dt), intent(in) :: axis

        integer ip

        allocate (send%map(axis%num_processors))
        allocate (recv%map(axis%num_processors))
        do ip = 0, axis%num_processors - 1
            send%map(ip + 1) = ip
            recv%map(ip + 1) = mod(axis%num_processors - ip, axis%num_processors)
        end do
        send%map = cshift(send%map, axis%rank)
        recv%map = cshift(recv%map, -axis%rank)

        return
    end subroutine

    ! ######################################################################
    ! ######################################################################
    subroutine tmpi_trp_forward_real(self, a, b)
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
        class(tmpi_transpose_dt), intent(in) :: self
        real(wp), intent(in) :: a(:)
        real(wp), intent(out) :: b(:)

        target b

        ! -----------------------------------------------------------------------
        integer(wi) size

#ifdef PROFILE_ON
        real(wp) time_loc_1, time_loc_2
#endif

        ! #######################################################################
#ifdef PROFILE_ON
        time_loc_1 = MPI_WTIME()
#endif
        if (trp_datatype_i == MPI_REAL4 .and. wp == dp) then
            size = self%size3d
            call c_f_pointer(c_loc(b), a_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi), b_wrk, shape=[size])
            a_wrk(1:size) = real(a(1:size), sp)
            call tmpi_trp_single(a_wrk, self%send, b_wrk, self%recv, &
                                 self%comm, self%size_block_processes, self%mode)
            b(1:size) = real(b_wrk(1:size), dp)
            nullify (a_wrk, b_wrk)
        else
            call tmpi_trp_double(a, self%send, b, self%recv, &
                                 self%comm, self%size_block_processes, self%mode)
        end if

#ifdef PROFILE_ON
        time_loc_2 = MPI_WTIME()
        ims_time_trans = ims_time_trans + (time_loc_2 - time_loc_1)
#endif

        return
    end subroutine

    ! ######################################################################
    ! ######################################################################
    subroutine tmpi_trp_backward_real(self, a, b)
        use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
        class(tmpi_transpose_dt), intent(in) :: self
        real(wp), intent(in) :: a(:)
        real(wp), intent(out) :: b(:)

        target b

        ! -----------------------------------------------------------------------
        integer(wi) size

#ifdef PROFILE_ON
        real(wp) time_loc_1, time_loc_2
#endif

        ! #######################################################################
#ifdef PROFILE_ON
        time_loc_1 = MPI_WTIME()
#endif
        if (trp_datatype_i == MPI_REAL4 .and. wp == dp) then
            size = self%size3d
            call c_f_pointer(c_loc(b), a_wrk, shape=[size])
            call c_f_pointer(c_loc(wrk_mpi), b_wrk, shape=[size])
            a_wrk(1:size) = real(a(1:size), sp)
            call tmpi_trp_single(a_wrk, self%recv, b_wrk, self%send, &
                                 self%comm, self%size_block_processes, self%mode)
            b(1:size) = real(b_wrk(1:size), dp)
            nullify (a_wrk, b_wrk)
        else
            call tmpi_trp_double(a, self%recv, b, self%send, &
                                 self%comm, self%size_block_processes, self%mode)
        end if

#ifdef PROFILE_ON
        time_loc_2 = MPI_WTIME()
        ims_time_trans = ims_time_trans + (time_loc_2 - time_loc_1)
#endif

        return
    end subroutine

    ! ######################################################################
    ! ######################################################################
    subroutine tmpi_trp_forward_complex(self, a, b)
        class(tmpi_transpose_dt), intent(in) :: self
        complex(wp), intent(in) :: a(:)
        complex(wp), intent(out) :: b(:)

        call tmpi_trp_complex(a, self%send, b, self%recv, &
                              self%comm, self%size_block_processes, self%mode)

        return
    end subroutine

    ! ######################################################################
    ! ######################################################################
    subroutine tmpi_trp_backward_complex(self, a, b)
        class(tmpi_transpose_dt), intent(in) :: self
        complex(wp), intent(in) :: a(:)
        complex(wp), intent(out) :: b(:)

        call tmpi_trp_complex(a, self%recv, b, self%send, &
                              self%comm, self%size_block_processes, self%mode)

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine tmpi_trp_double(in, send, out, recv, &
                               comm, step, mode)
        real(wp), intent(in) :: in(*)
        real(wp), intent(out) :: out(*)
        type(trp_mem_dt), intent(in) :: send, recv
        type(MPI_Comm), intent(in) :: comm
        integer(wi), intent(in) :: step
        integer, intent(in) :: mode

        ! -----------------------------------------------------------------------
        integer npro
        integer(wi) j, l, m, ns, nr, ips, ipr

        ! #######################################################################
        npro = size(send%disp(:))

        select case (mode)
        case (TLAB_MPI_TRP_ASYNCHRONOUS)
            do j = 1, npro, step
                l = 0
                do m = j, min(j + step - 1, npro)
                    ns = send%map(m) + 1; ips = ns - 1
                    nr = recv%map(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_ISEND(in(send%disp(ns) + 1), 1, send%type, ips, ims_tag, comm, request(l), ims_err)
                    l = l + 1
                    call MPI_IRECV(out(recv%disp(nr) + 1), 1, recv%type, ipr, ims_tag, comm, request(l), ims_err)
                end do
                call MPI_WAITALL(l, request, status, ims_err)
            end do

        case (TLAB_MPI_TRP_SENDRECV)
            do j = 1, npro, step
                do m = j, min(j + step - 1, npro)
                    ns = send%map(m) + 1; ips = ns - 1
                    nr = recv%map(m) + 1; ipr = nr - 1
                    call MPI_SENDRECV(in(send%disp(ns) + 1), 1, send%type, ips, ims_tag, &
                                      out(recv%disp(nr) + 1), 1, recv%type, ipr, ims_tag, comm, status(1), ims_err)
                end do
            end do

        case (TLAB_MPI_TRP_ALLTOALL)
            types_send(1:npro) = send%type
            types_recv(1:npro) = recv%type
            call MPI_ALLTOALLW(in, counts, send%disp*int(sizeof(1.0_wp)), types_send, &
                               out, counts, recv%disp*int(sizeof(1.0_wp)), types_recv, comm, ims_err)
            ! call MPI_ALLTOALLW(in, spread(1, 1, npro), send%disp*int(sizeof(1.0_wp)), spread(send%type, 1, npro), &
            !                    out, spread(1, 1, npro), recv%disp*int(sizeof(1.0_wp)), spread(recv%type, 1, npro), comm, ims_err)

        end select

        return
    end subroutine tmpi_trp_double

    !########################################################################
    !########################################################################
    subroutine tmpi_trp_single(in, send, out, recv, &
                               comm, step, mode)
        real(sp), intent(in) :: in(*)
        real(sp), intent(out) :: out(*)
        type(trp_mem_dt), intent(in) :: send, recv
        type(MPI_Comm), intent(in) :: comm
        integer(wi), intent(in) :: step
        integer, intent(in) :: mode

        ! -----------------------------------------------------------------------
        integer npro
        integer(wi) j, l, m, ns, nr, ips, ipr

        ! #######################################################################
        npro = size(send%disp(:))

        select case (mode)
        case (TLAB_MPI_TRP_ASYNCHRONOUS)
            do j = 1, npro, step
                l = 0
                do m = j, min(j + step - 1, npro)
                    ns = send%map(m) + 1; ips = ns - 1
                    nr = recv%map(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_ISEND(in(send%disp(ns) + 1), 1, send%type, ips, ims_tag, comm, request(l), ims_err)
                    l = l + 1
                    call MPI_IRECV(out(recv%disp(nr) + 1), 1, recv%type, ipr, ims_tag, comm, request(l), ims_err)
                end do
                call MPI_WAITALL(l, request, status, ims_err)
            end do

        case (TLAB_MPI_TRP_SENDRECV)
            do j = 1, npro, step
                do m = j, min(j + step - 1, npro)
                    ns = send%map(m) + 1; ips = ns - 1
                    nr = recv%map(m) + 1; ipr = nr - 1
                    call MPI_SENDRECV(in(send%disp(ns) + 1), 1, send%type, ips, ims_tag, &
                                      out(recv%disp(nr) + 1), 1, recv%type, ipr, ims_tag, comm, status(1), ims_err)
                end do
            end do

        case (TLAB_MPI_TRP_ALLTOALL)
            types_send(1:npro) = send%type
            types_recv(1:npro) = recv%type
            call MPI_ALLTOALLW(in, counts, send%disp*int(sizeof(1.0_sp)), types_send, &
                               out, counts, recv%disp*int(sizeof(1.0_sp)), types_recv, comm, ims_err)
            ! call MPI_ALLTOALLW(in, spread(1, 1, npro), send%disp*int(sizeof(1.0_wp)), spread(send%type, 1, npro), &
            !                    out, spread(1, 1, npro), recv%disp*int(sizeof(1.0_wp)), spread(recv%type, 1, npro), comm, ims_err)

        end select

        return
    end subroutine tmpi_trp_single

    !########################################################################
    !########################################################################
    subroutine tmpi_trp_complex(in, send, out, recv, &
                                comm, step, mode)
        complex(wp), intent(in) :: in(*)
        complex(wp), intent(out) :: out(*)
        type(trp_mem_dt), intent(in) :: send, recv
        type(MPI_Comm), intent(in) :: comm
        integer(wi), intent(in) :: step
        integer, intent(in) :: mode

        ! -----------------------------------------------------------------------
        integer npro
        integer(wi) j, l, m, ns, nr, ips, ipr

        ! #######################################################################
        npro = size(send%disp(:))

        select case (mode)
        case (TLAB_MPI_TRP_ASYNCHRONOUS)
            do j = 1, npro, step
                l = 0
                do m = j, min(j + step - 1, npro)
                    ns = send%map(m) + 1; ips = ns - 1
                    nr = recv%map(m) + 1; ipr = nr - 1
                    l = l + 1
                    call MPI_ISEND(in(send%disp(ns) + 1), 1, send%type, ips, ims_tag, comm, request(l), ims_err)
                    l = l + 1
                    call MPI_IRECV(out(recv%disp(nr) + 1), 1, recv%type, ipr, ims_tag, comm, request(l), ims_err)
                end do
                call MPI_WAITALL(l, request, status, ims_err)
            end do

        case (TLAB_MPI_TRP_SENDRECV)
            do j = 1, npro, step
                do m = j, min(j + step - 1, npro)
                    ns = send%map(m) + 1; ips = ns - 1
                    nr = recv%map(m) + 1; ipr = nr - 1
                    call MPI_SENDRECV(in(send%disp(ns) + 1), 1, send%type, ips, ims_tag, &
                                      out(recv%disp(nr) + 1), 1, recv%type, ipr, ims_tag, comm, status(1), ims_err)
                end do
            end do

        case (TLAB_MPI_TRP_ALLTOALL)
            types_send(1:npro) = send%type
            types_recv(1:npro) = recv%type
            call MPI_ALLTOALLW(in, counts, send%disp*int(2*sizeof(1.0_wp)), types_send, &
                               out, counts, recv%disp*int(2*sizeof(1.0_wp)), types_recv, comm, ims_err)
            ! call MPI_ALLTOALLW(in, spread(1, 1, npro), send%disp*int(sizeof(1.0_wp)), spread(send%type, 1, npro), &
            !                    out, spread(1, 1, npro), recv%disp*int(sizeof(1.0_wp)), spread(recv%type, 1, npro), comm, ims_err)

        end select

        return
    end subroutine tmpi_trp_complex

end module TLabMPI_Transpose_X
