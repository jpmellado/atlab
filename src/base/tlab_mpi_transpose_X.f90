#include "tlab_error.h"

! Circular transposition within directional communicators
! combined with local transpositions to avoid the need of MPI derived types
module TLabMPI_Transpose_X
    use mpi_f08
    use TLab_Constants, only: wp, dp, wi
    use TLab_Constants, only: lfile, efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLabMPI_VARS, only: TLAB_MPI_REAL_TYPE
    implicit none
    private

    public :: tmpi_transposeX_inner_dt  ! transpose along inner most index
    public :: tmpi_transposeX_outer_dt  ! transpose along outer most index

    ! -----------------------------------------------------------------------
    type trp_mem_dt
        type(MPI_Datatype) :: type              ! derived types
        integer(wi) :: count                    ! number of elements of type type to be transferred
        integer(wi), allocatable :: disp(:)     ! buffer displacements
        integer(wi), allocatable :: map(:)      ! processor mapping
    end type

    type :: tmpi_transposeX_dt
        ! sequence
        type(trp_mem_dt) :: send                ! send information
        type(trp_mem_dt) :: recv                ! recv information
        type(MPI_Comm) :: comm                  ! communicator
        integer :: size_block_processes
        integer(wi) :: nlines
    contains
        procedure :: initialize => tmpi_trpX_initialize
    end type tmpi_transposeX_dt

    ! -----------------------------------------------------------------------
    type, extends(tmpi_transposeX_dt) :: tmpi_transposeX_inner_dt
    contains
        private
        ! procedure :: tmpi_trpX_inner_fwd_real
        procedure :: tmpi_trpX_inner_fwd_real_overlap
        procedure :: tmpi_trpX_inner_bwd_real
        procedure :: tmpi_trpX_inner_fwd_complex
        procedure :: tmpi_trpX_inner_fwd_complex_inplace
        procedure :: tmpi_trpX_inner_bwd_complex
        procedure :: tmpi_trpX_inner_bwd_complex_inplace
        ! generic, public :: forward => tmpi_trpX_inner_fwd_real, tmpi_trpX_inner_fwd_complex, tmpi_trpX_inner_fwd_complex_inplace
        generic, public :: forward => tmpi_trpX_inner_fwd_real_overlap, tmpi_trpX_inner_fwd_complex, tmpi_trpX_inner_fwd_complex_inplace
        generic, public :: backward => tmpi_trpX_inner_bwd_real, tmpi_trpX_inner_bwd_complex, tmpi_trpX_inner_bwd_complex_inplace
    end type

    type, extends(tmpi_transposeX_dt) :: tmpi_transposeX_outer_dt
    contains
        private
        procedure, public :: forward => tmpi_trpX_outer_fwd_complex
        procedure, public :: backward => tmpi_trpX_outer_bwd_complex
    end type

    ! -----------------------------------------------------------------------
    ! Size of communication in explicit send/recv
    ! We assume that this will help to release some of the very heavy
    ! network load in transpositions on most systems
#ifdef HLRS_HAWK
    ! On hawk, we tested that 192 yields optimum performance;
    ! Blocking will thus only take effect in very large cases
    integer(wi) :: tmpi_sizeblock = 384
#else
    integer(wi) :: tmpi_sizeblock = 128
    ! integer(wi) :: trp_sizBlock_j=1e5   -- would essentially switch off the blocking
#endif

    ! -----------------------------------------------------------------------
    type(MPI_Status), allocatable :: status(:)
    type(MPI_Request), allocatable :: request(:)
    integer ims_tag

    integer ims_err

contains
    ! ######################################################################
    ! ######################################################################
    subroutine tmpi_trpX_initialize(self, nmax, npage, axis, locType, message)
        use TLabMPI_VARS, only: mpi_axis_dt
        class(tmpi_transposeX_dt), intent(out) :: self
        integer(wi), intent(in) :: npage, nmax
        type(mpi_axis_dt), intent(in) :: axis
        type(MPI_Datatype), intent(in), optional :: locType
        character(len=*), intent(in), optional :: message

        ! -----------------------------------------------------------------------
        integer(wi) i, locSize

        ! #######################################################################
        self%comm = axis%comm

        if (present(message)) &
            call TLab_Write_ASCII(lfile, 'Creating derived MPI types for '//trim(adjustl(message)))

        if (mod(npage, axis%num_processors) == 0) then
            self%nlines = npage/axis%num_processors
            allocate (self%send%disp(axis%num_processors), self%recv%disp(axis%num_processors))
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Ratio npage/npro not an integer.')
            call TLab_Stop(DNS_ERROR_PARPARTITION)
        end if

        ! Calculate array displacements in Forward Send/Receive
        self%send%disp(1) = 0
        self%recv%disp(1) = 0
        do i = 2, axis%num_processors
            self%send%disp(i) = self%send%disp(i - 1) + self%nlines*nmax
            self%recv%disp(i) = self%recv%disp(i - 1) + self%nlines*nmax
        end do

        ! Define types
        if (present(locType)) then
            self%send%type = locType
            self%recv%type = locType
        else
            self%send%type = TLAB_MPI_REAL_TYPE
            self%recv%type = TLAB_MPI_REAL_TYPE
        end if

        self%send%count = self%nlines*nmax
        self%recv%count = self%nlines*nmax

        ! -----------------------------------------------------------------------
        self%size_block_processes = tmpi_sizeblock

        ! -----------------------------------------------------------------------
        ! local PE mappings for explicit send/recv
        call explicit_mapping(self%send, self%recv, axis)

        ! -----------------------------------------------------------------------
        if (allocated(status)) then
            locSize = size(status)
            deallocate (status)
        else
            locSize = 0
        end if
        locSize = max(locSize, 2*max(tmpi_sizeblock, axis%num_processors))
        allocate (status(locSize))

        if (allocated(request)) then
            locSize = size(request)
            deallocate (request)
        else
            locSize = 0
        end if
        locSize = max(locSize, 2*max(tmpi_sizeblock, axis%num_processors))
        allocate (request(locSize))

        return
    end subroutine tmpi_trpX_initialize

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
    subroutine tmpi_trpX_inner_fwd_complex(self, a, b, wrk)
        use TLab_Transpose
        class(tmpi_transposeX_inner_dt), intent(in) :: self
        complex(wp), intent(in) :: a(:)
        complex(wp), intent(out) :: b(:)
        complex(wp), intent(inout) :: wrk(:)

        integer ib, nblocks
        integer(wi) count, nmax, ip

        nblocks = size(self%send%disp(:)) ! # of blocks = # of processors
        count = self%send%disp(2)   ! # all processors transfer same amount
        nmax = count/self%nlines ! size along direction

        do ib = 1, nblocks
            ip = (ib - 1)*count + 1
            call TLab_Transpose_Complex(a(ip:), nmax, self%nlines, nmax, wrk(ip:), self%nlines)
        end do

        call tmpi_trpX_complex(wrk, self%send, b, self%recv, &
                               self%comm, self%size_block_processes)

        return
    end subroutine

    subroutine tmpi_trpX_inner_fwd_complex_inplace(self, a, wrk)
        use TLab_Transpose
        class(tmpi_transposeX_inner_dt), intent(in) :: self
        complex(wp), intent(inout) :: a(:)
        complex(wp), intent(inout) :: wrk(:)

        integer ib, nblocks
        integer(wi) count, nmax, ip

        nblocks = size(self%send%disp(:)) ! # of blocks = # of processors
        count = self%send%disp(2)   ! # all processors transfer same amount
        nmax = count/self%nlines ! size along direction

        do ib = 1, nblocks
            ip = (ib - 1)*count + 1
            call TLab_Transpose_Complex(a(ip:), nmax, self%nlines, nmax, wrk(ip:), self%nlines)
        end do

        call tmpi_trpX_complex(wrk, self%send, a, self%recv, &
                               self%comm, self%size_block_processes)

        return
    end subroutine

    subroutine tmpi_trpX_inner_bwd_complex(self, a, b, wrk)
        use TLab_Transpose
        class(tmpi_transposeX_inner_dt), intent(in) :: self
        complex(wp), intent(in) :: a(:)
        complex(wp), intent(out) :: b(:)
        complex(wp), intent(inout) :: wrk(:)

        integer ib, nblocks
        integer(wi) count, nmax, ip

        call tmpi_trpX_complex(a, self%recv, wrk, self%send, &
                               self%comm, self%size_block_processes)

        nblocks = size(self%send%disp(:)) ! # of blocks = # of processors
        count = self%send%disp(2)   ! # all processors transfer same amount
        nmax = count/self%nlines ! size along direction

        do ib = 1, nblocks
            ip = (ib - 1)*count + 1
            call TLab_Transpose_Complex(wrk(ip:), self%nlines, nmax, self%nlines, b(ip:), nmax)
        end do

        return
    end subroutine

    subroutine tmpi_trpX_inner_bwd_complex_inplace(self, a, wrk)
        use TLab_Transpose
        class(tmpi_transposeX_inner_dt), intent(in) :: self
        complex(wp), intent(inout) :: a(:)
        complex(wp), intent(inout) :: wrk(:)

        integer ib, nblocks
        integer(wi) count, nmax, ip

        call tmpi_trpX_complex(a, self%recv, wrk, self%send, &
                               self%comm, self%size_block_processes)

        nblocks = size(self%send%disp(:)) ! # of blocks = # of processors
        count = self%send%disp(2)   ! # all processors transfer same amount
        nmax = count/self%nlines ! size along direction

        do ib = 1, nblocks
            ip = (ib - 1)*count + 1
            call TLab_Transpose_Complex(wrk(ip:), self%nlines, nmax, self%nlines, a(ip:), nmax)
        end do

        return
    end subroutine

    ! ######################################################################
    ! ######################################################################
    subroutine tmpi_trpX_inner_fwd_real(self, a, b, wrk)
        use TLab_Transpose
        class(tmpi_transposeX_inner_dt), intent(in) :: self
        real(wp), intent(in) :: a(:)
        real(wp), intent(out) :: b(:)
        real(wp), intent(inout) :: wrk(:)

        integer ib, nblocks
        integer(wi) count, nmax, ip

        ! Make first index last
        nblocks = size(self%send%disp(:)) ! # of blocks = # of processors
        count = self%send%disp(2)   ! # all processors transfer same amount
        nmax = count/self%nlines ! size along direction

        do ib = 1, nblocks
            ip = (ib - 1)*count + 1
            call TLab_Transpose_Real(a(ip:), nmax, self%nlines, nmax, wrk(ip:), self%nlines)
        end do

        call tmpi_trpX_real(wrk, self%send, b, self%recv, &
                            self%comm, self%size_block_processes)

        return
    end subroutine

    subroutine tmpi_trpX_inner_fwd_real_overlap(self, a, b, wrk)
        use TLab_Transpose
        class(tmpi_transposeX_inner_dt), intent(in) :: self
        real(wp), intent(in) :: a(:)
        real(wp), intent(out) :: b(:)
        real(wp), intent(inout) :: wrk(:)

        ! -----------------------------------------------------------------------
        integer npro, step
        integer(wi) j, l, m, ns, nr, ips, ipr
        integer(wi) count, nmax, ip

        ! #######################################################################
        npro = size(self%send%disp(:))
        step = self%size_block_processes

        ! Make first index last
        count = self%send%disp(2)   ! # all processors transfer same amount
        nmax = count/self%nlines ! size along direction

        do j = 1, npro, step
            l = 0
            do m = j, min(j + step - 1, npro)
                ns = self%send%map(m) + 1; ips = ns - 1
                nr = self%recv%map(m) + 1; ipr = nr - 1

                ! Make first index last
                ip = self%send%disp(ns) + 1
                call TLab_Transpose_Real(a(ip:), nmax, self%nlines, nmax, wrk(ip:), self%nlines)

                l = l + 1
                call MPI_ISEND(wrk(self%send%disp(ns) + 1), self%send%count, self%send%type, ips, ims_tag, self%comm, request(l), ims_err)
                l = l + 1
                call MPI_IRECV(b(self%recv%disp(nr) + 1), self%recv%count, self%recv%type, ipr, ims_tag, self%comm, request(l), ims_err)

            end do
            call MPI_WAITALL(l, request, status, ims_err)
        end do

        return
    end subroutine

    subroutine tmpi_trpX_inner_bwd_real(self, a, b, wrk)
        use TLab_Transpose
        class(tmpi_transposeX_inner_dt), intent(in) :: self
        real(wp), intent(in) :: a(:)
        real(wp), intent(out) :: b(:)
        real(wp), intent(inout) :: wrk(:)

        integer ib, nblocks
        integer(wi) count, nmax, ip

        call tmpi_trpX_real(a, self%recv, wrk, self%send, &
                            self%comm, self%size_block_processes)

        ! Make last index first
        nblocks = size(self%send%disp(:)) ! # of blocks = # of processors
        count = self%send%disp(2)   ! # all processors transfer same amount
        nmax = count/self%nlines ! size along direction

        do ib = 1, nblocks
            ip = (ib - 1)*count + 1
            call TLab_Transpose_Real(wrk(ip:), self%nlines, nmax, self%nlines, b(ip:), nmax)
        end do

        return
    end subroutine

    ! ######################################################################
    ! ######################################################################
    subroutine tmpi_trpX_outer_fwd_complex(self, a, b, wrk)
        class(tmpi_transposeX_outer_dt), intent(in) :: self
        complex(wp), intent(in) :: a(:)
        complex(wp), intent(out) :: b(:)
        complex(wp), intent(inout) :: wrk(:)

        integer ib, nblocks
        integer(wi) count, nmax, ip1, ip2

        nblocks = size(self%send%disp(:)) ! # of blocks = # of processors
        count = self%send%disp(2)   ! # all processors transfer same amount
        nmax = count/self%nlines ! size along direction

        do ib = 1, nblocks
            ip1 = (ib - 1)*self%nlines + 1
            ip2 = (ib - 1)*count + 1
            call reduce_fwd(a(ip1:), self%nlines, nblocks, nmax, wrk(ip2:))
        end do

        call tmpi_trpX_complex(wrk, self%send, b, self%recv, &
                               self%comm, self%size_block_processes)

        return
    end subroutine

    subroutine reduce_fwd(a, nlines, nblocks, nmax, b)
        integer(wi), intent(in) :: nlines, nblocks, nmax
        complex(wp), intent(in) :: a(nlines*nblocks, *)
        complex(wp), intent(out) :: b(nlines, nmax)

        integer n

        do n = 1, nmax
            b(1:nlines, n) = a(1:nlines, n)
        end do

        return
    end subroutine

    subroutine tmpi_trpX_outer_bwd_complex(self, a, b, wrk)
        class(tmpi_transposeX_outer_dt), intent(in) :: self
        complex(wp), intent(in) :: a(:)
        complex(wp), intent(out) :: b(:)
        complex(wp), intent(inout) :: wrk(:)

        integer ib, nblocks
        integer(wi) count, nmax, ip1, ip2

        call tmpi_trpX_complex(a, self%recv, wrk, self%send, &
                               self%comm, self%size_block_processes)

        nblocks = size(self%send%disp(:)) ! # of blocks = # of processors
        count = self%send%disp(2)   ! # all processors transfer same amount
        nmax = count/self%nlines ! size along direction

        do ib = 1, nblocks
            ip1 = (ib - 1)*self%nlines + 1
            ip2 = (ib - 1)*count + 1
            call reduce_bwd(wrk(ip2:), self%nlines, nblocks, nmax, b(ip1:))
        end do

        return
    end subroutine

    subroutine reduce_bwd(a, nlines, nblocks, nmax, b)
        integer(wi), intent(in) :: nlines, nblocks, nmax
        complex(wp), intent(in) :: a(nlines, nmax)
        complex(wp), intent(out) :: b(nlines*nblocks, *)

        integer n

        do n = 1, nmax
            b(1:nlines, n) = a(1:nlines, n)
        end do

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine tmpi_trpX_real(in, send, out, recv, &
                              comm, step)
        real(wp), intent(in) :: in(*)
        real(wp), intent(out) :: out(*)
        type(trp_mem_dt), intent(in) :: send, recv
        type(MPI_Comm), intent(in) :: comm
        integer(wi), intent(in) :: step

        ! -----------------------------------------------------------------------
        integer npro
        integer(wi) j, l, m, ns, nr, ips, ipr

        ! #######################################################################
        npro = size(send%disp(:))

        do j = 1, npro, step
            l = 0
            do m = j, min(j + step - 1, npro)
                ns = send%map(m) + 1; ips = ns - 1
                nr = recv%map(m) + 1; ipr = nr - 1
                l = l + 1
                call MPI_ISEND(in(send%disp(ns) + 1), send%count, send%type, ips, ims_tag, comm, request(l), ims_err)
                l = l + 1
                call MPI_IRECV(out(recv%disp(nr) + 1), recv%count, recv%type, ipr, ims_tag, comm, request(l), ims_err)
            end do
            call MPI_WAITALL(l, request, status, ims_err)
        end do

        return
    end subroutine tmpi_trpX_real

    subroutine tmpi_trpX_complex(in, send, out, recv, &
                                 comm, step)
        complex(wp), intent(in) :: in(*)
        complex(wp), intent(out) :: out(*)
        type(trp_mem_dt), intent(in) :: send, recv
        type(MPI_Comm), intent(in) :: comm
        integer(wi), intent(in) :: step

        ! -----------------------------------------------------------------------
        integer npro
        integer(wi) j, l, m, ns, nr, ips, ipr

        ! #######################################################################
        npro = size(send%disp(:))

        do j = 1, npro, step
            l = 0
            do m = j, min(j + step - 1, npro)
                ns = send%map(m) + 1; ips = ns - 1
                nr = recv%map(m) + 1; ipr = nr - 1
                l = l + 1
                call MPI_ISEND(in(send%disp(ns) + 1), send%count, send%type, ips, ims_tag, comm, request(l), ims_err)
                l = l + 1
                call MPI_IRECV(out(recv%disp(nr) + 1), recv%count, recv%type, ipr, ims_tag, comm, request(l), ims_err)
            end do
            call MPI_WAITALL(l, request, status, ims_err)
        end do

        return
    end subroutine tmpi_trpX_complex

end module
