#include "tlab_error.h"

! Parallel diagonal dominant algorithm based on splitting recurrence

module Thomas_Parallel
    use TLab_Constants, only: wp, wi, small_wp, roundoff_wp
    use TLab_Constants, only: efile, wfile, lfile, fmt_r
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
#ifdef USE_MPI
    use TLabMPI_VARS, only: mpi_axis_dt
#endif
    use Thomas_Split
    implicit none
    private

    public :: thomas_parallel_dt

    public :: ThomasSplit_3_Reduce_Serial       ! For testing the algorithm in serial
    type, public :: data_dt
        real(wp), pointer :: p(:, :)
    end type data_dt

    ! -----------------------------------------------------------------------
    type, extends(thomas_split_dt) :: thomas_parallel_dt
        real(wp), allocatable :: y(:, :)
        !
        logical :: circulant = .true.
        integer :: block_id
        integer(wi) :: nmin, nmax
#ifdef USE_MPI
        type(mpi_axis_dt) mpi
#endif
    contains
        procedure :: initialize => thomas_initialize_dt
#ifdef USE_MPI
        procedure :: reduce => thomas_reduce_mpi_dt
#endif
    end type thomas_parallel_dt

contains
    !########################################################################
    !########################################################################
    subroutine thomas_initialize_dt(self, lhs, points, block_id, circulant)
        class(thomas_parallel_dt), intent(out) :: self
        real(wp), intent(in) :: lhs(:, :)
        integer(wi), intent(in) :: points(:)   ! sequence of splitting points in ascending order
        integer, intent(in) :: block_id
        logical, intent(in) :: circulant

        ! -------------------------------------------------------------------
        integer nblocks, npoints
        integer nsize, nmin, nmax, j
        character(len=32) str

        !########################################################################
        npoints = size(points)
        nblocks = npoints + 1
        nsize = size(lhs, 1)

        ! define block in global grid
        !
        !           point=1     point=2                 point=m
        !  +-----------+-----------+-----------+-----------+-----------+
        !     block=1     block=2                            block=m+1
        !
        nmin = 1
        nmax = nsize
        if (block_id > 1) nmin = points(block_id - 1) + 1
        if (block_id < nblocks) nmax = points(block_id)
        nsize = nmax - nmin + 1

        call self%initialize_base(lhs(1:nsize, :))

        self%block_id = block_id
        self%circulant = circulant
        self%nmin = nmin
        self%nmax = nmax

        if (allocated(self%y)) deallocate (self%y)
        if (self%circulant) then
            allocate (self%y(1:nsize, 0:npoints), source=0.0_wp)
        else
            allocate (self%y(1:nsize, 1:npoints), source=0.0_wp)
        end if

        select case (size(lhs, 2))
        case (3)
            call ThomasSplit_3_Initialize(self, lhs, points)

            ! Calculate truncation error; using 1. block as reference
            if (self%block_id == 1) then
                j = self%block_id
                write (str, fmt_r) abs(self%y(1, j)/self%y(self%nmax - self%nmin + 1, j))
                call TLab_Write_ASCII(wfile, 'Truncation error in splitting algorithm equal to '//trim(adjustl(str))//'.')
            end if

        case default
            call TLab_Write_ASCII(efile, __FILE__//'Only tridiagonal case implemented in splitting algorithm.')
            call TLab_Stop(DNS_ERROR_THOMAS)

        end select

        return
    end subroutine thomas_initialize_dt

    !########################################################################
    !########################################################################
    subroutine ThomasSplit_3_Initialize(self, lhs, points)
        class(thomas_parallel_dt), intent(inout) :: self
        real(wp), intent(in) :: lhs(:, :)
        integer(wi), intent(in) :: points(:)   ! sequence of splitting points in ascending order

        ! -------------------------------------------------------------------
        integer j_min, j, npoints, nblocks
        integer(wi) nsize
        integer(wi) p, p_plus_1, p_loc

        real(wp), allocatable :: z_loc(:)
        real(wp), allocatable :: lhs_loc(:, :)
        real(wp), allocatable :: lu_loc(:, :)
        integer, allocatable :: points_extended(:)
        real(wp) beta_loc, gamma_loc

        !########################################################################
        npoints = size(points)
        nblocks = npoints + 1

        ! temporary arrays to calculate z_j
        nsize = size(lhs, 1)

        if (allocated(lhs_loc)) deallocate (lhs_loc)    ! block diagonal matrix
        allocate (lhs_loc, source=lhs)

        if (allocated(lu_loc)) deallocate (lu_loc)      ! lu decomposition of block diagonal matrix
        allocate (lu_loc, mold=lhs)

        if (allocated(z_loc)) deallocate (z_loc)
        allocate (z_loc(nsize))

        ! -------------------------------------------------------------------
        j_min = 1
        if (self%circulant) j_min = 0   ! handling of circulant case

        p_loc = 1                       ! index of starting subarray
        do j = j_min, npoints
            if (j == 0) then            ! handling of circulant case
                p = nsize
            else
                p = points(j)
            end if
            p_plus_1 = mod(p, nsize) + 1

            lu_loc(p_loc:, 1:1) = lhs_loc(p_loc:, 1:1)      ! recover original system for current block
            lu_loc(p_loc:, 2:3) = lhs_loc(p_loc:, 2:3)
            call Thomas_3_Split_InPlace(L=lu_loc(p_loc:, 1:1), &
                                        U=lu_loc(p_loc:, 2:3), &
                                        z=z_loc(p_loc:), &
                                        index=p - p_loc + 1)

            self%y(:, j) = z_loc(self%nmin:self%nmax)

            ! store block matrix A for next splitting
            lhs_loc(p, 2) = lhs_loc(p, 2) - lhs_loc(p, 3)                           ! b_n - c_n
            lhs_loc(p_plus_1, 2) = lhs_loc(p_plus_1, 2) - lhs_loc(p_plus_1, 1)      ! b_1 - a_1

            if (j > j_min) then
                beta_loc = z_loc(p_loc)
                self%y(:, j) = self%y(:, j) + beta_loc*self%y(:, j - 1)

                if (self%circulant) then
                    gamma_loc = z_loc(nsize)
                    self%y(:, j) = self%y(:, j) + gamma_loc*self%y(:, 0)
                end if

            end if

            p_loc = p_plus_1

            if (j == 0) call decay_index(z_loc)

        end do

        ! -------------------------------------------------------------------
        ! block matrix Am and LU decomposition
        self%L(:, :) = lu_loc(self%nmin:self%nmax, 1:1)
        self%U(:, :) = lu_loc(self%nmin:self%nmax, 2:3)

        return
    end subroutine ThomasSplit_3_Initialize

#ifdef USE_MPI
    !########################################################################
    !########################################################################
    subroutine thomas_reduce_mpi_dt(self, f, alpha, tmp)
        use mpi_f08
#ifdef PROFILE_ON
        use TLabMPI_VARS, only: ims_time_trans
#endif
        class(thomas_parallel_dt), intent(in) :: self
        real(wp), intent(inout) :: f(:, :)
        real(wp), intent(inout) :: alpha(:)         ! auxiliary memory space for local alpha
        real(wp), intent(inout) :: tmp(:)           ! auxiliary memory space for all alphas
        integer(wi) n, nsize, nlines
        integer(wi) nblocks

        integer ims_err
        integer source, dest, tag

#ifdef PROFILE_ON
        real(wp) time_loc_1, time_loc_2
#endif

        !########################################################################
        ! Assume circulant matrix and need alpha_0

        nblocks = self%mpi%num_processors
        nlines = size(f, 1)
        nsize = size(f, 2)              ! Assume all blocks have same size

#ifdef PROFILE_ON
        time_loc_1 = MPI_WTIME()
#endif

        if (self%mpi%num_processors == 1 .and. self%circulant) then
            call Thomas_3_Split_Reduce(self%L, &
                                       self%U, &
                                       self%y(:, 0), &
                                       f, alpha, size(self%L, 1))
            return
        end if

        ! -------------------------------------------------------------------
        ! pass x(:,1) to previous block and calculate local coefficient
#define xp(j) alpha(j)
        dest = mod(self%mpi%rank - 1 + self%mpi%num_processors, self%mpi%num_processors)
        source = mod(self%mpi%rank + 1, self%mpi%num_processors)
        tag = 0
        call MPI_Sendrecv(f(:, 1), nlines, MPI_REAL8, dest, tag, &
                          xp(:), nlines, MPI_REAL8, source, tag, &
                          self%mpi%comm, MPI_STATUS_IGNORE, ims_err)
        alpha(:) = f(:, nsize) + xp(:)
!
#undef xp

        ! -------------------------------------------------------------------
        ! Truncated algorithm
        ! Pass alpha to following block
        dest = mod(self%mpi%rank + 1, self%mpi%num_processors)
        source = mod(self%mpi%rank - 1 + self%mpi%num_processors, self%mpi%num_processors)
        tag = 1
        call MPI_Sendrecv(alpha, nlines, MPI_REAL8, dest, tag, &
                          tmp, nlines, MPI_REAL8, source, tag, &
                          self%mpi%comm, MPI_STATUS_IGNORE, ims_err)

#ifdef PROFILE_ON
        time_loc_2 = MPI_WTIME()
        ims_time_trans = ims_time_trans + (time_loc_2 - time_loc_1)
#endif

        ! Update solution
        do n = 1, nsize
            f(:, n) = f(:, n) + alpha(:)*self%y(n, dest) &
                      + tmp(:)*self%y(n, self%block_id - 1)
        end do

        ! -------------------------------------------------------------------
        ! Full algorithm
        ! Distribute coefficients
        ! call MPI_Allgather(alpha, nlines, MPI_REAL8, &
        !                    tmp, nlines, MPI_REAL8, self%communicator, ims_err)

        ! ! Update solution
        ! do m = 1, nblocks
        !     mm = mod(m, nblocks)
        !     do n = 1, nsize
        !         f(:, n) = f(:, n) + tmp(:, m)*self%y(n, mm)
        !     end do
        ! end do

        return
    end subroutine thomas_reduce_mpi_dt
#endif

    !########################################################################
    !########################################################################
    subroutine ThomasSplit_3_Reduce_Serial(self, f)
        type(thomas_parallel_dt), intent(in) :: self(:)
        type(data_dt), intent(inout) :: f(:)

        integer(wi) n, nsize, nlines
        integer(wi) m, mmax
        integer(wi) k, nblocks, k_plus_1
        real(wp), allocatable :: xp(:, :), alpha(:, :)

        !########################################################################
        nblocks = size(self)
        nlines = size(f(1)%p, 1)

        ! if (nblocks == 1) then
        !     call Thomas_3_Split_Reduce(self(1)%L, &
        !                                self(1)%U, &
        !                                self(1)%y(:, 0), &
        !                                f(1)%p(:, :), alpha, size(self(1)%L, 1))
        !     return
        ! end if

        ! pass x(:,1) to previous block
        allocate (xp(nlines, nblocks))
        do k = nblocks, 1, -1                   ! loop over blocks
            k_plus_1 = mod(k + nblocks, nblocks) + 1
            xp(:, k) = f(k_plus_1)%p(:, 1)      ! solution at left boundary of block k+1
        end do

        !  calculate local alpha
        allocate (alpha(nlines, nblocks))       ! the idea is that each block needs only one alpha
        do k = nblocks, 1, -1                   ! loop over blocks
            nsize = size(f(k)%p, 2)
            alpha(:, k) = f(k)%p(:, nsize) + xp(:, k)
        end do

        ! send alpha to all blocks

        ! calculate truncated solution
        do k = nblocks, 1, -1                   ! loop over blocks
            nsize = size(f(k)%p, 2)

            ! Truncated
            ! m = k
            m = mod(k, nblocks)
            do n = 1, nsize
                f(k)%p(:, n) = f(k)%p(:, n) + alpha(:, k)*self(k)%y(n, m)
            end do
            m = mod(k - 2 + nblocks, nblocks) + 1
            do n = 1, nsize
                f(k)%p(:, n) = f(k)%p(:, n) + alpha(:, m)*self(k)%y(n, k - 1)
            end do

        end do

        ! calculate full solution
        ! do k = nblocks, 1, -1                   ! loop over blocks
        !     nsize = size(f(k)%p, 2)

        !     if (self(k)%circulant) then
        !         mmax = nblocks
        !     else
        !         mmax = nblocks - 1
        !     end if
        !     do m = 1, mmax
        !         do n = 1, nsize
        !             f(k)%p(:, n) = f(k)%p(:, n) + alpha(:, m)*self(k)%y(n, m)
        !         end do
        !     end do

        ! end do

        deallocate (xp, alpha)

        return
    end subroutine ThomasSplit_3_Reduce_Serial

end module Thomas_Parallel
