#include "tlab_error.h"

! Splitting algorithm of linear solver

module Thomas_Split
    use TLab_Constants, only: wp, wi, small_wp, roundoff_wp
    use TLab_Constants, only: efile, wfile, lfile, fmt_r
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
#ifdef USE_MPI
    use TLabMPI_VARS, only: mpi_axis_dt
#endif
    use Thomas
    implicit none
    private

    public :: thomas_split_dt

    public :: Thomas_3_Split_InPlace

    public :: ThomasSplit_3_Reduce_Serial       ! For testing the algorithm in serial
    type, public :: data_dt
        real(wp), pointer :: p(:, :)
    end type data_dt

    ! -----------------------------------------------------------------------
    type, extends(thomas_base_dt) :: thomas_split_dt
        real(wp), allocatable :: y(:, :)
        !
        logical :: circulant = .true.
        integer :: block_id
        integer(wi) :: nmin, nmax
        real(wp) :: alpha(2) = [0.0_wp, 0.0_wp]
#ifdef USE_MPI
        type(mpi_axis_dt) mpi
#endif
    contains
        procedure :: initialize => thomas_initialize_dt
#ifdef USE_MPI
        procedure :: reduce => thomas_reduce_mpi_dt
#endif
    end type thomas_split_dt

contains
    !########################################################################
    !########################################################################
    subroutine thomas_initialize_dt(self, lhs, points, block_id, circulant)
        class(thomas_split_dt), intent(out) :: self
        real(wp), intent(in) :: lhs(:, :)
        integer(wi), intent(in) :: points(:)   ! sequence of splitting points in ascending order
        integer, intent(in) :: block_id
        logical, intent(in) :: circulant

        ! -------------------------------------------------------------------
        integer nblocks, nsize, nmin, nmax

        !########################################################################
        ! Number of coefficients are nblocks-1 for the tridiagonal case
        ! and we add one for the circulant case, which is managed by last block
        nblocks = size(points)
        nsize = size(lhs, 1)

        nmin = points(mod(block_id - 2 + nblocks, nblocks) + 1)
        nmin = mod(nmin, nsize) + 1
        nmax = points(block_id)
        nsize = nmax - nmin + 1

        call self%initialize_base(lhs(1:nsize, :))

        self%block_id = block_id
        self%circulant = circulant
        self%nmin = nmin
        self%nmax = points(block_id)

        if (allocated(self%y)) deallocate (self%y)
        allocate (self%y(1:nsize, nblocks), source=0.0_wp)

        select case (size(lhs, 2))
        case (3)
            call ThomasSplit_3_Initialize(self, lhs, points)

        case default
            call TLab_Write_ASCII(efile, __FILE__//'Only tridiagonal case implemented in splitting algorithm.')
            call TLab_Stop(DNS_ERROR_THOMAS)

        end select

        return
    end subroutine thomas_initialize_dt

    !########################################################################
    !########################################################################
    subroutine ThomasSplit_3_Initialize(self, lhs, points)
        class(thomas_split_dt), intent(inout) :: self
        real(wp), intent(in) :: lhs(:, :)
        integer(wi), intent(in) :: points(:)   ! sequence of splitting points in ascending order

        ! -------------------------------------------------------------------
        integer nblocks, m
        integer(wi) nsize
        integer(wi) p, p_plus_1, p_loc

        real(wp), allocatable :: z_loc(:)
        real(wp), allocatable :: lhs_loc(:, :)
        real(wp), allocatable :: lu_loc(:, :)
        real(wp) alpha_0(2)
        real(wp) alpha_previous(2), beta_loc, gamma_loc
        character(len=32) str

        !########################################################################
        nblocks = size(points)

        ! temporary arrays to calculate z_j
        nsize = size(lhs, 1)

        if (allocated(lhs_loc)) deallocate (lhs_loc)    ! block diagonal matrix
        allocate (lhs_loc, source=lhs)

        if (allocated(lu_loc)) deallocate (lu_loc)      ! lu decompsition of block diagonal matrix
        allocate (lu_loc, mold=lhs)

        if (allocated(z_loc)) deallocate (z_loc)
        allocate (z_loc(nsize))

        if (self%circulant) then
            m = nblocks                                 ! last block has the circulant coefficient

            lu_loc(:, 1:1) = lhs_loc(:, 1:1)            ! recover original system for current block
            lu_loc(:, 2:3) = lhs_loc(:, 2:3)
            call Thomas_3_Split_InPlace(L=lu_loc(:, 1:1), &
                                        U=lu_loc(:, 2:3), &
                                        z=z_loc, &
                                        index=points(m))

            self%y(:, m) = z_loc(self%nmin:self%nmax)

            alpha_0(1) = lu_loc(1, 1)            ! a1
            alpha_0(2) = lu_loc(points(m), 3)    ! cn
            if (self%block_id == m) then
                self%alpha(1:2) = alpha_0(1:2)
            end if

            ! store block matrix A
            lhs_loc(1, 2) = lhs_loc(1, 2) - lhs_loc(nsize, 3)
            lhs_loc(nsize, 2) = lhs_loc(nsize, 2) - lhs_loc(1, 1)

            ! Calculate decay index
            call decay_index(z_loc)
            write (str, fmt_r) abs(z_loc(self%nmax - self%nmin + 1)/z_loc(1))
            call TLab_Write_ASCII(wfile, 'Truncation error in splitting algorithm equal to '//trim(adjustl(str))//'.')

        end if

        p_loc = 1
        do m = 1, nblocks - 1                               ! loop over remaining coefficients
            lu_loc(p_loc:, 1:1) = lhs_loc(p_loc:, 1:1)      ! recover original system for current block
            lu_loc(p_loc:, 2:3) = lhs_loc(p_loc:, 2:3)
            if (m > 1) z_loc(:p_loc - 1) = 0.0_wp
            call Thomas_3_Split_InPlace(L=lu_loc(p_loc:, 1:1), &
                                        U=lu_loc(p_loc:, 2:3), &
                                        z=z_loc(p_loc:), &
                                        index=points(m) - p_loc + 1)

            self%y(:, m) = z_loc(self%nmin:self%nmax)

            if (self%circulant) then
                gamma_loc = alpha_0(1)*z_loc(nsize) + alpha_0(2)*z_loc(1)
                self%y(:, m) = self%y(:, m) + gamma_loc*self%y(:, nblocks)
            end if

            if (m > 1) then
                p = points(m - 1)
                p_plus_1 = p + 1
                beta_loc = alpha_previous(1)*z_loc(p) + &
                           alpha_previous(2)*z_loc(p_plus_1)
                self%y(:, m) = self%y(:, m) + &
                               beta_loc*self%y(:, m - 1)
            end if

            alpha_previous(1) = lu_loc(points(m) + 1, 1)    ! a_p_plus_1
            alpha_previous(2) = lu_loc(points(m), 3)        ! c_p
            if (self%block_id == m) then
                self%alpha(:) = alpha_previous(:)
            end if

            ! store block matrix A
            lhs_loc(points(m), 2) = lhs_loc(points(m), 2) - lhs_loc(points(m) + 1, 1)
            lhs_loc(points(m) + 1, 2) = lhs_loc(points(m) + 1, 2) - lhs_loc(points(m), 3)

            p_loc = points(m) + 1

        end do

        ! -------------------------------------------------------------------
        ! block matrix Am and LU decomposition
        self%L(:, :) = lu_loc(self%nmin:self%nmax, 1:1)
        self%U(:, :) = lu_loc(self%nmin:self%nmax, 2:3)

        return
    end subroutine ThomasSplit_3_Initialize

    !########################################################################
    !########################################################################
    subroutine Thomas_3_Split_InPlace(L, U, z, index)
        real(wp), intent(inout) :: L(:, :), U(:, :)
        real(wp), intent(out) :: z(1, size(L, 1))
        integer, intent(in) :: index

        ! -----------------------------------------------------------------------
        integer nmax, p, p_plus_1
        real(wp) alpha(2), delta

        ! #######################################################################
        nmax = size(L, 1)

        p = mod(index - 1, nmax) + 1                ! in circulant cases, this n
        p_plus_1 = mod(p + 1 - 1, nmax) + 1         ! in circulant cases, this 1

#define a(i) L(i,1)
#define b(i) U(i,1)
#define c(i) U(i,2)

        ! Start definition of alpha; in circulant cases, this is a1 and cn
        alpha(1) = a(p_plus_1)
        alpha(2) = c(p)

        ! Generate matrix A1
        b(p) = b(p) - a(p_plus_1)
        a(p_plus_1) = 0.0_wp
        b(p_plus_1) = b(p_plus_1) - c(p)
        c(p) = 0.0_wp

        ! call Thomas3_FactorLU_InPlace(L, U)
        call Thomas_FactorLU_InPlace(L, U)

        ! Generate vector z1
        z(1, :) = 0.0_wp
        z(1, p) = 1.0_wp
        z(1, p_plus_1) = 1.0_wp

        call Thomas3_SolveL(L, z)
        call Thomas3_SolveU(U, z)

        ! Calculate normalized alpha coefficients
        delta = 1.0_wp + alpha(1)*z(1, p) + alpha(2)*z(1, p_plus_1)
        if (abs(delta) < small_wp) then
            call TLab_Write_ASCII(efile, __FILE__//'. Singular matrix M.')
            call TLab_Stop(DNS_ERROR_THOMAS)
        end if

        a(p_plus_1) = -alpha(1)/delta
        c(p) = -alpha(2)/delta

        ! -------------------------------------------------------------------
        ! Calculate decay index
        ! call decay_index(z(1, p_plus_1:))!, n_smw_decay)

#undef a
#undef b
#undef c

        return
    end subroutine Thomas_3_Split_InPlace

    !########################################################################
    !########################################################################
    subroutine decay_index(z, index)
        use TLab_Constants, only: roundoff_wp
        real(wp), intent(in) :: z(:)
        integer, intent(out), optional :: index

        integer n
        character(len=32) str

        do n = 2, size(z)
            if (abs(z(n)/z(1)) < roundoff_wp) exit
            ! print *, abs(z(n)/z(1)
        end do
        write (str, *) n
        call TLab_Write_ASCII(lfile, 'Decay to round-off in SMW algorithm in '//trim(adjustl(str))//' indexes.')

        if (present(index)) index = n

        return
    end subroutine

#ifdef USE_MPI
    !########################################################################
    !########################################################################
    subroutine thomas_reduce_mpi_dt(self, f, alpha, tmp)
        use mpi_f08
        class(thomas_split_dt), intent(in) :: self
        real(wp), intent(inout) :: f(:, :)
        real(wp), intent(inout) :: alpha(:)         ! auxiliary memory space for local alpha
        real(wp), intent(inout) :: tmp(:)           ! auxiliary memory space
        integer(wi) n, nsize, nlines
        integer(wi) nblocks

        integer ims_err
        integer source, dest, tag

        !########################################################################
        ! Assume circulant matrix and need alpha_0

        nblocks = self%mpi%num_processors
        nlines = size(f, 1)
        nsize = size(f, 2)              ! Assume all blocks have same size

        ! -------------------------------------------------------------------
        ! pass x(:,1) to previous block and calculate local coefficient
#define xp(j) alpha(j)
        dest = mod(self%mpi%rank - 1 + self%mpi%num_processors, self%mpi%num_processors)
        source = mod(self%mpi%rank + 1, self%mpi%num_processors)
        tag = 0
        call MPI_Sendrecv(f(:, 1), nlines, MPI_REAL8, dest, tag, &
                          xp(:), nlines, MPI_REAL8, source, tag, &
                          self%mpi%comm, MPI_STATUS_IGNORE, ims_err)
        alpha(:) = self%alpha(1)*f(:, nsize) + self%alpha(2)*xp(:)

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

        ! Update solution
        do n = 1, nsize
            f(:, n) = f(:, n) + alpha(:)*self%y(n, self%block_id) &
                      + tmp(:)*self%y(n, source + 1)
        end do

        ! -------------------------------------------------------------------
        ! Full algorithm
        ! Distribute coefficients
        ! call MPI_Allgather(alpha, nlines, MPI_REAL8, &
        !                    tmp, nlines, MPI_REAL8, self%communicator, ims_err)

        ! ! Update solution
        ! do m = 1, nblocks
        !     do n = 1, nsize
        !         f(:, n) = f(:, n) + tmp(:, m)*self%y(n, m)
        !     end do
        ! end do

        return
    end subroutine thomas_reduce_mpi_dt
#endif

    !########################################################################
    !########################################################################
    subroutine ThomasSplit_3_Reduce_Serial(self, f)
        type(thomas_split_dt), intent(in) :: self(:)
        type(data_dt), intent(inout) :: f(:)

        integer(wi) n, nsize, nlines
        integer(wi) m, mmax
        integer(wi) k, nblocks, k_plus_1
        real(wp), allocatable :: xp(:, :), alpha(:, :)

        !########################################################################
        nblocks = size(self)
        nlines = size(f(1)%p, 1)

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
            alpha(:, k) = self(k)%alpha(1)*f(k)%p(:, nsize) + self(k)%alpha(2)*xp(:, k)
        end do

        ! send alpha to all blocks

        ! calculate solution
        do k = nblocks, 1, -1                   ! loop over blocks
            nsize = size(f(k)%p, 2)
            if (self(k)%circulant) then
                mmax = nblocks
            else
                mmax = nblocks - 1
            end if

            ! Truncated
            m = k
            do n = 1, nsize
                f(k)%p(:, n) = f(k)%p(:, n) + alpha(:, m)*self(k)%y(n, m)
            end do
            m = mod(k - 2 + nblocks, nblocks) + 1
            do n = 1, nsize
                f(k)%p(:, n) = f(k)%p(:, n) + alpha(:, m)*self(k)%y(n, m)
            end do

            ! Full
            ! do m = 1, mmax
            !     do n = 1, nsize
            !         f(k)%p(:, n) = f(k)%p(:, n) + alpha(:, m)*self(k)%y(n, m)
            !     end do
            ! end do

        end do

        deallocate (xp, alpha)

        return
    end subroutine ThomasSplit_3_Reduce_Serial

end module Thomas_Split
