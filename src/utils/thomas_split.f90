#include "tlab_error.h"

! Splitting algorithm of linear solver

module Thomas_Split
    use TLab_Constants, only: wp, wi, small_wp, roundoff_wp
    use TLab_Constants, only: efile, wfile, fmt_r
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
#ifdef USE_MPI
    use TLabMPI_VARS, only: mpi_axis_dt
#endif
    use Thomas
    implicit none
    private

    public :: thomas_split_dt

    public :: ThomasSplit_3_Reduce_Serial       ! For testing the algorithm in serial
    type, public :: data_dt
        real(wp), pointer :: p(:, :)
    end type data_dt

    ! -----------------------------------------------------------------------
    ! type, extends(thomas_dt) :: thomas_split_dt
    type :: thomas_split_dt
        real(wp), allocatable :: L(:, :)
        real(wp), allocatable :: U(:, :)
        procedure(thomas_ice), pointer, nopass :: ptr_solveL
        procedure(thomas_ice), pointer, nopass :: ptr_solveU
        !
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
        procedure :: solveL => thomas_solveL_dt
        procedure :: solveU => thomas_solveU_dt
#ifdef USE_MPI
        procedure :: reduce => thomas_reduce_mpi_dt
#endif
    end type thomas_split_dt

    abstract interface
        subroutine thomas_ice(A, f)
            import wp
            real(wp), intent(in) :: A(:, :)
            real(wp), intent(inout) :: f(:, :)          ! RHS and solution
        end subroutine thomas_ice
    end interface

contains
    !########################################################################
    !########################################################################
    subroutine thomas_initialize_dt(self, lhs, points, block_id, circulant)
        class(thomas_split_dt), intent(out) :: self
        real(wp), intent(inout) :: lhs(:, :)
        integer(wi), intent(in) :: points(:)   ! sequence of splitting points in ascending order
        integer, intent(in) :: block_id
        logical, intent(in) :: circulant

        ! -------------------------------------------------------------------
        integer ndl, nblocks, nsize, k

        !########################################################################
        self%block_id = block_id
        self%circulant = circulant

        nblocks = size(points)
        ! Number of coefficients are nblocks-1 for the tridiagonal case
        ! and we add one for the circulant case, which is managed by last block

        nsize = size(lhs, 1)

        k = self%block_id
        self%nmin = points(mod(k - 2 + nblocks, nblocks) + 1)
        self%nmin = mod(self%nmin, nsize) + 1
        self%nmax = points(k)

        ! -------------------------------------------------------------------
        ! memory management
        nsize = self%nmax - self%nmin + 1

        ndl = size(lhs, 2)
        allocate (self%L(1:nsize, 1:ndl/2))
        allocate (self%U(1:nsize, ndl/2 + 1:ndl))

        if (allocated(self%y)) deallocate (self%y)
        allocate (self%y(1:nsize, nblocks))
        self%y(:, :) = 0.0_wp

        select case (ndl)
        case (3)
            call ThomasSplit_3_Initialize(self, lhs, points)
        end select

        select case (ndl)
        case (3)
            self%ptr_solveL => Thomas3_SolveL
            self%ptr_solveU => Thomas3_SolveU
        end select

        return
    end subroutine thomas_initialize_dt

    !########################################################################
    !########################################################################
    subroutine ThomasSplit_3_Initialize(self, lhs, points)
        class(thomas_split_dt), intent(inout) :: self
        real(wp), intent(inout) :: lhs(:, :)
        integer(wi), intent(in) :: points(:)   ! sequence of splitting points in ascending order

        ! -------------------------------------------------------------------
        integer nblocks, m
        integer(wi) n, nsize
        integer(wi) p, p_plus_1
        character(len=32) str

        real(wp), allocatable :: lhs_loc(:, :), zloc(:)
        real(wp) alpha_0(2), alpha(2)
        real(wp) delta
        real(wp) alpha_previous(2), beta_loc, gamma_loc

        !########################################################################
        nblocks = size(points)

        ! temporary arrays to calculate z_j
        nsize = size(lhs, 1)
        allocate (lhs_loc(nsize, 3), zloc(nsize))

        if (self%circulant) then
            m = nblocks             ! last block has the circulant coefficient

            call Splitting(L=lhs(:, 1:1), &
                           U=lhs(:, 2:3), &
                           p=points(m), &               ! index of splitting point
                           p_plus_1=1, &
                           alpha=alpha_0, &
                           z=zloc)
            if (self%block_id == m) then
                self%alpha(1:2) = alpha_0(1:2)
            end if
            self%y(:, m) = zloc(self%nmin:self%nmax)

            ! Calculate decay index
            do n = 2, nsize
                if (abs(zloc(n)/zloc(1)) < roundoff_wp) exit
                ! print *, n, abs(zloc(n)/zloc(1))
            end do
            write (str, *) n
            call TLab_Write_ASCII(wfile, 'Decay to round-off in selfting algorithm in '//trim(adjustl(str))//' indexes.')
            write (str, fmt_r) abs(zloc(self%nmax - self%nmin + 1)/zloc(1))
            call TLab_Write_ASCII(wfile, 'Truncation error in splitting algorithm equal to '//trim(adjustl(str))//'.')

        end if

        do m = 1, nblocks - 1       ! loop over remaining coefficients
            call Splitting(L=lhs(:, 1:1), &
                           U=lhs(:, 2:3), &
                           p=points(m), &               ! index of splitting point
                           p_plus_1=points(m) + 1, &
                           alpha=alpha, &
                           z=zloc)

            if (self%block_id == m) then
                self%alpha(1:2) = alpha(1:2)
            end if

            self%y(:, m) = zloc(self%nmin:self%nmax)
            if (self%circulant) then
                gamma_loc = alpha_0(1)*zloc(nsize) + alpha_0(2)*zloc(1)
                self%y(:, m) = self%y(:, m) + gamma_loc*self%y(:, nblocks)
            end if
            if (m > 1) then
                p = points(m - 1)
                p_plus_1 = p + 1
                beta_loc = alpha_previous(1)*zloc(p) + &
                           alpha_previous(2)*zloc(p_plus_1)
                self%y(:, m) = self%y(:, m) + &
                               beta_loc*self%y(:, m - 1)
            end if
            alpha_previous(:) = alpha(:)

        end do

        ! -------------------------------------------------------------------
        ! block matrix Am and LU decomposition
        self%L(:, :) = lhs_loc(self%nmin:self%nmax, 1:1)
        self%U(:, :) = lhs_loc(self%nmin:self%nmax, 2:3)

        deallocate (lhs_loc, zloc)

        return

    contains
        subroutine Splitting(L, U, p, p_plus_1, alpha, z)
            real(wp), intent(inout) :: L(:, :), U(:, :)
            integer(wi), intent(in) :: p, p_plus_1
            real(wp), intent(inout) :: alpha(2)
            real(wp), intent(inout) :: z(1, size(L, 1))

#define a(i) L(i,1)
#define b(i) U(i,1)
#define c(i) U(i,2)

            ! Start definition of alpha
            alpha(1) = a(p_plus_1)
            alpha(2) = c(p)

            ! Generate matrix Ak
            b(p) = b(p) - a(p_plus_1)
            a(p_plus_1) = 0.0_wp
            b(p_plus_1) = b(p_plus_1) - c(p)
            c(p) = 0.0_wp

            ! Generate vector zk
            lhs_loc(:, 1) = a(:)
            lhs_loc(:, 2) = b(:)
            lhs_loc(:, 3) = c(:)
            call Thomas_FactorLU_InPlace(lhs_loc(:, 1:1), &
                                         lhs_loc(:, 2:3))

            z(1, :) = 0.0_wp
            z(1, p) = 1.0_wp
            z(1, p_plus_1) = 1.0_wp
            call Thomas3_SolveL(lhs_loc(:, 1:1), z)
            call Thomas3_SolveU(lhs_loc(:, 2:3), z)

            ! Complete definition of alpha
            delta = 1.0_wp + alpha(1)*z(1, p) + alpha(2)*z(1, p_plus_1)
            if (abs(delta) < small_wp) then
                call TLab_Write_ASCII(efile, __FILE__//'. Singular matrix M.')
                call TLab_Stop(DNS_ERROR_THOMAS)
            end if
            alpha(1) = -alpha(1)/delta
            alpha(2) = -alpha(2)/delta

            return
        end subroutine Splitting

    end subroutine ThomasSplit_3_Initialize

    !########################################################################
    !########################################################################
    subroutine thomas_solveL_dt(self, f)
        class(thomas_split_dt), intent(in) :: self
        real(wp), intent(inout) :: f(:, :)

        call self%ptr_solveL(self%L, f)

        return
    end subroutine

    subroutine thomas_solveU_dt(self, f)
        class(thomas_split_dt), intent(in) :: self
        real(wp), intent(inout) :: f(:, :)

        call self%ptr_solveU(self%U, f)

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
