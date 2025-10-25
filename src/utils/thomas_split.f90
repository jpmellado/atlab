#include "tlab_error.h"

! Matrix splitting

module Thomas_Split
    use TLab_Constants, only: wp, wi, small_wp, roundoff_wp
    use TLab_Constants, only: efile, lfile, fmt_r
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
#ifdef USE_MPI
    use mpi_f08, only: MPI_Comm
#endif
    use Thomas
    implicit none
    private

    public :: Thomas3_Split_Initialize
    public :: Thomas3_Split_Solve_Serial
#ifdef USE_MPI
    public :: Thomas3_Split_Solve_MPI
#endif

    type, public :: thomas3_split_dt
#ifdef USE_MPI
        type(MPI_Comm) communicator
        integer rank, n_ranks
#endif
        integer :: block_id
        integer(wi) :: nmin, nmax
        logical :: circulant = .true.
        real(wp), allocatable :: lhs(:, :)
        real(wp), allocatable :: y(:, :)
        real(wp) :: alpha(2) = [0.0_wp, 0.0_wp]
    end type thomas3_split_dt

    type, public :: data_dt                     ! for validating the algorithm in serial
        real(wp), pointer :: p(:, :)
    end type data_dt

contains
    !########################################################################
    !########################################################################
    subroutine Thomas3_Split_Initialize(L, U, points, split)
        real(wp), intent(inout) :: L(:, :), U(:, :)
        integer(wi), intent(in) :: points(:)   ! sequence of splitting points in ascending order
        type(thomas3_split_dt), intent(inout) :: split

        ! -------------------------------------------------------------------
        integer nblocks, m, k
        integer(wi) n, nsize
        integer(wi) p, p_plus_1
        character(len=32) str

        ! real(wp), allocatable :: aloc(:), bloc(:), cloc(:), zloc(:)
        real(wp), allocatable :: lhs_loc(:, :), zloc(:)
        real(wp) alpha_0(2), alpha(2)
        real(wp) delta
        real(wp) alpha_previous(2), beta_loc, gamma_loc

#define a(i) L(i,1)
#define b(i) U(i,1)
#define c(i) U(i,2)

        !########################################################################
        nblocks = size(points)
        ! Number of coefficients are nblocks-1 for the tridiagonal case
        ! and we add one for the circulant case, which is managed by last block

        nsize = size(a(:))

        k = split%block_id
        split%nmin = points(mod(k - 2 + nblocks, nblocks) + 1)
        split%nmin = mod(split%nmin, nsize) + 1
        split%nmax = points(k)

        ! -------------------------------------------------------------------
        ! memory management
        nsize = split%nmax - split%nmin + 1

        if (allocated(split%lhs)) deallocate (split%lhs)
        allocate (split%lhs(1:nsize, 3))
        split%lhs(:, :) = 0.0_wp

        if (allocated(split%y)) deallocate (split%y)
        allocate (split%y(1:nsize, nblocks))
        split%y(:, :) = 0.0_wp

        ! -------------------------------------------------------------------
        ! temporary arrays to calculate z_j
        nsize = size(a(:))
        ! allocate (aloc(nsize), bloc(nsize), cloc(nsize), zloc(nsize))
        allocate (lhs_loc(nsize, 3), zloc(nsize))

        if (split%circulant) then
            m = nblocks             ! last block has the circulant coefficient
            p = points(m)           ! index of splitting point
            p_plus_1 = 1

            ! call Splitting(a, b, c, p, p_plus_1, alpha_0, zloc)
            call Splitting(L, U, p, p_plus_1, alpha_0, zloc)

            if (split%block_id == m) then
                split%alpha(1:2) = alpha_0(1:2)
            end if
            split%y(:, m) = zloc(split%nmin:split%nmax)

            ! Calculate decay index
            do n = 2, nsize
                if (abs(zloc(n)/zloc(1)) < roundoff_wp) exit
                ! print *, n, abs(zloc(n)/zloc(1))
            end do
            write (str, *) n
            call TLab_Write_ASCII(lfile, 'Decay to round-off in splitting algorithm in '//trim(adjustl(str))//' indexes.')
            write (str, fmt_r) abs(zloc(split%nmax - split%nmin + 1)/zloc(1))
            call TLab_Write_ASCII(lfile, 'Truncation error in splitting algorithm equal to '//trim(adjustl(str))//'.')
        end if

        do m = 1, nblocks - 1       ! loop over remaining coefficients
            p = points(m)           ! index of splitting point
            p_plus_1 = p + 1

            ! call Splitting(a, b, c, p, p_plus_1, alpha, zloc)
            call Splitting(L, U, p, p_plus_1, alpha, zloc)

            if (split%block_id == m) then
                split%alpha(1:2) = alpha(1:2)
            end if

            split%y(:, m) = zloc(split%nmin:split%nmax)
            if (split%circulant) then
                gamma_loc = alpha_0(1)*zloc(nsize) + alpha_0(2)*zloc(1)
                split%y(:, m) = split%y(:, m) + gamma_loc*split%y(:, nblocks)
            end if
            if (m > 1) then
                p = points(m - 1)
                p_plus_1 = p + 1
                beta_loc = alpha_previous(1)*zloc(p) + alpha_previous(2)*zloc(p_plus_1)
                split%y(:, m) = split%y(:, m) + beta_loc*split%y(:, m - 1)
            end if
            alpha_previous(:) = alpha(:)

        end do

        ! -------------------------------------------------------------------
        ! block matrix Am and LU decomposition
        ! split%lhs(:, 1) = aloc(split%nmin:split%nmax)
        ! split%lhs(:, 2) = bloc(split%nmin:split%nmax)
        ! split%lhs(:, 3) = cloc(split%nmin:split%nmax)
        split%lhs(:, 1) = lhs_loc(split%nmin:split%nmax, 1)
        split%lhs(:, 2) = lhs_loc(split%nmin:split%nmax, 2)
        split%lhs(:, 3) = lhs_loc(split%nmin:split%nmax, 3)

        deallocate (lhs_loc, zloc)

        return

    contains
        ! subroutine Splitting(a, b, c, p, p_plus_1, alpha, z)
        !     real(wp), intent(inout) :: a(:), b(:), c(:)
        subroutine Splitting(L, U, p, p_plus_1, alpha, z)
            real(wp), intent(inout) :: L(:, :), U(:, :)
            integer(wi), intent(in) :: p, p_plus_1
            real(wp), intent(inout) :: alpha(2)
            real(wp), intent(inout) :: z(1, size(L, 1))

            ! Start definition of alpha
            alpha(1) = a(p_plus_1)
            alpha(2) = c(p)

            ! Generate matrix Ak
            b(p) = b(p) - a(p_plus_1)
            a(p_plus_1) = 0.0_wp
            b(p_plus_1) = b(p_plus_1) - c(p)
            c(p) = 0.0_wp

            ! Generate vector zk
            ! aloc(:) = a(:); bloc(:) = b(:); cloc(:) = c(:)
            lhs_loc(:, 1) = a(:)
            lhs_loc(:, 2) = b(:)
            lhs_loc(:, 3) = c(:)
            ! call Thomas3_FactorLU(size(aloc), aloc, bloc, cloc)
            call Thomas_FactorLU_InPlace(lhs_loc(:, 1:1), &
                                         lhs_loc(:, 2:3))

            z(1, :) = 0.0_wp
            z(1, p) = 1.0_wp
            z(1, p_plus_1) = 1.0_wp
            ! call Thomas3_SolveLU(size(aloc), 1, aloc, bloc, cloc, z(:))
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

    end subroutine Thomas3_Split_Initialize

    !########################################################################
    !########################################################################
    subroutine Thomas3_Split_Solve_Serial(split, f)
        type(thomas3_split_dt), intent(in) :: split(:)
        type(data_dt), intent(inout) :: f(:)

        integer(wi) n, nsize, nlines
        integer(wi) m, mmax
        integer(wi) k, nblocks, k_plus_1
        real(wp), allocatable :: xp(:, :), alpha(:, :)

        !########################################################################
        nblocks = size(split)
        nlines = size(f(1)%p, 1)

        ! Solving block system Am in each block
        ! this part involves lhs and could be out of this routine.
        ! this routine is then called Thomas3_Split_Reduction
        do k = 1, nblocks                       ! loop over blocks
            ! nsize = split(k)%nmax - split(k)%nmin + 1
            nsize = size(f(k)%p, 2)
            ! call Thomas3_SolveLU(nsize, nlines, &
            !                      split(k)%lhs(:, 1), split(k)%lhs(:, 2), split(k)%lhs(:, 3), f(k)%p(:, :))
            call Thomas3_SolveL(split(k)%lhs(:, 1:1), f(k)%p(:, :))
            call Thomas3_SolveU(split(k)%lhs(:, 2:3), f(k)%p(:, :))
        end do

        ! -------------------------------------------------------------------
        ! Reduction step

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
            alpha(:, k) = split(k)%alpha(1)*f(k)%p(:, nsize) + split(k)%alpha(2)*xp(:, k)
        end do

        ! send alpha to all blocks

        ! calculate solution
        do k = nblocks, 1, -1                   ! loop over blocks
            nsize = size(f(k)%p, 2)
            if (split(k)%circulant) then
                mmax = nblocks
            else
                mmax = nblocks - 1
            end if

            ! Truncated
            m = k
            do n = 1, nsize
                f(k)%p(:, n) = f(k)%p(:, n) + alpha(:, m)*split(k)%y(n, m)
            end do
            m = mod(k - 2 + nblocks, nblocks) + 1
            do n = 1, nsize
                f(k)%p(:, n) = f(k)%p(:, n) + alpha(:, m)*split(k)%y(n, m)
            end do

            ! Full
            ! do m = 1, mmax
            !     do n = 1, nsize
            !         f(k)%p(:, n) = f(k)%p(:, n) + alpha(:, m)*split(k)%y(n, m)
            !     end do
            ! end do

        end do

        deallocate (xp, alpha)

        return
    end subroutine Thomas3_Split_Solve_Serial

#ifdef USE_MPI
    !########################################################################
    !########################################################################
    subroutine Thomas3_Split_Solve_MPI(split, f, alpha, tmp)
        use mpi_f08

        type(thomas3_split_dt), intent(in) :: split
        real(wp), intent(inout) :: f(:, :)
        real(wp), intent(inout) :: alpha(:)         ! auxiliary memory space for local alpha
        real(wp), intent(inout) :: tmp(:)           ! auxiliary memory space
        integer(wi) n, nsize, nlines
        integer(wi) m
        integer(wi) k, nblocks

        integer ims_err
        integer source, dest, tag

        !########################################################################
        nblocks = split%n_ranks
        nlines = size(f, 1)
        nsize = size(f, 2)              ! Assume all blocks have same size

        !########################################################################
        ! Solving block system Am in each block
        ! call Thomas3_SolveLU(nsize, nlines, &
        !                      split%lhs(:, 1), split%lhs(:, 2), split%lhs(:, 3), f(:, :))
        call Thomas3_SolveL(split%lhs(:, 1:1), f)
        call Thomas3_SolveU(split%lhs(:, 2:3), f)

        !########################################################################
        ! Reduction step
        ! Assume circulant matrix and need alpha_0

        ! -------------------------------------------------------------------
        ! pass x(:,1) to previous block and calculate local coefficient
#define xp(j) alpha(j)

        dest = mod(split%rank - 1 + split%n_ranks, split%n_ranks)
        source = mod(split%rank + 1, split%n_ranks)
        tag = 0
        call MPI_Sendrecv(f(:, 1), nlines, MPI_REAL8, dest, tag, &
                          xp(:), nlines, MPI_REAL8, source, tag, &
                          split%communicator, MPI_STATUS_IGNORE, ims_err)
        alpha(:) = split%alpha(1)*f(:, nsize) + split%alpha(2)*xp(:)

#undef xp

        ! -------------------------------------------------------------------
        ! Truncated algorithm
        ! Pass alpha to following block
        dest = mod(split%rank + 1, split%n_ranks)
        source = mod(split%rank - 1 + split%n_ranks, split%n_ranks)
        tag = 1
        call MPI_Sendrecv(alpha, nlines, MPI_REAL8, dest, tag, &
                          tmp, nlines, MPI_REAL8, source, tag, &
                          split%communicator, MPI_STATUS_IGNORE, ims_err)

        ! Update solution
        do n = 1, nsize
            f(:, n) = f(:, n) + alpha(:)*split%y(n, split%block_id) &
                      + tmp(:)*split%y(n, source + 1)
        end do

        ! -------------------------------------------------------------------
        ! Full algorithm
        ! Distribute coefficients
        ! call MPI_Allgather(alpha, nlines, MPI_REAL8, &
        !                    tmp, nlines, MPI_REAL8, split%communicator, ims_err)

        ! ! Update solution
        ! do m = 1, nblocks
        !     do n = 1, nsize
        !         f(:, n) = f(:, n) + tmp(:, m)*split%y(n, m)
        !     end do
        ! end do

        return
    end subroutine Thomas3_Split_Solve_MPI
#endif

end module Thomas_Split
