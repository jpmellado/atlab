#include "tlab_error.h"

! Matrix splitting

module Thomas3_Split
    use TLab_Constants, only: wp, wi, small_wp, efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Thomas3
    implicit none
    private

    public :: Thomas3_Split_Initialize_Global
    public :: Thomas3_Split_Solve_Global

    type, public :: matrix_split_dt
        integer :: np                           ! number of splitting points
        logical :: periodic = .false.
        integer(wi), allocatable :: points(:)   ! sequence of splitting points in ascending order
        real(wp), allocatable :: z(:, :)
        real(wp), allocatable :: alpha(:, :)
        real(wp), allocatable :: beta(:)
        real(wp), allocatable :: gamma(:)
    end type matrix_split_dt

    type, public :: thomas3_split_dt
        integer :: block
        integer(wi) :: pmin, pmax
        logical :: periodic = .false.
        real(wp), allocatable :: a(:), b(:), c(:)
        real(wp), allocatable :: z(:, :)
        real(wp), allocatable :: alpha(:, :)
        real(wp), allocatable :: beta(:)
        real(wp), allocatable :: gamma(:)
    end type thomas3_split_dt

contains
    !########################################################################
    !########################################################################
    subroutine Thomas3_Split_Initialize_Global(a, b, c, points, split)
        real(wp), intent(inout) :: a(:), b(:), c(:)
        integer(wi), intent(in) :: points(:)   ! sequence of splitting points in ascending order
        type(matrix_split_dt), intent(inout) :: split

        ! -------------------------------------------------------------------
        integer(wi) nmax, m, mmax, mmin, p, pmax, p_plus_1, k

        real(wp), allocatable :: aloc(:), bloc(:), cloc(:)
        real(wp) delta

        !########################################################################
        nmax = size(a)
        mmax = size(points)

        mmin = 1
        if (split%periodic) mmin = 0

        ! splitting data
        if (allocated(split%points)) deallocate (split%points)
        allocate (split%points(0:mmax + 1))

        split%points(1:mmax) = points(1:mmax)
        split%points(0) = nmax; split%points(mmax + 1) = nmax      ! extended
        print *, points(:)

        if (allocated(split%z)) deallocate (split%z)
        allocate (split%z(nmax, mmin:mmax))

        if (allocated(split%alpha)) deallocate (split%alpha)
        allocate (split%alpha(2, mmin:mmax))

        if (allocated(split%beta)) deallocate (split%beta)
        allocate (split%beta(mmax))

        if (allocated(split%gamma)) deallocate (split%gamma)
        allocate (split%gamma(mmax))

        ! -------------------------------------------------------------------
        ! temporary arrays to calculate z
        allocate (aloc(nmax), bloc(nmax), cloc(nmax))

        do m = mmin, mmax               ! loop over splitting points
            p = split%points(m)         ! index of splitting point
            p_plus_1 = mod(p + 1, nmax)

            ! Start definition of alpha
            split%alpha(1, m) = a(p_plus_1)
            split%alpha(2, m) = c(p)

            ! Generate matrix Ak
            b(p) = b(p) - a(p_plus_1)
            a(p_plus_1) = 0.0_wp
            b(p_plus_1) = b(p_plus_1) - c(p)
            c(p) = 0.0_wp

            ! Generate vector zk
            aloc(:) = a(:); bloc(:) = b(:); cloc(:) = c(:)
            call Thomas3_LU(nmax, aloc, bloc, cloc)

            split%z(:, m) = 0.0_wp
            split%z(p, m) = 1.0_wp
            split%z(p_plus_1, m) = 1.0_wp
            call Thomas3_Solve(nmax, 1, aloc, bloc, cloc, split%z(:, m))

            ! Complete definition of alpha
            delta = 1.0_wp + split%alpha(1, m)*split%z(p, m) + split%alpha(2, m)*split%z(p_plus_1, m)
            if (abs(delta) < small_wp) then
                call TLab_Write_ASCII(efile, __FILE__//'. Singular matrix M.')
                call TLab_Stop(DNS_ERROR_THOMAS)
            end if
            split%alpha(1, m) = -split%alpha(1, m)/delta
            split%alpha(2, m) = -split%alpha(2, m)/delta

        end do

        deallocate (aloc, bloc, cloc)

        ! Calculate beta
        do m = 1, mmax - 1
            p = split%points(m)
            split%beta(m) = split%alpha(1, m)*split%z(p, m + 1) + split%alpha(2, m)*split%z(p + 1, m + 1)
        end do

        ! Calculate gamma
        if (split%periodic) then
            do m = 1, mmax
                split%gamma(m) = split%alpha(1, 0)*split%z(nmax, m) + split%alpha(2, 0)*split%z(1, m)
            end do
        end if

        ! -------------------------------------------------------------------
        ! LU decomposition of block matrix Am
        do k = 1, mmax + 1                  ! loop over blocks
            p = mod(split%points(k - 1), nmax)
            p = p + 1
            pmax = split%points(k)
            call Thomas3_LU(pmax - p + 1, a(p:pmax), b(p:pmax), c(p:pmax))
        end do

        return
    end subroutine Thomas3_Split_Initialize_Global

    !########################################################################
    !########################################################################
    subroutine Thomas3_Split_Solve_Global(a, b, c, split, f)
        real(wp), intent(in) :: a(:), b(:), c(:)
        type(matrix_split_dt), intent(in) :: split
        real(wp), intent(inout) :: f(:, :)          ! forcing and solution

        integer(wi) n, nmax, m, mmin, mmax, p, pmax, p_plus_1, k, len
        real(wp), allocatable :: xp(:), alpha(:, :)

        !########################################################################
        len = size(f, 1)
        nmax = size(f, 2)
        mmax = size(split%points) - 2

        allocate (xp(len))
        allocate (alpha(len, 0:mmax))

        ! -------------------------------------------------------------------
        ! Solving block system Am in m+1 block
        do k = 1, mmax + 1                  ! loop over blocks
            p = mod(split%points(k - 1), nmax)
            p = p + 1
            pmax = split%points(k)
            call Thomas3_Solve(pmax - p + 1, len, &
                               a(p:pmax), b(p:pmax), c(p:pmax), f(:, p:pmax))
        end do

        ! -------------------------------------------------------------------
        ! Collect solution at block boundaries and calculate alphas
        ! written in terms of blocks to mimic parallel computing
        do k = 1, mmax + 1
            m = mod(k, mmax + 1)            ! block k calculates alpha_k
            !                                 block m+1 calculates alpha_0 for circulant case
            p = split%points(m)
            p_plus_1 = mod(p + 1, nmax)
            xp(:) = f(:, p_plus_1)          ! this info comes from block k+1
            alpha(:, m) = split%alpha(1, m)*f(:, p) + split%alpha(2, m)*xp(:)
        end do

        do k = mmax, 1, -1
            m = k                           ! block k calculates alpha k
            if (m == mmax) cycle            ! the last alpha has no beta correction
            alpha(:, m) = alpha(:, m) + split%beta(m)*alpha(:, m + 1)
        end do

        if (split%periodic) then            ! block m+1 handles circulant case
            do m = 1, mmax
                alpha(:, 0) = alpha(:, 0) + split%gamma(m)*alpha(:, m)
            end do
        end if

        ! -------------------------------------------------------------------
        ! calculate solution in m+1 blocks
        mmin = 1
        if (split%periodic) mmin = 0        ! if circulant, down to alpha_0

        do k = mmax + 1, 1, -1              ! loop over blocks
            p = mod(split%points(k - 1), nmax)
            p = p + 1
            pmax = split%points(k)
            do m = min(k, mmax), mmin, -1   ! each time I have an alpha, update
                !                             block k+1 needs same alphas as block k
                do n = p, pmax
                    f(:, n) = f(:, n) + alpha(:, m)*split%z(n, m)
                end do
            end do
        end do

        deallocate (xp, alpha)

        return
    end subroutine Thomas3_Split_Solve_Global

end module Thomas3_Split
