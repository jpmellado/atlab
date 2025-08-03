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
    public :: Thomas3_Split_Solve

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
        real(wp) :: alpha(2)
        real(wp) :: beta
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
        integer(wi) gsize, m, mmax, p, p_plus_1, k
        integer(wi) nmin, nmax, nsize

        real(wp), allocatable :: aloc(:), bloc(:), cloc(:)
        real(wp) delta

        !########################################################################
        gsize = size(a)
        mmax = size(points)

        ! splitting data
        if (allocated(split%points)) deallocate (split%points)
        allocate (split%points(0:mmax + 1))

        split%points(1:mmax) = points(1:mmax)
        split%points(mmax + 1) = gsize      ! extended for circulant case
        split%points(0) = gsize; 
        print *, points(:)

        if (allocated(split%z)) deallocate (split%z)
        allocate (split%z(gsize, mmax + 1))

        if (allocated(split%alpha)) deallocate (split%alpha)
        allocate (split%alpha(2, mmax + 1))

        if (allocated(split%beta)) deallocate (split%beta)
        allocate (split%beta(mmax))

        if (allocated(split%gamma)) deallocate (split%gamma)
        allocate (split%gamma(mmax))

        ! -------------------------------------------------------------------
        ! Calculate alpha and z

        ! temporary arrays to calculate z
        allocate (aloc(gsize), bloc(gsize), cloc(gsize))

        if (split%periodic) then        ! block mmax+1 has the circulant case
            m = mmax + 1
            p = gsize
            p_plus_1 = 1
            call Splitting(a, b, c, p, p_plus_1, split%alpha(:, m), split%z(:, m))
        end if

        do m = 1, mmax                  ! loop over splitting points
            p = split%points(m)         ! index of splitting point
            p_plus_1 = mod(p + 1, gsize)
            call Splitting(a, b, c, p, p_plus_1, split%alpha(:, m), split%z(:, m))

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
                split%gamma(m) = split%alpha(1, mmax + 1)*split%z(gsize, m) + split%alpha(2, mmax + 1)*split%z(1, m)
            end do
        end if

        ! -------------------------------------------------------------------
        ! LU decomposition of block matrix Am
        do k = 1, mmax + 1              ! loop over blocks
            nmin = mod(split%points(k - 1), gsize)
            nmin = nmin + 1
            nmax = split%points(k)
            nsize = nmax - nmin + 1
            call Thomas3_LU(nsize, a(nmin:nmax), b(nmin:nmax), c(nmin:nmax))
        end do

        return

    contains
        subroutine Splitting(a, b, c, p, p_plus_1, alpha, z)
            real(wp), intent(inout) :: a(:), b(:), c(:)
            integer(wi), intent(in) :: p, p_plus_1
            real(wp), intent(inout) :: alpha(2)
            real(wp), intent(inout) :: z(:)

            ! Start definition of alpha
            alpha(1) = a(p_plus_1)
            alpha(2) = c(p)

            ! Generate matrix Ak
            b(p) = b(p) - a(p_plus_1)
            a(p_plus_1) = 0.0_wp
            b(p_plus_1) = b(p_plus_1) - c(p)
            c(p) = 0.0_wp

            ! Generate vector zk
            aloc(:) = a(:); bloc(:) = b(:); cloc(:) = c(:)
            call Thomas3_LU(gsize, aloc, bloc, cloc)

            z(:) = 0.0_wp
            z(p) = 1.0_wp
            z(p_plus_1) = 1.0_wp
            call Thomas3_Solve(gsize, 1, aloc, bloc, cloc, z(:))

            ! Complete definition of alpha
            delta = 1.0_wp + alpha(1)*z(p) + alpha(2)*z(p_plus_1)
            if (abs(delta) < small_wp) then
                call TLab_Write_ASCII(efile, __FILE__//'. Singular matrix M.')
                call TLab_Stop(DNS_ERROR_THOMAS)
            end if
            alpha(1) = -alpha(1)/delta
            alpha(2) = -alpha(2)/delta

            return
        end subroutine Splitting

    end subroutine Thomas3_Split_Initialize_Global

    !########################################################################
    !########################################################################
    subroutine Thomas3_Split_Solve_Global(a, b, c, split, f)
        real(wp), intent(in) :: a(:), b(:), c(:)
        type(matrix_split_dt), intent(in) :: split
        real(wp), intent(inout) :: f(:, :)          ! forcing and solution

        integer(wi) n, nmin, nmax, nsize, nlines
        integer(wi) m, mmax, nblocks
        integer(wi) p, k
        integer(wi) gsize
        real(wp), allocatable :: tmp(:), xp(:, :), alpha(:, :), alpha_0(:, :)

        !########################################################################
        nlines = size(f, 1)
        mmax = size(split%points) - 2
        nblocks = mmax + 1

        gsize = size(f, 2)

        allocate (tmp(nlines))
        allocate (xp(nlines, nblocks))          ! the idea is that each block needs only one alpha
        allocate (alpha(nlines, nblocks))       ! the idea is that each block needs only one alpha
        allocate (alpha_0(nlines, nblocks))     ! the idea is that each block needs only one alpha

        ! -------------------------------------------------------------------
        ! Solving block system Am in each block
        do k = 1, nblocks                       ! loop over blocks
            nmin = mod(split%points(k - 1), gsize)
            nmin = nmin + 1
            nmax = split%points(k)
            nsize = nmax - nmin + 1
            call Thomas3_Solve(nsize, nlines, &
                               a(nmin:nmax), b(nmin:nmax), c(nmin:nmax), f(:, nmin:nmax))
        end do

        ! -------------------------------------------------------------------
        ! Reduction step
        ! pass x(:,1) to previous block
        do k = nblocks, 1, -1                   ! loop over blocks
            p = mod(split%points(k), gsize) + 1
            xp(:, k) = f(:, p)                  ! solution at left boundary of block k+1
        end do

        do k = nblocks, 1, -1                   ! loop over blocks
            m = k
            ! m = mod(k, nblocks)                 ! block k calculates alpha_k
            !                                     block m+1 calculates alpha_0 for circulant case
            nmax = split%points(k)
            alpha(:, k) = split%alpha(1, m)*f(:, nmax) + split%alpha(2, m)*xp(:, k)

            if (k <= nblocks - 2) then          ! the last 2 have no beta correction
                tmp(:) = alpha(:, m + 1)        ! this info comes from block k+1
                alpha(:, k) = alpha(:, k) + split%beta(k)*tmp(:)
            end if

            if (k <= nblocks - 1) then          ! update solution
                ! send alpha_k to all PEs with larger rank
                do m = k, nblocks
                    nmin = mod(split%points(m - 1), gsize)
                    nmin = nmin + 1
                    nmax = split%points(m)
                    do n = nmin, nmax
                        f(:, n) = f(:, n) + alpha(:, k)*split%z(n, k)
                    end do
                end do
            end if

            if (split%periodic) then         ! accumulate alpha_0
                if (k <= nblocks - 1) then
                    tmp(:) = alpha_0(:, k + 1)  ! this info comes from block k+1
                    alpha_0(:, k) = tmp(:) + split%gamma(k)*alpha(:, k)
                else
                    alpha_0(:, k) = alpha(:, k)
                end if
            end if

        end do

        if (split%periodic) then
            ! broadcast alpha_0(:) from block 1 into tmp of all ranks
            tmp(:) = alpha_0(:, 1)
            do k = nblocks, 1, -1              ! loop over blocks
                nmin = mod(split%points(k - 1), gsize)
                nmin = nmin + 1
                nmax = split%points(k)
                do n = nmin, nmax
                    f(:, n) = f(:, n) + tmp(:)*split%z(n, nblocks)
                end do
            end do
        end if

        deallocate (xp, tmp, alpha, alpha_0)

        return
    end subroutine Thomas3_Split_Solve_Global

    !########################################################################
    !########################################################################
    subroutine Thomas3_Split_Solve(split, f)
        type(thomas3_split_dt), intent(in) :: split(:)
        real(wp), intent(inout) :: f(:, :)          ! forcing and solution

        integer(wi) n, nmin, nmax, nsize, nlines
        integer(wi) m, mmax, nblocks
        integer(wi) p, k
        real(wp), allocatable :: tmp(:), xp(:, :), alpha(:, :), alpha_0(:, :)

        !########################################################################
        nlines = size(f, 1)
        nblocks = nblocks
        mmax = nblocks - 1

        allocate (tmp(nlines))
        allocate (xp(nlines, nblocks))          ! the idea is that each block needs only one alpha
        allocate (alpha(nlines, nblocks))       ! the idea is that each block needs only one alpha
        allocate (alpha_0(nlines, nblocks))     ! the idea is that each block needs only one alpha

        ! -------------------------------------------------------------------
        ! Solving block system Am in each block
        do k = 1, nblocks                       ! loop over blocks
            nmin = split(k)%pmin
            nmax = split(k)%pmax
            nsize = nmax - nmin + 1
            call Thomas3_Solve(nsize, nlines, &
                               split(k)%a(:), split(k)%b(:), split(k)%c(:), f(:, nmin:nmax))
        end do

        ! -------------------------------------------------------------------
        ! Reduction step
        ! pass x(:,1) to previous block
        do k = nblocks, 1, -1                   ! loop over blocks
            p = split(mod(k, nblocks) + 1)%pmin
            xp(:, k) = f(:, p)                  ! solution at left boundary of block k+1
        end do

        do k = nblocks, 1, -1                   ! loop over blocks
            nmax = split(k)%pmax
            alpha(:, k) = split(k)%alpha(1)*f(:, nmax) + split(k)%alpha(2)*xp(:, k)

            if (k <= nblocks - 2) then          ! the last 2 have no beta correction
                tmp(:) = alpha(:, k + 1)        ! this info comes from block k+1
                alpha(:, k) = alpha(:, k) + split(k)%beta*tmp(:)
            end if

            if (k <= nblocks - 1) then          ! update solution
                ! send alpha_k to all PEs with larger rank
                do m = k, nblocks
                    nmin = split(m)%pmin
                    nmax = split(m)%pmax
                    do n = nmin, nmax
                        f(:, n) = f(:, n) + alpha(:, k)*split(k)%z(n, k)
                    end do
                end do
            end if

            if (split(k)%periodic) then         ! accumulate alpha_0
                if (k <= nblocks - 1) then
                    tmp(:) = alpha_0(:, k + 1)  ! this info comes from block k+1
                    alpha_0(:, k) = tmp(:) + split(k)%gamma*alpha(:, k)
                else
                    alpha_0(:, k) = alpha(:, k)
                end if
            end if

        end do

        if (split(k)%periodic) then
            ! broadcast alpha_0(:) from block 1 into tmp of all ranks
            tmp(:) = alpha_0(:, 1)
            do k = nblocks, 1, -1               ! loop over blocks
                nmin = split(k)%pmin
                nmax = split(k)%pmax
                do n = nmin, nmax
                    f(:, n) = f(:, n) + tmp(:)*split(k)%z(n, 0)
                end do

            end do
        end if

        deallocate (tmp, alpha, alpha_0)

        return
    end subroutine Thomas3_Split_Solve

end module Thomas3_Split
