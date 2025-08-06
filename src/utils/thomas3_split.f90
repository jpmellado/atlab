#include "tlab_error.h"

! Matrix splitting

module Thomas3_Split
    use TLab_Constants, only: wp, wi, small_wp, efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Thomas3
    implicit none
    private

    public :: Thomas3_Split_Initialize
    public :: Thomas3_Split_Solve

    type, public :: thomas3_split_dt
        integer :: block_id
        integer(wi) :: nmin, nmax
        logical :: circulant = .true.
        real(wp), allocatable :: a(:), b(:), c(:)
        real(wp), allocatable :: z(:, :)
        real(wp) :: alpha(2)
        real(wp) :: beta
        real(wp) :: gamma
    end type thomas3_split_dt

    type, public :: data_dt                     ! for testing algorithm in serial
        real(wp), pointer :: p(:, :)
    end type data_dt

contains
    !########################################################################
    !########################################################################
    subroutine Thomas3_Split_Initialize(a, b, c, points, split)
        real(wp), intent(inout) :: a(:), b(:), c(:)
        integer(wi), intent(in) :: points(:)   ! sequence of splitting points in ascending order
        type(thomas3_split_dt), intent(inout) :: split

        ! -------------------------------------------------------------------
        integer nblocks, m
        integer(wi) nsize
        integer(wi) p, p_plus_1

        real(wp), allocatable :: aloc(:), bloc(:), cloc(:), zloc(:)
        real(wp) alpha_0(2), alpha(2)
        real(wp) delta

        !########################################################################
        nblocks = size(points)
        ! Number of coefficients are nblocks-1 for the tridiagonal case
        ! and we add one for the circulant case, which is managed by last block

        ! -------------------------------------------------------------------
        ! memory management
        nsize = split%nmax - split%nmin + 1

        if (allocated(split%a)) deallocate (split%a)
        allocate (split%a(1:nsize))
        split%a(:) = 0.0_wp
        if (allocated(split%b)) deallocate (split%b)
        allocate (split%b(1:nsize))
        split%b(:) = 0.0_wp
        if (allocated(split%c)) deallocate (split%c)
        allocate (split%c(1:nsize))
        split%c(:) = 0.0_wp

        if (allocated(split%z)) deallocate (split%z)
        allocate (split%z(1:nsize, nblocks))
        split%z(:, :) = 0.0_wp

        ! -------------------------------------------------------------------
        ! temporary arrays to calculate z_j
        nsize = size(a)
        allocate (aloc(nsize), bloc(nsize), cloc(nsize), zloc(nsize))

        if (split%circulant) then
            m = nblocks             ! last block has the circulant coefficient
            p = points(m)           ! index of splitting point
            p_plus_1 = 1

            call Splitting(a, b, c, p, p_plus_1, alpha_0, zloc)

            if (split%block_id == m) then
                split%alpha(1:2) = alpha_0(1:2)
            end if
            split%z(:, m) = zloc(split%nmin:split%nmax)
        end if

        do m = 1, nblocks - 1       ! loop over remaining coefficients
            p = points(m)           ! index of splitting point
            p_plus_1 = p + 1

            call Splitting(a, b, c, p, p_plus_1, alpha, zloc)

            if (split%block_id == m) then
                split%alpha(1:2) = alpha(1:2)
                split%gamma = alpha_0(1)*zloc(nsize) + alpha_0(2)*zloc(1)
            end if
            split%z(:, m) = zloc(split%nmin:split%nmax)

            if (split%block_id == m - 1) then
                p = points(split%block_id)
                p_plus_1 = p + 1
                split%beta = split%alpha(1)*zloc(p) + split%alpha(2)*zloc(p_plus_1)
            end if

        end do

        deallocate (aloc, bloc, cloc, zloc)

        ! -------------------------------------------------------------------
        ! block matrix Am and LU decomposition
        split%a(:) = a(split%nmin:split%nmax)
        split%b(:) = b(split%nmin:split%nmax)
        split%c(:) = c(split%nmin:split%nmax)

        nsize = split%nmax - split%nmin + 1
        call Thomas3_LU(nsize, split%a(:), split%b(:), split%c(:))

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
            call Thomas3_LU(size(aloc), aloc, bloc, cloc)

            z(:) = 0.0_wp
            z(p) = 1.0_wp
            z(p_plus_1) = 1.0_wp
            call Thomas3_Solve(size(aloc), 1, aloc, bloc, cloc, z(:))

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

    end subroutine Thomas3_Split_Initialize

    !########################################################################
    !########################################################################
    subroutine Thomas3_Split_Solve(split, f)
        type(thomas3_split_dt), intent(in) :: split(:)
        type(data_dt), intent(in) :: f(:)

        integer(wi) n, nsize, nlines
        integer(wi) m, mmax, nblocks
        integer(wi) p, k
        real(wp), allocatable :: tmp(:), xp(:, :), alpha(:, :), alpha_0(:, :)

        !########################################################################
        nlines = size(f(1)%p, 1)
        nblocks = size(split)
        mmax = nblocks - 1

        allocate (tmp(nlines))
        allocate (xp(nlines, nblocks))          ! the idea is that each block needs only one alpha
        allocate (alpha(nlines, nblocks))       ! the idea is that each block needs only one alpha
        allocate (alpha_0(nlines, nblocks))     ! the idea is that each block needs only one alpha

        ! -------------------------------------------------------------------
        ! Solving block system Am in each block
        do k = 1, nblocks                       ! loop over blocks
            nsize = split(k)%nmax - split(k)%nmin + 1
            ! nsize = size(f(k)%p, 2)
            call Thomas3_Solve(nsize, nlines, &
                               split(k)%a(:), split(k)%b(:), split(k)%c(:), f(k)%p(:, :))
        end do

        ! -------------------------------------------------------------------
        ! Reduction step
        ! pass x(:,1) to previous block
        do k = nblocks, 1, -1                   ! loop over blocks
            p = mod(k + nblocks, nblocks) + 1
            xp(:, k) = f(p)%p(:, 1)             ! solution at left boundary of block k+1
        end do

        do k = nblocks, 1, -1                   ! loop over blocks
            m = k
            nsize = split(k)%nmax - split(k)%nmin + 1
            ! nsize = size(f(k)%p, 2)
            alpha(:, k) = split(m)%alpha(1)*f(k)%p(:, nsize) + split(m)%alpha(2)*xp(:, k)

            if (k <= nblocks - 2) then          ! the last 2 have no beta correction
                tmp(:) = alpha(:, m + 1)        ! this info comes from block k+1
                alpha(:, k) = alpha(:, k) + split(k)%beta*tmp(:)
            end if

            if (k <= nblocks - 1) then          ! update solution
                ! send alpha_k to all PEs with larger rank
                do m = k, nblocks
                    nsize = split(m)%nmax - split(m)%nmin + 1
                    do n = 1, nsize
                        f(m)%p(:, n) = f(m)%p(:, n) + alpha(:, k)*split(m)%z(n, k)
                    end do
                end do
            end if

            if (split(k)%circulant) then        ! accumulate alpha_0
                if (k <= nblocks - 1) then
                    tmp(:) = alpha_0(:, k + 1)  ! this info comes from block k+1
                    alpha_0(:, k) = tmp(:) + split(k)%gamma*alpha(:, k)
                else
                    alpha_0(:, k) = alpha(:, k)
                end if
            end if

        end do

        if (split(nblocks)%circulant) then
            ! broadcast alpha_0(:) from block 1 into tmp of all ranks
            tmp(:) = alpha_0(:, 1)
            do k = nblocks, 1, -1               ! loop over blocks
                nsize = split(k)%nmax - split(k)%nmin + 1
                ! nsize = size(f(k)%p, 2)
                do n = 1, nsize
                    f(k)%p(:, n) = f(k)%p(:, n) + tmp(:)*split(k)%z(n, nblocks)
                end do
            end do
        end if

        deallocate (tmp, alpha, alpha_0)

        return
    end subroutine Thomas3_Split_Solve

end module Thomas3_Split
