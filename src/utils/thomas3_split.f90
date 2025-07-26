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

    type, public :: matrix_split_dt
        integer :: np                           ! number of splitting points
        integer(wi), allocatable :: points(:)   ! sequence of splitting points in ascending order
        real(wp), allocatable :: alpha(:, :)
        real(wp), allocatable :: z(:, :)
    end type matrix_split_dt

contains
    !########################################################################
    !########################################################################
    subroutine Thomas3_Split_Initialize(a, b, c, points, split)
        real(wp), intent(inout) :: a(:), b(:), c(:)
        integer(wi), intent(in) :: points(:)   ! sequence of splitting points in ascending order
        type(matrix_split_dt), intent(inout) :: split

        ! -------------------------------------------------------------------
        integer(wi) nmax, m, mmax, p, pmax

        real(wp), allocatable :: aloc(:), bloc(:), cloc(:)
        real(wp) delta

        !########################################################################
        nmax = size(a)
        mmax = size(points)

        ! splitting data
        if (allocated(split%points)) deallocate (split%points)
        allocate (split%points(mmax))
        split%points(:) = points(:)

        if (allocated(split%alpha)) deallocate (split%alpha)
        allocate (split%alpha(2, mmax))

        if (allocated(split%z)) deallocate (split%z)
        allocate (split%z(nmax, mmax))

        ! -------------------------------------------------------------------
        ! temporary arrays to calculate z
        allocate (aloc(nmax), bloc(nmax), cloc(nmax))

        do m = 1, mmax              ! loop over splitting points
            p = points(m)           ! index of splitting point

            ! Start definition of alpha
            split%alpha(1, m) = a(p + 1)
            split%alpha(2, m) = c(p)

            ! Generate matrix A1
            b(p) = b(p) - a(p + 1)
            a(p + 1) = 0.0_wp
            b(p + 1) = b(p + 1) - c(p)
            c(p) = 0.0_wp

            ! Generate vector z
            aloc(:) = a(:); bloc(:) = b(:); cloc(:) = c(:)
            call Thomas3_LU(nmax, aloc, bloc, cloc)

            split%z(:, m) = 0.0_wp
            split%z(p:p + 1, m) = 1.0_wp
            call Thomas3_Solve(nmax, 1, aloc, bloc, cloc, split%z(:, m))

            ! Complete definition of alpha
            delta = 1.0_wp + split%alpha(1, m)*split%z(p, m) + split%alpha(2, m)*split%z(p + 1, m)
            if (abs(delta) < small_wp) then
                call TLab_Write_ASCII(efile, __FILE__//'. Singular matrix M.')
                call TLab_Stop(DNS_ERROR_THOMAS)
            end if
            split%alpha(1, m) = -split%alpha(1, m)/delta
            split%alpha(2, m) = -split%alpha(2, m)/delta

        end do

        deallocate (aloc, bloc, cloc)

        ! -------------------------------------------------------------------
        ! LU decomposition of block matrix Am
        ! print *, a
        ! print *, b
        ! print *, c
        p = 1
        do m = 1, mmax
            pmax = points(m)
            call Thomas3_LU(pmax - p + 1, a(p:pmax), b(p:pmax), c(p:pmax))
            p = pmax + 1
        end do
        call Thomas3_LU(nmax - p + 1, a(p:nmax), b(p:nmax), c(p:nmax))

        return
    end subroutine Thomas3_Split_Initialize

    !########################################################################
    !########################################################################
    subroutine Thomas3_Split_Solve(a, b, c, split, f, wrk)
        real(wp), intent(in) :: a(:), b(:), c(:)
        type(matrix_split_dt), intent(in) :: split
        real(wp), intent(inout) :: f(:, :)          ! forcing and solution
        real(wp), intent(inout) :: wrk(:)

        integer(wi) n, nmax, m, mmax, p, pmax, len

        !########################################################################
        len = size(f, 1)
        nmax = size(f, 2)
        mmax = size(split%points)

        ! -------------------------------------------------------------------
        ! Solving block system Am
        p = 1
        do m = 1, mmax
            pmax = split%points(m)
            call Thomas3_Solve(pmax - p + 1, len, &
                               a(p:pmax), b(p:pmax), c(p:pmax), f(:, p:pmax))
            p = pmax + 1
        end do
        call Thomas3_Solve(nmax - p + 1, len, a(p:nmax), b(p:nmax), c(p:nmax), f(:, p:nmax))

        ! -------------------------------------------------------------------
        ! Obtain x
        m = mmax
        p = split%points(m)
        wrk(:) = split%alpha(1, m)*f(:, p) + split%alpha(2, m)*f(:, p + 1)
        do n = 1, nmax
            f(:, n) = f(:, n) + wrk(:)*split%z(n, m)
        end do

        ! do m = mmax - 1, 1, -1
        ! end do

        return
    end subroutine Thomas3_Split_Solve

end module Thomas3_Split
