module TLab_CacheBlock
    use TLab_Constants, only: wp, wi
    implicit none
    private

    public :: tlab_cache_reduce_z
    public :: tlab_cache_spread_add_z
    public :: tlab_cache_reduce_y
    public :: tlab_cache_spread_add_y

contains
    !########################################################################
    !########################################################################
    subroutine tlab_cache_reduce_z(a, nlines, nca, nmax, b)
        integer(wi), intent(in) :: nlines, nca, nmax
        real(wp), intent(in) :: a(nca, *)
        real(wp), intent(out) :: b(nlines, nmax)

        integer n

        do n = 1, nmax
            b(1:nlines, n) = a(1:nlines, n)
        end do

        return
    end subroutine

    subroutine tlab_cache_spread_add_z(a, nlines, ncb, nmax, b)
        integer(wi), intent(in) :: nlines, ncb, nmax
        real(wp), intent(in) :: a(nlines, nmax)
        real(wp), intent(out) :: b(ncb, *)

        integer n

        do n = 1, nmax
            b(1:nlines, n) = b(1:nlines, n) + a(1:nlines, n)
        end do

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine tlab_cache_reduce_y(a, nx, ny, nz, i, k, nlines, b)
        real(wp), intent(in) :: a(nx, ny, nz)
        integer, intent(in) :: nx, ny, nz, i, k, nlines
        real(wp), intent(out) :: b(*)

        integer ip, jj, kk

        ip = 0
        do jj = 1, ny
            do kk = k, k + nlines - 1
                ip = ip + 1
                b(ip) = a(i, jj, kk)
            end do
        end do

        return
    end subroutine

    subroutine tlab_cache_spread_add_y(a, nx, ny, nz, i, k, nlines, b)
        real(wp), intent(in) :: a(*)
        integer, intent(in) :: nx, ny, nz, i, k, nlines
        real(wp), intent(inout) :: b(nx, ny, nz)

        integer ip, jj, kk

        ip = 0
        do jj = 1, ny
            do kk = k, k + nlines - 1
                ip = ip + 1
                b(i, jj, kk) = b(i, jj, kk) + a(ip)
            end do
        end do

        return
    end subroutine

end module
