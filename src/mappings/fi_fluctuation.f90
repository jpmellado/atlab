!########################################################################
! Reynolds fluctuations of array a
!########################################################################
subroutine FI_Fluctuation_InPlace(nx, ny, nz, a)
    use TLab_Constants, only: wp, wi
    use Averages, only: AVG_IK
    implicit none

    integer(wi), intent(IN) :: nx, ny, nz
    real(wp), intent(INOUT) :: a(nx, ny, nz)

    ! -------------------------------------------------------------------
    real(wp) dummy
    integer(wi) k

    ! ###################################################################
    do k = 1, nz
        dummy = AVG_IK(nx, ny, nz, k, a)
        a(:, :, k) = a(:, :, k) - dummy
    end do

    return
end subroutine FI_Fluctuation_InPlace

!########################################################################
! Reynolds fluctuations of array a
!########################################################################
subroutine FI_Fluctuation_OutPlace(nx, ny, nz, a, b)
    use TLab_Constants, only: wp, wi
    use Averages, only: AVG_IK
    implicit none

    integer(wi), intent(IN) :: nx, ny, nz
    real(wp), intent(in) :: a(nx, ny, nz)
    real(wp), intent(out) :: b(nx, ny, nz)

    ! -------------------------------------------------------------------
    real(wp) dummy
    integer(wi) k

    ! ###################################################################
    do k = 1, nz
        dummy = AVG_IK(nx, ny, nz, k, a)
        b(:, :, k) = a(:, :, k) - dummy
    end do

    return
end subroutine FI_Fluctuation_OutPlace
