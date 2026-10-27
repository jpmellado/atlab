module TLab_Pointers
    use TLab_Constants, only: wp
    implicit none

    type :: pointers_dt
        sequence
        character(len=32) :: tag
        real(wp), pointer :: field(:)
    end type pointers_dt

    real(wp), pointer :: u(:) => null()
    real(wp), pointer :: v(:) => null()
    real(wp), pointer :: w(:) => null()
    real(wp), pointer :: e(:) => null()
    real(wp), pointer :: rho(:) => null()
    real(wp), pointer :: p(:) => null()
    real(wp), pointer :: T(:) => null()
    real(wp), pointer :: vis(:) => null()

    real(wp), pointer :: tmp1(:) => null()
    real(wp), pointer :: tmp2(:) => null()
    real(wp), pointer :: tmp3(:) => null()
    real(wp), pointer :: tmp4(:) => null()
    real(wp), pointer :: tmp5(:) => null()
    real(wp), pointer :: tmp6(:) => null()
    real(wp), pointer :: tmp7(:) => null()
    real(wp), pointer :: tmp8(:) => null()
    real(wp), pointer :: tmp9(:) => null()

end module TLab_Pointers

! ###################################################################
! ###################################################################
module TLab_Pointers_3D
    use TLab_Constants, only: wp
    implicit none

    type :: pointers3d_dt
        sequence
        character(len=32) :: tag
        real(wp), pointer :: field(:, :, :)
    end type pointers3d_dt

    real(wp), pointer :: u(:, :, :) => null()
    real(wp), pointer :: v(:, :, :) => null()
    real(wp), pointer :: w(:, :, :) => null()

    real(wp), pointer :: p_q(:, :, :, :) => null()
    real(wp), pointer :: p_s(:, :, :, :) => null()
    real(wp), pointer :: p_wrk1d(:, :) => null()
    real(wp), pointer :: p_wrk2d(:, :, :) => null()
    real(wp), pointer :: p_wrk3d(:, :, :) => null()

    real(wp), pointer :: tmp1(:, :, :) => null()
    real(wp), pointer :: tmp2(:, :, :) => null()
    real(wp), pointer :: tmp3(:, :, :) => null()
    real(wp), pointer :: tmp4(:, :, :) => null()
    real(wp), pointer :: tmp5(:, :, :) => null()
    real(wp), pointer :: tmp6(:, :, :) => null()
    real(wp), pointer :: tmp7(:, :, :) => null()
    real(wp), pointer :: tmp8(:, :, :) => null()
    real(wp), pointer :: tmp9(:, :, :) => null()

end module TLab_Pointers_3D

! ###################################################################
! ###################################################################
module TLab_Pointers_2D
    use TLab_Constants, only: wp
    implicit none

    real(wp), pointer :: pxy_q(:, :, :) => null()
    real(wp), pointer :: pxy_s(:, :, :) => null()

    real(wp), pointer :: pxy_wrk3d(:, :) => null()
    real(wp), pointer :: pxz_wrk3d(:, :) => null()
    real(wp), pointer :: pyz_wrk3d(:, :) => null()

    real(wp), pointer :: pxy_wrk2d(:, :) => null()
    real(wp), pointer :: pxz_wrk2d(:, :) => null()
    real(wp), pointer :: pyz_wrk2d(:, :) => null()

end module TLab_Pointers_2D

! ###################################################################
! ###################################################################
module TLab_Pointers_C
    use TLab_Constants, only: wp
    implicit none

    complex(wp), pointer :: c_wrk3d(:) => null()

end module TLab_Pointers_C

