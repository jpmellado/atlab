!########################################################################
!#
!# Evolution equations per unit volume, nonlinear term in convective form and the
!# viscous term explicit: 9 2nd order + 9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!#
!# It is written such that u and v transposes are calculated first for the
!# Ox and Oy momentum equations, stored in tmp4 and tmp5 and then used as needed.
!# This saves 2 transpositions.
!# Includes the scalar to benefit from the same reduction
!#
!########################################################################
subroutine NSE_Anelastic_PerVolume(hq, hs, dte, remove_divergence)
    use TLab_Constants, only: wp, wi, BCS_NN
    use TLab_Memory, only: imax, jmax, kmax, isize_field, inb_flow, inb_scal
    use TLab_Pointers, only: u, v, w, tmp1, tmp2, tmp3, tmp4
    use TLab_Arrays, only: s
    use OPR_Partial
    use NSE_Burgers
    use OPR_Elliptic, only: OPR_Poisson
    implicit none

    real(wp), intent(out) :: hq(isize_field, inb_flow)
    real(wp), intent(out) :: hs(isize_field, inb_scal)
    real(wp), intent(in) :: dte
    logical, intent(in) :: remove_divergence

    ! -----------------------------------------------------------------------
    integer(wi) is

    ! #######################################################################
    ! Diffusion and advection terms
    ! #######################################################################
    call NSE_AddBurgers_PerVolume_Z(0, imax, jmax, kmax, w, hq(:, 3), tmp1, rhou_out=tmp3)          ! store rho w in tmp3
    call NSE_AddBurgers_PerVolume_Z(0, imax, jmax, kmax, u, hq(:, 1), tmp1, rhou_in=tmp3)
    call NSE_AddBurgers_PerVolume_Z(0, imax, jmax, kmax, v, hq(:, 2), tmp1, rhou_in=tmp3)
    do is = 1, inb_scal
        call NSE_AddBurgers_PerVolume_Z(is, imax, jmax, kmax, s(:, is), hs(:, is), tmp1, rhou_in=tmp3)
    end do

    call NSE_AddBurgers_PerVolume_X(0, imax, jmax, kmax, u, hq(:, 1), tmp1, tmp3)                   ! store rho u transposed in tmp3
    call NSE_AddBurgers_PerVolume_X(0, imax, jmax, kmax, v, hq(:, 2), tmp1, tmp2, rhou_in=tmp3)
    call NSE_AddBurgers_PerVolume_X(0, imax, jmax, kmax, w, hq(:, 3), tmp1, tmp2, rhou_in=tmp3)
    do is = 1, inb_scal
        call NSE_AddBurgers_PerVolume_X(is, imax, jmax, kmax, s(:, is), hs(:, is), tmp1, tmp2, rhou_in=tmp3)
    end do

    call NSE_AddBurgers_PerVolume_Y(0, imax, jmax, kmax, v, hq(:, 2), tmp1, tmp3)                   ! store rho v transposed in tmp3
    call NSE_AddBurgers_PerVolume_Y(0, imax, jmax, kmax, u, hq(:, 1), tmp1, tmp2, rhou_in=tmp3)
    call NSE_AddBurgers_PerVolume_Y(0, imax, jmax, kmax, w, hq(:, 3), tmp1, tmp2, rhou_in=tmp3)
    do is = 1, inb_scal
        call NSE_AddBurgers_PerVolume_Y(is, imax, jmax, kmax, s(:, is), hs(:, is), tmp1, tmp2, rhou_in=tmp3)
    end do

    ! #######################################################################
    ! Pressure term
    ! #######################################################################
    ! Forcing term
    if (remove_divergence) then ! remove residual divergence
        call Add_Residual_Divergence(hq)

        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, tmp4, tmp1)
        call OPR_Partial_Y(OPR_P1_ADD, imax, jmax, kmax, tmp3, tmp4, tmp1)
        call OPR_Partial_X(OPR_P1_ADD, imax, jmax, kmax, tmp2, tmp4, tmp1) ! forcing term in tmp1

    else
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, hq(:, 3), tmp1)
        call OPR_Partial_Y(OPR_P1_ADD, imax, jmax, kmax, hq(:, 2), tmp2, tmp1)
        call OPR_Partial_X(OPR_P1_ADD, imax, jmax, kmax, hq(:, 1), tmp2, tmp1)

    end if

    ! Solution of Poisson equation: pressure in tmp1
    call OPR_Poisson(imax, jmax, kmax, BCS_NN, tmp1, tmp2, tmp3, &
                     bcs_hb=hq(1:imax*jmax, 3), &                               ! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
                     bcs_ht=hq(isize_field - imax*jmax + 1:isize_field, 3))

    ! Add pressure gradient
    call OPR_Partial_X(OPR_P1_SUBTRACT, imax, jmax, kmax, tmp1, tmp2, hq(:, 1))
    call OPR_Partial_Y(OPR_P1_SUBTRACT, imax, jmax, kmax, tmp1, tmp2, hq(:, 2))
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, tmp1, tmp2)
    hq(:, 3) = hq(:, 3) - tmp2(:)

    return

contains
    subroutine Add_Residual_Divergence(hq)
        use TLab_Pointers_3D, only: p_q, pxy_tmp2 => tmp2, pxy_tmp3 => tmp3, pxy_tmp4 => tmp4
        use Thermo_Anelastic, only: rbackground
        real(wp), intent(out) :: hq(imax, jmax, kmax, inb_flow)

        integer(wi) k
        real(wp) dummy

        dummy = 1.0_wp/dte
        do k = 1, kmax
            pxy_tmp2(:, :, k) = hq(:, :, k, 1) + p_q(:, :, k, 1)*dummy*rbackground(k)
            pxy_tmp3(:, :, k) = hq(:, :, k, 2) + p_q(:, :, k, 2)*dummy*rbackground(k)
            pxy_tmp4(:, :, k) = hq(:, :, k, 3) + p_q(:, :, k, 3)*dummy*rbackground(k)
        end do

        return
    end subroutine

end subroutine NSE_Anelastic_PerVolume
