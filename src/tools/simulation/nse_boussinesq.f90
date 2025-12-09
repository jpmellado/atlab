!########################################################################
!#
!# Evolution equations, nonlinear term in convective form and the
!# viscous term explicit: 9 2nd order + 9 1st order derivatives.
!# Pressure term requires 3 1st order derivatives
!#
!# It is written such that u and v transposes are calculated first for the
!# Ox and Oy momentum equations, stored in tmp4 and tmp5 and then used as needed.
!# This saves 2 transpositions.
!# Includes the scalar to benefit from the same reduction
!#
!########################################################################
subroutine NSE_Boussinesq()
    use TLab_Constants, only: wp, wi, BCS_NN
    use TLab_Memory, only: imax, jmax, kmax, inb_flow, inb_scal
    use TLab_Arrays, only: s
    use TLab_Pointers, only: u, v, w, tmp1, tmp2, tmp3
    use DNS_Arrays
    use TimeMarching, only: dte, remove_divergence
    use BoundaryConditions
    use OPR_Partial
    use NSE_Burgers
    use OPR_Elliptic, only: OPR_Poisson

    implicit none

    ! -----------------------------------------------------------------------
    integer(wi) iq, is
    integer ibc
    real(wp) dummy

    ! #######################################################################
    ! Preliminaries for Scalar BC
    ! (flow BCs initialized below as they are used for pressure in between)
    ! #######################################################################
    BcsScalKmin%ref = 0.0_wp ! default is no-slip (dirichlet)
    BcsScalKmax%ref = 0.0_wp

    ! #######################################################################
    ! Diffusion and advection terms
    ! #######################################################################
    call NSE_AddBurgers_PerVolume_Z(0, imax, jmax, kmax, w, hq(:, 3), tmp1, rhou_in=w)
    call NSE_AddBurgers_PerVolume_Z(0, imax, jmax, kmax, u, hq(:, 1), tmp1, rhou_in=w)
    call NSE_AddBurgers_PerVolume_Z(0, imax, jmax, kmax, v, hq(:, 2), tmp1, rhou_in=w)
    do is = 1, inb_scal
        call NSE_AddBurgers_PerVolume_Z(is, imax, jmax, kmax, s(:, is), hs(:, is), tmp1, rhou_in=w)
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
        dummy = 1.0_wp/dte
        tmp2(:) = hq(:, 3) + w(:)*dummy
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, tmp2, tmp1)
        tmp2(:) = hq(:, 2) + v(:)*dummy
        call OPR_Partial_Y(OPR_P1_ADD, imax, jmax, kmax, tmp2, tmp3, tmp1)
        tmp2(:) = hq(:, 1) + u(:)*dummy
        call OPR_Partial_X(OPR_P1_ADD, imax, jmax, kmax, tmp2, tmp3, tmp1) ! forcing term in tmp1

    else
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, hq(:, 3), tmp1)
        call OPR_Partial_Y(OPR_P1_ADD, imax, jmax, kmax, hq(:, 2), tmp2, tmp1)
        call OPR_Partial_X(OPR_P1_ADD, imax, jmax, kmax, hq(:, 1), tmp2, tmp1)

    end if

    ! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
    BcsFlowKmin%ref(:, :, 3) = p_hq(:, :, 1, 3)
    BcsFlowKmax%ref(:, :, 3) = p_hq(:, :, kmax, 3)

    ! Solution of Poisson equation: pressure in tmp1
    call OPR_Poisson(imax, jmax, kmax, BCS_NN, tmp1, tmp2, tmp3, BcsFlowKmin%ref(:, :, 3), BcsFlowKmax%ref(:, :, 3))

    ! Add pressure gradient
    call OPR_Partial_X(OPR_P1_SUBTRACT, imax, jmax, kmax, tmp1, tmp2, hq(:, 1))
    call OPR_Partial_Y(OPR_P1_SUBTRACT, imax, jmax, kmax, tmp1, tmp2, hq(:, 2))
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, tmp1, tmp2)
    hq(:, 3) = hq(:, 3) - tmp2(:)

    ! #######################################################################
    ! Boundary conditions
    ! #######################################################################
    BcsFlowKmin%ref = 0.0_wp ! default is no-slip (dirichlet)
    BcsFlowKmax%ref = 0.0_wp ! Scalar BCs initialized at start of routine

    do iq = 1, inb_flow
        ibc = 0
        if (BcsFlowKmin%type(iq) == DNS_BCS_Neumann) ibc = ibc + 1
        if (BcsFlowKmax%type(iq) == DNS_BCS_Neumann) ibc = ibc + 2
        if (ibc > 0) then
            call BCS_Neumann_Z(ibc, imax*jmax, kmax, hq(:, iq), &
                               BcsFlowKmin%ref(:, :, iq), BcsFlowKmax%ref(:, :, iq))
        end if

        p_hq(:, :, 1, iq) = BcsFlowKmin%ref(:, :, iq)
        p_hq(:, :, kmax, iq) = BcsFlowKmax%ref(:, :, iq)

    end do

    do is = 1, inb_scal
        ibc = 0
        if (BcsScalKmin%type(is) == DNS_BCS_Neumann) ibc = ibc + 1
        if (BcsScalKmax%type(is) == DNS_BCS_Neumann) ibc = ibc + 2
        if (ibc > 0) then
            call BCS_Neumann_Z(ibc, imax*jmax, kmax, hs(:, is), &
                               BcsScalKmin%ref(:, :, is), BcsScalKmax%ref(:, :, is))
        end if

        p_hs(:, :, 1, is) = BcsScalKmin%ref(:, :, is)
        p_hs(:, :, kmax, is) = BcsScalKmax%ref(:, :, is)

    end do

    return
end subroutine NSE_Boussinesq
