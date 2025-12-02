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
subroutine NSE_Anelastic()
    use TLab_Constants, only: wp, wi, BCS_NN
    use TLab_Memory, only: imax, jmax, kmax, inb_flow, inb_scal
    use TLab_Arrays, only: s
    use TLab_Pointers, only: u, v, w, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7
    use DNS_Arrays
    use TimeMarching, only: dte, remove_divergence
    use Thermo_Anelastic, only: rbackground, ribackground, Thermo_Anelastic_Weight_InPlace
    use Thermo_Anelastic, only: Thermo_Anelastic_Weight_OutPlace, Thermo_Anelastic_Weight_Subtract
    use BoundaryConditions
    use OPR_Partial
    use NSE_Burgers
    use OPR_Elliptic

    implicit none

    ! -----------------------------------------------------------------------
    integer(wi) iq, is
    integer ibc
    real(wp) dummy

    ! #######################################################################
    ! Preliminaries for Scalar BC
    ! (flow BCs initialized below as they are used for pressure in between)
    ! #######################################################################
    BcsScalKmin%ref(:, :, :) = 0.0_wp ! default is no-slip (dirichlet)
    BcsScalKmax%ref(:, :, :) = 0.0_wp

    ! Keep the old tendency of the scalar at the boundary to be used in dynamic BCs
    if (any(BcsScalKmin%SfcType(1:inb_scal) == DNS_SFC_LINEAR) .or. &
        any(BcsScalKmax%SfcType(1:inb_scal) == DNS_SFC_LINEAR)) then
        do is = 1, inb_scal
            if (BcsScalKmin%SfcType(is) == DNS_SFC_LINEAR) BcsScalKmin%ref(:, :, is) = p_hs(:, :, 1, is)
            if (BcsScalKmax%SfcType(is) == DNS_SFC_LINEAR) BcsScalKmax%ref(:, :, is) = p_hs(:, :, kmax, is)
        end do
    end if

    ! #######################################################################
    ! Diffusion and advection terms
    ! #######################################################################
    ! Diagonal terms in horizontal directions and transposed velocity arrays
    call NSE_Burgers_X(0, imax, jmax, kmax, u, tmp1, tmp4)                      ! store u transposed in tmp4
    call NSE_Burgers_Y(0, imax, jmax, kmax, v, tmp2, tmp5)                      ! store v transposed in tmp5

    ! Ox momentum equation
    call NSE_Burgers_Y(0, imax, jmax, kmax, u, tmp6, tmp7, tmp5)                ! tmp5 contains v transposed
    call NSE_Burgers_Z(0, imax, jmax, kmax, u, tmp3, w)
    hq(:, 1) = hq(:, 1) + tmp1(:) + tmp6(:) + tmp3(:)

    ! Oy momentum equation
    call NSE_Burgers_X(0, imax, jmax, kmax, v, tmp6, tmp7, tmp4)                ! tmp4 contains u transposed
    call NSE_Burgers_Z(0, imax, jmax, kmax, v, tmp3, w)
    hq(:, 2) = hq(:, 2) + tmp6(:) + tmp2(:) + tmp3(:)

    ! Oz momentum equation
    call NSE_Burgers_X(0, imax, jmax, kmax, w, tmp1, tmp6, tmp4)                ! tmp4 contains u transposed
    call NSE_Burgers_Y(0, imax, jmax, kmax, w, tmp2, tmp6, tmp5)                ! tmp5 contains v transposed
    call NSE_Burgers_Z(0, imax, jmax, kmax, w, tmp3, w)
    hq(:, 3) = hq(:, 3) + tmp1(:) + tmp2(:) + tmp3(:)

    ! Scalar equations
    do is = 1, inb_scal
        call NSE_Burgers_X(is, imax, jmax, kmax, s(:, is), tmp1, tmp6, tmp4)    ! tmp4 contains u transposed
        call NSE_Burgers_Y(is, imax, jmax, kmax, s(:, is), tmp2, tmp6, tmp5)    ! tmp5 contains v transposed
        call NSE_Burgers_Z(is, imax, jmax, kmax, s(:, is), tmp3, w)
        hs(:, is) = hs(:, is) + tmp1(:) + tmp2(:) + tmp3(:)

    end do

    ! #######################################################################
    ! Pressure term
    ! #######################################################################
    ! Forcing term
    if (remove_divergence) then ! remove residual divergence
        dummy = 1.0_wp/dte
        tmp2(:) = hq(:, 1) + u(:)*dummy
        tmp3(:) = hq(:, 2) + v(:)*dummy
        tmp4(:) = hq(:, 3) + w(:)*dummy

        call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, rbackground, tmp2)
        call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, rbackground, tmp3)
        call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, rbackground, tmp4)

    else
        call Thermo_Anelastic_Weight_OutPlace(imax, jmax, kmax, rbackground, hq(:, 1), tmp2)
        call Thermo_Anelastic_Weight_OutPlace(imax, jmax, kmax, rbackground, hq(:, 2), tmp3)
        call Thermo_Anelastic_Weight_OutPlace(imax, jmax, kmax, rbackground, hq(:, 3), tmp4)

    end if

    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, tmp2, tmp1)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, tmp3, tmp2)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, tmp4, tmp3)
    tmp1(:) = tmp1(:) + tmp2(:) + tmp3(:) ! forcing term in tmp1

    ! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
    BcsFlowKmin%ref(:, :, 3) = p_hq(:, :, 1, 3)*rbackground(1)
    BcsFlowKmax%ref(:, :, 3) = p_hq(:, :, kmax, 3)*rbackground(kmax)

    ! Solution of Poisson equation: pressure in tmp1
    call OPR_Poisson(imax, jmax, kmax, BCS_NN, tmp1, tmp2, tmp3, BcsFlowKmin%ref(:, :, 3), BcsFlowKmax%ref(:, :, 3))

    ! Add pressure gradient
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, tmp1, tmp2)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, tmp1, tmp3)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, tmp1, tmp4)
    call Thermo_Anelastic_Weight_Subtract(imax, jmax, kmax, ribackground, tmp2, hq(:, 1))
    call Thermo_Anelastic_Weight_Subtract(imax, jmax, kmax, ribackground, tmp3, hq(:, 2))
    call Thermo_Anelastic_Weight_Subtract(imax, jmax, kmax, ribackground, tmp4, hq(:, 3))

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

        ! if (BcsScalKmin%type(is) /= DNS_SFC_STATIC .or. &
        !     BcsScalKmax%type(is) /= DNS_SFC_STATIC) then
        !     call BCS_SURFACE_Z(is, s, hs, tmp1, tmp2)
        ! end if

        p_hs(:, :, 1, is) = BcsScalKmin%ref(:, :, is)
        p_hs(:, :, kmax, is) = BcsScalKmax%ref(:, :, is)

    end do

    return
end subroutine NSE_Anelastic

!########################################################################
!########################################################################
subroutine NSE_Anelastic_PerVolume()
    use TLab_Constants, only: wp, wi, BCS_NN
    use TLab_Memory, only: imax, jmax, kmax, inb_flow, inb_scal
    use TLab_Arrays, only: s
    use TLab_Pointers, only: u, v, w, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7
    use TLab_Pointers_3D, only: p_q, pxy_tmp2 => tmp2, pxy_tmp3 => tmp3, pxy_tmp4 => tmp4
    use DNS_Arrays
    use TimeMarching, only: dte, remove_divergence
    use Thermo_Anelastic, only: rbackground
    use BoundaryConditions
    use OPR_Partial
    use NSE_Burgers_PerVolume
    use OPR_Elliptic, only: OPR_Poisson

    implicit none

    ! -----------------------------------------------------------------------
    integer(wi) iq, is, k
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
    ! Diagonal terms in horizontal directions and (transposed) velocity times density arrays
    call NSE_Burgers_PerVolume_X(0, imax, jmax, kmax, u, tmp1, tmp4)                   ! store rho u transposed in tmp4
    call NSE_Burgers_PerVolume_Y(0, imax, jmax, kmax, v, tmp2, tmp5)                   ! store rho v transposed in tmp5
    call NSE_Burgers_PerVolume_Z_Add(0, imax, jmax, kmax, w, hq(:, 3), tmp3, rhou_out=tmp6)          ! store rho w in tmp6

    ! Ox momentum equation
    call NSE_Burgers_PerVolume_Y(0, imax, jmax, kmax, u, tmp7, tmp3, rhou_in=tmp5)     ! tmp5 contains rho v transposed
    call NSE_Burgers_PerVolume_Z_Add(0, imax, jmax, kmax, u, hq(:, 1), tmp3, rhou_in=tmp6)
    hq(:, 1) = hq(:, 1) + tmp1(:) + tmp7(:)

    ! Oy momentum equation
    call NSE_Burgers_PerVolume_X(0, imax, jmax, kmax, v, tmp1, tmp3, rhou_in=tmp4)     ! tmp4 contains rho u transposed
    call NSE_Burgers_PerVolume_Z_Add(0, imax, jmax, kmax, v, hq(:, 2), tmp3, rhou_in=tmp6)
    hq(:, 2) = hq(:, 2) + tmp1(:) + tmp2(:)

    ! Oz momentum equation
    call NSE_Burgers_PerVolume_X(0, imax, jmax, kmax, w, tmp1, tmp3, rhou_in=tmp4)     ! tmp4 contains rho u transposed
    call NSE_Burgers_PerVolume_Y(0, imax, jmax, kmax, w, tmp2, tmp3, rhou_in=tmp5)     ! tmp5 contains rho v transposed
    hq(:, 3) = hq(:, 3) + tmp1(:) + tmp2(:)

    ! Scalar equations
    do is = 1, inb_scal
        call NSE_Burgers_PerVolume_X(is, imax, jmax, kmax, s(:, is), tmp1, tmp7, rhou_in=tmp4)    ! tmp4 contains u transposed
        call NSE_Burgers_PerVolume_Y(is, imax, jmax, kmax, s(:, is), tmp2, tmp7, rhou_in=tmp5)    ! tmp5 contains v transposed
        call NSE_Burgers_PerVolume_Z_Add(is, imax, jmax, kmax, s(:, is), hs(:, is), tmp3, rhou_in=tmp6)
        hs(:, is) = hs(:, is) + tmp1(:) + tmp2(:)

    end do

    ! #######################################################################
    ! Pressure term
    ! #######################################################################
    ! Forcing term
    if (remove_divergence) then ! remove residual divergence
        dummy = 1.0_wp/dte
        do k = 1, kmax
            pxy_tmp2(:, :, k) = p_hq(:, :, k, 1) + p_q(:, :, k, 1)*dummy*rbackground(k)
            pxy_tmp3(:, :, k) = p_hq(:, :, k, 2) + p_q(:, :, k, 2)*dummy*rbackground(k)
            pxy_tmp4(:, :, k) = p_hq(:, :, k, 3) + p_q(:, :, k, 3)*dummy*rbackground(k)
        end do

        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, tmp2, tmp1)
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, tmp3, tmp2)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, tmp4, tmp3)

    else
        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, hq(:, 1), tmp1)
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, hq(:, 2), tmp2)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, hq(:, 3), tmp3)

    end if
    tmp1(:) = tmp1(:) + tmp2(:) + tmp3(:) ! forcing term in tmp1

    ! Neumman BCs in d/dy(p) s.t. v=0 (no-penetration)
    BcsFlowKmin%ref(:, :, 3) = p_hq(:, :, 1, 3)
    BcsFlowKmax%ref(:, :, 3) = p_hq(:, :, kmax, 3)

    ! Solution of Poisson equation: pressure in tmp1
    call OPR_Poisson(imax, jmax, kmax, BCS_NN, tmp1, tmp2, tmp3, BcsFlowKmin%ref(:, :, 3), BcsFlowKmax%ref(:, :, 3))

    ! Add pressure gradient
    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, tmp1, tmp2)
    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, tmp1, tmp3)
    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, tmp1, tmp4)
    hq(:, 1) = hq(:, 1) - tmp2(:)
    hq(:, 2) = hq(:, 2) - tmp3(:)
    hq(:, 3) = hq(:, 3) - tmp4(:)

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
            call BCS_Neumann_Z_PerVolume(ibc, imax*jmax, kmax, hq(:, iq), &
                                         BcsFlowKmin%ref(:, :, iq), BcsFlowKmax%ref(:, :, iq))
        end if

        p_hq(:, :, 1, iq) = BcsFlowKmin%ref(:, :, iq)*rbackground(1)
        p_hq(:, :, kmax, iq) = BcsFlowKmax%ref(:, :, iq)*rbackground(kmax)

    end do

    do is = 1, inb_scal
        ibc = 0
        if (BcsScalKmin%type(is) == DNS_BCS_Neumann) ibc = ibc + 1
        if (BcsScalKmax%type(is) == DNS_BCS_Neumann) ibc = ibc + 2
        if (ibc > 0) then
            call BCS_Neumann_Z_PerVolume(ibc, imax*jmax, kmax, hs(:, is), &
                                         BcsScalKmin%ref(:, :, is), BcsScalKmax%ref(:, :, is))
        end if

        p_hs(:, :, 1, is) = BcsScalKmin%ref(:, :, is)*rbackground(1)
        p_hs(:, :, kmax, is) = BcsScalKmax%ref(:, :, is)*rbackground(kmax)

    end do

    return
end subroutine NSE_Anelastic_PerVolume
