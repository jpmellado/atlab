#include "tlab_error.h"

module Diagnostics
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Pointers
    use TLab_Memory, only: imax, jmax, kmax
    use NavierStokes
    implicit none
    private

    public :: Diagnose_ScalarGradientEquation
    public :: Diagnose_StrainEquation
    public :: Diagnose_EnstrophyEquation
    public :: Diagnose_Thetas_Anelastic
    public :: Diagnose_Moisture_Anelastic

contains
    !########################################################################
    !########################################################################
    subroutine Diagnose_ScalarGradientEquation(is, vars)
        use TLab_Arrays, only: s
        use FI_GRADIENT_EQN
        integer, intent(in) :: is
        type(pointers_dt), allocatable, intent(out) :: vars(:)

        integer ifield, nfield
        real(wp) diff

        ! #######################################################################
        nfield = 5
        allocate (vars(nfield))

        ! Accumulate vars
        ifield = 0

        call FI_GRADIENT_PRODUCTION(imax, jmax, kmax, s(:, is), u, v, w, &
                                    tmp1, tmp2, tmp3, tmp4, tmp5, tmp6)
        ifield = ifield + 1; vars(ifield)%field => tmp1; vars(ifield)%tag = 'ProductionMsGiGjSij'

        if (nse_diffusion == EQNS_NONE) then; diff = 0.0_wp
        else; diff = visc/schmidt(is)
        end if
        call FI_GRADIENT_DIFFUSION(imax, jmax, kmax, s(:, is), &
                                   tmp2, tmp3, tmp4, tmp5, tmp6, tmp7)

        tmp2(:) = tmp2(:)*diff
        ifield = ifield + 1; vars(ifield)%field => tmp2; vars(ifield)%tag = 'DiffusionNuGiLapGi'

        call FI_GRADIENT(imax, jmax, kmax, s(:, is), tmp3, tmp4)
        ifield = ifield + 1; vars(ifield)%field => tmp3; vars(ifield)%tag = 'GiGi'

        tmp4(:) = log(tmp3(:))
        ifield = ifield + 1; vars(ifield)%field => tmp4; vars(ifield)%tag = 'LnGiGi'; 
        tmp5(:) = tmp1(:)/tmp3(:)
        ifield = ifield + 1; vars(ifield)%field => tmp5; vars(ifield)%tag = 'StrainAMsNiNjSij'

        ! Check
        if (nfield /= ifield) then
            call TLab_Write_ASCII(efile, __FILE__//'. Array space nfield incorrect.')
            call TLab_Stop(DNS_ERROR_WRKSIZE)
        end if

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine Diagnose_EnstrophyEquation(vars)
        use TLab_Memory, only: isize_field
        use TLab_Arrays, only: txc, s, wrk3d, wrk1d
        use OPR_Partial
        use FI_VORTICITY_EQN
        use FI_VECTORCALCULUS
        use Gravity
        type(pointers_dt), allocatable, intent(out) :: vars(:)

        integer ifield, nfield

        ! #######################################################################
        nfield = 8
        allocate (vars(nfield))

        ! Accumulate vars
        ifield = 0

        if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
            ! txc(:, 4) = 0.0_wp; txc(:, 5) = 0.0_wp; txc(:, 6) = 0.0_wp

            select case (nse_eqns)
            case (DNS_EQNS_ANELASTIC)
                ! call Thermo_Anelastic_BUOYANCY(imax, jmax, kmax, s, wrk3d)

            case (DNS_EQNS_BOUSSINESQ)
                wrk1d(1:kmax, 1) = bbackground(1:kmax)
                bbackground(1:kmax) = 0.0_wp
                wrk3d(1:isize_field) = 0.0_wp
                call Gravity_AddSource(gravityProps, imax, jmax, kmax, s, wrk3d, gravityProps%vector(3))
                bbackground(1:kmax) = wrk1d(1:kmax, 1)
            end select
            s(1:isize_field, 1) = wrk3d(1:isize_field)

            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, s, txc(1, 4))
            txc(:, 4) = -txc(:, 4)
            txc(:, 5) = 0.0_wp
            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, s, txc(1, 6))

        else
            call FI_VORTICITY_BAROCLINIC(imax, jmax, kmax, rho, p, txc(1, 4), txc(1, 3), txc(1, 7))
        end if
        ! result vector in txc1, txc2, txc3
        call FI_CURL(imax, jmax, kmax, u, v, w, txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 7))
        ! scalar product, store in txc8
        txc(1:isize_field, 8) = txc(1:isize_field, 1)*txc(1:isize_field, 4) &
                                + txc(1:isize_field, 2)*txc(1:isize_field, 5) + txc(1:isize_field, 3)*txc(1:isize_field, 6)

        call FI_VORTICITY_PRODUCTION(imax, jmax, kmax, u, v, w, txc(1, 1), &
                                     txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))

        call FI_VORTICITY_DIFFUSION(imax, jmax, kmax, u, v, w, txc(1, 2), &
                                    txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7))
        txc(1:isize_field, 2) = visc*txc(1:isize_field, 2)

        call FI_VORTICITY(imax, jmax, kmax, u, v, w, txc(1, 3), txc(1, 4), txc(1, 5))

        call FI_INVARIANT_P(imax, jmax, kmax, u, v, w, txc(1, 4), txc(1, 5))

        txc(1:isize_field, 5) = txc(1:isize_field, 4)*txc(1:isize_field, 3) ! -w^2 div(u)
        txc(1:isize_field, 4) = txc(1:isize_field, 1)/txc(1:isize_field, 3) ! production rate
        txc(1:isize_field, 6) = log(txc(1:isize_field, 3))                  ! ln(w^2)

        ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'WiWi'
        ifield = ifield + 1; vars(ifield)%field => txc(:, 6); vars(ifield)%tag = 'LnWiWi'
        ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'ProductionWiWjSij'
        ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'DiffusionNuWiLapWi'
        ifield = ifield + 1; vars(ifield)%field => txc(:, 5); vars(ifield)%tag = 'DilatationMsWiWiDivU'
        ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'Baroclinic'
        ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'RateANiNjSij'

        ! Check
        if (nfield /= ifield) then
            call TLab_Write_ASCII(efile, __FILE__//'. Array space nfield incorrect.')
            call TLab_Stop(DNS_ERROR_WRKSIZE)
        end if

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine Diagnose_StrainEquation(vars)
        use TLab_Arrays, only: txc, q, s
        use FI_STRAIN_EQN
        use NSE_Pressure
        type(pointers_dt), allocatable, intent(out) :: vars(:)

        integer ifield, nfield

        ! #######################################################################
        nfield = 5
        allocate (vars(nfield))

        ! Accumulate vars
        ifield = 0

        if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
            call NSE_Pressure_Incompressible(q, s, txc(:, 1), txc(:, 2), txc(:, 5), txc(:, 6))
            call FI_STRAIN_PRESSURE(imax, jmax, kmax, u, v, w, tmp1, &
                                    tmp2, tmp3, tmp4, tmp5, tmp6)
        else
            call FI_STRAIN_PRESSURE(imax, jmax, kmax, u, v, w, p, &
                                    tmp2, tmp3, tmp4, tmp5, tmp6)
        end if
        tmp1(:) = 2.0_wp*tmp2(:)
        ifield = ifield + 1; vars(5)%field => tmp1; vars(ifield)%tag = 'Pressure2SijPij'

        call FI_STRAIN_PRODUCTION(imax, jmax, kmax, u, v, w, &
                                  tmp2, tmp3, tmp4, tmp5, tmp6, tmp7)
        tmp2(:) = 2.0_wp*tmp2(:)
        ifield = ifield + 1; vars(3)%field => tmp2; vars(ifield)%tag = 'ProductionMs2SijSjkS_ki'

        call FI_STRAIN_DIFFUSION(imax, jmax, kmax, u, v, w, &
                                 tmp3, tmp4, tmp5, tmp6, tmp7, tmp8)
        tmp3(:) = 2.0_wp*visc*tmp3(:)
        ifield = ifield + 1; vars(4)%field => tmp3; vars(ifield)%tag = 'DiffusionNuSijLapSij'

        call FI_STRAIN(imax, jmax, kmax, u, v, w, tmp4, tmp5, tmp6)
        tmp4(:) = 2.0_wp*tmp4(:)
        tmp5(:) = log(tmp4(:))
        ifield = ifield + 1; vars(1)%field => tmp4; vars(ifield)%tag = '2SijSij'
        ifield = ifield + 1; vars(2)%field => tmp5; vars(ifield)%tag = 'Ln2SijSij'

        ! Check
        if (nfield /= ifield) then
            call TLab_Write_ASCII(efile, __FILE__//'. Array space nfield incorrect.')
            call TLab_Stop(DNS_ERROR_WRKSIZE)
        end if

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine Diagnose_Fluxes(vars)
        type(pointers_dt), allocatable, intent(out) :: vars(:)

        integer ifield, nfield

        ! #######################################################################
        nfield = 4
        allocate (vars(nfield))

        ! Accumulate vars
        ifield = 0

        ! Check
        if (nfield /= ifield) then
            call TLab_Write_ASCII(efile, __FILE__//'. Array space nfield incorrect.')
            call TLab_Stop(DNS_ERROR_WRKSIZE)
        end if

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine Diagnose_Thetas_Anelastic(vars)
        use TLab_Arrays, only: wrk3d, s
        use Thermo_Anelastic
        type(pointers_dt), allocatable, intent(out) :: vars(:)

        integer ifield, nfield

        ! #######################################################################
        nfield = 4
        allocate (vars(nfield))

        ! Accumulate vars
        ifield = 0

        call Thermo_Anelastic_Theta(imax, jmax, kmax, s, tmp2, wrk3d)
        ifield = ifield + 1; vars(1)%field => tmp2; vars(ifield)%tag = 'Theta'

        call Thermo_Anelastic_ThetaV(imax, jmax, kmax, s, tmp3, wrk3d)
        ifield = ifield + 1; vars(1)%field => tmp3; vars(ifield)%tag = 'ThetaV'

        call Thermo_Anelastic_ThetaE(imax, jmax, kmax, s, tmp4, wrk3d)
        ifield = ifield + 1; vars(1)%field => tmp4; vars(ifield)%tag = 'ThetaE'

        call Thermo_Anelastic_ThetaL(imax, jmax, kmax, s, tmp5, wrk3d)
        ifield = ifield + 1; vars(1)%field => tmp5; vars(ifield)%tag = 'ThetaL'

        ! Check
        if (nfield /= ifield) then
            call TLab_Write_ASCII(efile, __FILE__//'. Array space nfield incorrect.')
            call TLab_Stop(DNS_ERROR_WRKSIZE)
        end if

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine Diagnose_Moisture_Anelastic(vars)
        use TLab_Arrays, only: wrk3d, s
        use Thermo_Anelastic
        type(pointers_dt), allocatable, intent(out) :: vars(:)

        integer ifield, nfield

        ! #######################################################################
        nfield = 4
        allocate (vars(nfield))

        ! Accumulate vars
        ifield = 0

        call Thermo_Anelastic_RH(imax, jmax, kmax, s, tmp1, wrk3d)
        ifield = ifield + 1; vars(1)%field => tmp1; vars(ifield)%tag = 'RelativeHumidity'

        ! call Thermo_Anelastic_Weight_DewPoint(imax, jmax, kmax, s, tmp2, wrk3d)
        ! ifield = ifield + 1; vars(1)%field => tmp2; vars(ifield)%tag = 'Dewpoint'

        ! Saturation pressure

        ! Saturation specific humidity

        ! Check
        if (nfield /= ifield) then
            call TLab_Write_ASCII(efile, __FILE__//'. Array space nfield incorrect.')
            call TLab_Stop(DNS_ERROR_WRKSIZE)
        end if

        return
    end subroutine

end module Diagnostics
