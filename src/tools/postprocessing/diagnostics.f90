#include "tlab_error.h"

module Diagnostics
    use TLab_Constants, only: wp, wi, small_wp
    use TLab_Constants, only: efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Pointers
    use TLab_Memory, only: imax, jmax, kmax
    use OPR_Partial
    use NavierStokes
    implicit none
    private

    public :: Diagnose_Thermodynamics
    public :: Diagnose_Buoyancy
    public :: Diagnose_PressureForce
    public :: Diagnose_PressurePartition
    public :: Diagnose_ScalarGradientEquation
    public :: Diagnose_StrainEquation
    public :: Diagnose_EnstrophyEquation
    public :: Diagnose_PotentialEnstrophy
    public :: Diagnose_Thetas_Anelastic
    public :: Diagnose_Energies_Anelastic
    public :: Diagnose_Moisture_Anelastic

contains
    !########################################################################
    !########################################################################
    subroutine Diagnose_Thermodynamics(vars)
        use TLab_Arrays, only: wrk3d, s
        use Thermodynamics
        use Thermo_Anelastic
        type(pointers_dt), allocatable, intent(out) :: vars(:)

        integer ifield, nfield

        ! #######################################################################
        select case (imode_thermo)
        case (THERMO_TYPE_ANELASTIC)
            nfield = 2
            allocate (vars(nfield))

            ! Accumulate vars
            ifield = 0

            call Thermo_Anelastic_Rho(imax, jmax, kmax, s, tmp1, wrk3d)
            ifield = ifield + 1; vars(ifield)%field => tmp1; vars(ifield)%tag = 'Density'

            call Thermo_Anelastic_T(imax, jmax, kmax, s, tmp1)
            ifield = ifield + 1; vars(ifield)%field => tmp1; vars(ifield)%tag = 'Temperature'

            ! Pressure in anelastic is simply p_background

        case (THERMO_TYPE_COMPRESSIBLE)
            nfield = 4
            allocate (vars(nfield))

            ! Accumulate vars
            ifield = 0

            ifield = ifield + 1; vars(ifield)%field => e; vars(ifield)%tag = 'InternalEnergy'
            ifield = ifield + 1; vars(ifield)%field => rho; vars(ifield)%tag = 'Density'
            ifield = ifield + 1; vars(ifield)%field => p; vars(ifield)%tag = 'Pressure'
            ifield = ifield + 1; vars(ifield)%field => T; vars(ifield)%tag = 'Temperature'

        end select

        ! Check
        if (nfield /= ifield) then
            call TLab_Write_ASCII(efile, __FILE__//'. Array space nfield incorrect.')
            call TLab_Stop(DNS_ERROR_WRKSIZE)
        end if

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine Diagnose_Buoyancy(vars)
        use TLab_Arrays, only: wrk1d, s
        use Gravity
        use Thermo_Anelastic, only: Thermo_Anelastic_AddBuoyancy
        type(pointers_dt), allocatable, intent(out) :: vars(:)

        integer ifield, nfield
        real(wp) dummy

        ! #######################################################################
        nfield = 5
        allocate (vars(nfield))

        ! Accumulate vars
        ifield = 0

        dummy = 1.0_wp/froude
        tmp1 = 0.0_wp
        select case (nse_eqns)
        case (DNS_EQNS_BOUSSINESQ)
            wrk1d(1:kmax, 1) = bbackground(1:kmax)
            bbackground(1:kmax) = 0.0_wp
            call Gravity_AddSource(gravityProps, imax, jmax, kmax, s, tmp1, dummy)
            bbackground(1:kmax) = wrk1d(1:kmax, 1)

        case (DNS_EQNS_ANELASTIC)
            call Thermo_Anelastic_AddBuoyancy(imax, jmax, kmax, s, tmp1, dummy)

        end select
        ifield = ifield + 1; vars(ifield)%field => tmp1; vars(ifield)%tag = 'Buoyancy'

        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, tmp1, tmp2)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, u, tmp3)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, v, tmp4)
        tmp2(:) = abs(tmp2(:))/(tmp3(:)**2 + tmp4(:)**2 + small_wp)
        ifield = ifield + 1; vars(ifield)%field => tmp2; vars(ifield)%tag = 'GradientRi'

        ! buoyancy flux along Oz
        tmp3(:) = tmp1(:)*w(:)
        ifield = ifield + 1; vars(ifield)%field => tmp3; vars(ifield)%tag = 'Fwb'

        ! buoyancy fluctuation
        call FI_Fluctuation_OutPlace(imax, jmax, kmax, tmp1, tmp4)
        ifield = ifield + 1; vars(ifield)%field => tmp4; vars(ifield)%tag = 'bPrime'

        ! covariance; turbulent buoyancy flux
        call FI_Fluctuation_OutPlace(imax, jmax, kmax, w, tmp5)
        tmp5(:) = tmp5(:)*tmp4(:)
        ifield = ifield + 1; vars(ifield)%field => tmp5; vars(ifield)%tag = 'Cwb'

        ! Using buoyancy to calculate density
        ! wrk1d(1:kmax, 1) = bbackground(1:kmax)
        ! bbackground(1:kmax) = 0.0_wp
        ! dummy = 1.0_wp/froude
        ! txc(1:isize_field, 1) = 1.0_wp
        ! call Gravity_AddSource(gravityProps, imax, jmax, kmax, s, txc(:, 1), dummy)
        ! ! txc(1:isize_field, 1) = txc(1:isize_field, 1)/froude + 1.0_wp
        ! bbackground(1:kmax) = wrk1d(1:kmax, 1)
        ! call Write_Visuals(plot_file, txc(:, 1:1))

        ! Check
        if (nfield /= ifield) then
            call TLab_Write_ASCII(efile, __FILE__//'. Array space nfield incorrect.')
            call TLab_Stop(DNS_ERROR_WRKSIZE)
        end if

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine Diagnose_PressureForce(vars)
        use TLab_Arrays, only: q, s, txc
        use NSE_Pressure
        type(pointers_dt), allocatable, intent(out) :: vars(:)

        integer ifield, nfield

        ! #######################################################################
        nfield = 5
        allocate (vars(nfield))

        ! Accumulate vars
        ifield = 0

        select case (nse_eqns)
        case (DNS_EQNS_COMPRESSIBLE)
            tmp1 = q(:, 6)

        case (DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC)
            call NSE_Pressure_Incompressible(q, s, txc(:, 1), txc(:, 2), txc(:, 5), txc(:, 6))

        end select
        ifield = ifield + 1; vars(ifield)%field => tmp1; vars(ifield)%tag = 'Pressure'

        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, tmp1, tmp2)
        tmp2(:) = -tmp2(:)*u(:)
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, tmp1, tmp3)
        tmp2(:) = tmp2(:) - tmp3(:)*v(:)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, tmp1, tmp3)
        tmp2(:) = tmp2(:) - tmp3(:)*w(:)
        ifield = ifield + 1; vars(ifield)%field => tmp2; vars(ifield)%tag = 'PressureGradientPower'

        call FI_Fluctuation_OutPlace(imax, jmax, kmax, tmp1, tmp7)      ! p fluctuation

        call FI_Fluctuation_OutPlace(imax, jmax, kmax, u, tmp4)         ! u fluctuation
        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, tmp4, tmp3)
        tmp3(:) = tmp3(:)*tmp7(:)
        ifield = ifield + 1; vars(ifield)%field => tmp3; vars(ifield)%tag = 'PressureStrainX'

        call FI_Fluctuation_OutPlace(imax, jmax, kmax, v, tmp5)         ! v fluctuation
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, tmp5, tmp4)
        tmp4(:) = tmp4(:)*tmp7(:)
        ifield = ifield + 1; vars(ifield)%field => tmp4; vars(ifield)%tag = 'PressureStrainY'

        call FI_Fluctuation_OutPlace(imax, jmax, kmax, w, tmp6)         ! w fluctuation
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, tmp6, tmp5)
        tmp5(:) = tmp5(:)*tmp7(:)
        ifield = ifield + 1; vars(ifield)%field => tmp5; vars(ifield)%tag = 'PressureStrainZ'

        ! Check
        if (nfield /= ifield) then
            call TLab_Write_ASCII(efile, __FILE__//'. Array space nfield incorrect.')
            call TLab_Stop(DNS_ERROR_WRKSIZE)
        end if

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine Diagnose_PressurePartition(vars)
        use TLab_Arrays, only: q, s, txc
        use NSE_Pressure
        type(pointers_dt), allocatable, intent(out) :: vars(:)

        integer ifield, nfield

        ! #######################################################################
        nfield = 2
        allocate (vars(nfield))

        ! Accumulate vars
        ifield = 0

        txc(:, 2:4) = 0.0_wp
        call NSE_Pressure_Incompressible(txc(:, 2), s, txc(:, 1), txc(:, 5), txc(:, 8), txc(:, 9))
        ifield = ifield + 1; vars(ifield)%field => tmp1; vars(ifield)%tag = 'PressureHydrostatic'

        select case (nse_eqns)
        case (DNS_EQNS_COMPRESSIBLE)
            tmp2 = q(:, 6)

        case (DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC)
            call NSE_Pressure_Incompressible(q, s, txc(:, 2), txc(:, 3), txc(:, 6), txc(:, 7))

        end select

        tmp2(:) = tmp2(:) - tmp1(:)
        ifield = ifield + 1; vars(ifield)%field => tmp2; vars(ifield)%tag = 'PressureHydrodynamic'

        ! Check
        if (nfield /= ifield) then
            call TLab_Write_ASCII(efile, __FILE__//'. Array space nfield incorrect.')
            call TLab_Stop(DNS_ERROR_WRKSIZE)
        end if

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine Diagnose_ScalarGradientEquation(is, vars)
        use TLab_Arrays, only: s
        use FI_VECTORCALCULUS
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
        use FI_VORTICITY_EQN
        use FI_VECTORCALCULUS
        use Gravity
        type(pointers_dt), allocatable, intent(out) :: vars(:)

        integer ifield, nfield

        ! #######################################################################
        nfield = 7
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
    subroutine Diagnose_PotentialEnstrophy(vars)
        use TLab_Arrays, only: s
        use Gravity
        use Thermo_Anelastic
        use FI_VORTICITY_EQN
        use FI_VECTORCALCULUS
        type(pointers_dt), allocatable, intent(out) :: vars(:)

        integer ifield, nfield
        real(wp) dummy

        ! #######################################################################
        nfield = 2
        allocate (vars(nfield))

        ! Accumulate vars
        ifield = 0

        call FI_VORTICITY(imax, jmax, kmax, u, v, w, tmp1, tmp2, tmp3)
        tmp1(:) = log10(tmp1(:) + small_wp)
        ifield = ifield + 1; vars(ifield)%field => tmp1; vars(ifield)%tag = 'LogEnstrophy'

        call FI_CURL(imax, jmax, kmax, u, v, w, tmp2, tmp3, tmp4, tmp5) ! vorticity
        dummy = 1.0_wp/froude                                           ! grad b
        tmp5(:) = 0.0_wp
        select case (nse_eqns)
        case (DNS_EQNS_BOUSSINESQ)
            ! wrk1d(1:kmax, 1) = bbackground(1:kmax)
            ! bbackground(1:kmax) = 0.0_wp
            call Gravity_AddSource(gravityProps, imax, jmax, kmax, s, tmp5, dummy)
            ! bbackground(1:kmax) = wrk1d(1:kmax, 1)

        case (DNS_EQNS_ANELASTIC)
            call Thermo_Anelastic_AddBuoyancy(imax, jmax, kmax, s, tmp5, dummy)

        end select
        call OPR_Partial_X(OPR_P1, imax, jmax, kmax, tmp5, tmp6)
        tmp2(:) = tmp2(:)*tmp6(:)
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, tmp5, tmp6)
        tmp2(:) = tmp2(:) + tmp3(:)*tmp6(:)
        call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, tmp5, tmp6)
        tmp2(:) = tmp2(:) + tmp4(:)*tmp6(:)
        tmp2(:) = log10(tmp2(:) + small_wp)
        ifield = ifield + 1; vars(ifield)%field => tmp2; vars(ifield)%tag = 'LogPotentialEnstrophy'

        ! To be checked
        ! call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
        ! txc(1:isize_field, 6) = txc(1:isize_field, 1)*txc(1:isize_field, 1) &
        !                         + txc(1:isize_field, 2)*txc(1:isize_field, 2) &
        !                         + txc(1:isize_field, 3)*txc(1:isize_field, 3) ! Enstrophy
        ! call OPR_Partial_X(OPR_P1, imax, jmax, kmax, s(1, 1), txc(1, 4))
        ! txc(1:isize_field, 1) = txc(1:isize_field, 1)*txc(1:isize_field, 4)
        ! txc(1:isize_field, 5) = txc(1:isize_field, 4)*txc(1:isize_field, 4) ! norm grad b
        ! call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, s(1, 1), txc(1, 4))
        ! txc(1:isize_field, 1) = txc(1:isize_field, 1) + txc(1:isize_field, 2)*txc(1:isize_field, 4)
        ! txc(1:isize_field, 5) = txc(1:isize_field, 5) + txc(1:isize_field, 4)*txc(1:isize_field, 4) ! norm grad b
        ! call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, s(1, 1), txc(1, 4))
        ! txc(1:isize_field, 1) = txc(1:isize_field, 1) + txc(1:isize_field, 3)*txc(1:isize_field, 4)
        ! txc(1:isize_field, 5) = txc(1:isize_field, 5) + txc(1:isize_field, 4)*txc(1:isize_field, 4) ! norm grad b

        ! txc(1:isize_field, 5) = sqrt(txc(1:isize_field, 5) + small_wp)
        ! txc(1:isize_field, 6) = sqrt(txc(1:isize_field, 6) + small_wp)
        ! txc(1:isize_field, 2) = txc(1:isize_field, 1)/(txc(1:isize_field, 5)*txc(1:isize_field, 6)) ! Cosine of angle between 2 vectors

        ! ifield = ifield + 1; vars_old(ifield)%field => txc(:, 1); vars_old(ifield)%tag = 'PV'
        ! ifield = ifield + 1; vars_old(ifield)%field => txc(:, 2); vars_old(ifield)%tag = 'Cos'

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
        ifield = ifield + 1; vars(ifield)%field => tmp1; vars(ifield)%tag = 'Pressure2SijPij'

        call FI_STRAIN_PRODUCTION(imax, jmax, kmax, u, v, w, &
                                  tmp2, tmp3, tmp4, tmp5, tmp6, tmp7)
        tmp2(:) = 2.0_wp*tmp2(:)
        ifield = ifield + 1; vars(ifield)%field => tmp2; vars(ifield)%tag = 'ProductionMs2SijSjkS_ki'

        call FI_STRAIN_DIFFUSION(imax, jmax, kmax, u, v, w, &
                                 tmp3, tmp4, tmp5, tmp6, tmp7, tmp8)
        tmp3(:) = 2.0_wp*visc*tmp3(:)
        ifield = ifield + 1; vars(ifield)%field => tmp3; vars(ifield)%tag = 'DiffusionNuSijLapSij'

        call FI_STRAIN(imax, jmax, kmax, u, v, w, tmp4, tmp5, tmp6)
        tmp4(:) = 2.0_wp*tmp4(:)
        tmp5(:) = log(tmp4(:))
        ifield = ifield + 1; vars(ifield)%field => tmp4; vars(ifield)%tag = '2SijSij'
        ifield = ifield + 1; vars(ifield)%field => tmp5; vars(ifield)%tag = 'Ln2SijSij'

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
        nfield = 2
        allocate (vars(nfield))

        ! Accumulate vars
        ifield = 0

        call Thermo_Anelastic_Theta(imax, jmax, kmax, s, tmp1, wrk3d)
        ifield = ifield + 1; vars(ifield)%field => tmp1; vars(ifield)%tag = 'Theta'

        call Thermo_Anelastic_ThetaV(imax, jmax, kmax, s, tmp2, wrk3d)
        ifield = ifield + 1; vars(ifield)%field => tmp2; vars(ifield)%tag = 'ThetaV'

        ! Check
        if (nfield /= ifield) then
            call TLab_Write_ASCII(efile, __FILE__//'. Array space nfield incorrect.')
            call TLab_Stop(DNS_ERROR_WRKSIZE)
        end if

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine Diagnose_Energies_Anelastic(vars)
        use TLab_Arrays, only: wrk3d, s
        use Thermo_Anelastic
        type(pointers_dt), allocatable, intent(out) :: vars(:)

        integer ifield, nfield

        ! #######################################################################
        nfield = 3
        allocate (vars(nfield))

        ! Accumulate vars
        ifield = 0

        call Thermo_Anelastic_ThetaE(imax, jmax, kmax, s, tmp1, wrk3d)
        ifield = ifield + 1; vars(ifield)%field => tmp1; vars(ifield)%tag = 'ThetaE'

        call Thermo_Anelastic_ThetaL(imax, jmax, kmax, s, tmp2, wrk3d)
        ifield = ifield + 1; vars(ifield)%field => tmp2; vars(ifield)%tag = 'ThetaL'

        call Thermo_Anelastic_MSE(imax, jmax, kmax, s, tmp3)
        ifield = ifield + 1; vars(ifield)%field => tmp3; vars(ifield)%tag = 'MoistStaticEnergy'

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
        use Thermo_Base, only: Thermo_Psat_Polynomial
        use Thermo_Anelastic
        use Thermo_AirWater, only: inb_scal_T
        type(pointers_dt), allocatable, intent(out) :: vars(:)

        integer ifield, nfield

        ! #######################################################################
        nfield = 2
        allocate (vars(nfield))

        ! Accumulate vars
        ifield = 0

        call Thermo_Anelastic_RH(imax, jmax, kmax, s, tmp1, wrk3d)
        ifield = ifield + 1; vars(ifield)%field => tmp1; vars(ifield)%tag = 'RelativeHumidity'

        call Thermo_Psat_Polynomial(imax*jmax*kmax, s(:, inb_scal_T), tmp2)
        ifield = ifield + 1; vars(ifield)%field => tmp2; vars(ifield)%tag = 'SaturationPressure'

        ! Saturation specific humidity?

        ! Check
        if (nfield /= ifield) then
            call TLab_Write_ASCII(efile, __FILE__//'. Array space nfield incorrect.')
            call TLab_Stop(DNS_ERROR_WRKSIZE)
        end if

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine Diagnose_Template(vars)
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

end module Diagnostics
