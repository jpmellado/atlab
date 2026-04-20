#include "tlab_error.h"

module Diagnostics
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Pointers
    use TLab_Arrays, only: s
    use TLab_Memory, only: imax, jmax, kmax
    use NavierStokes
    use FI_GRADIENT_EQN
    implicit none
    private

    public :: Diagnose_ScalarGradientEquation

contains
    !########################################################################
    !########################################################################
    subroutine Diagnose_EnstrophyEquation(vars)
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
    subroutine Diagnose_StrainEquation(vars)
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
    subroutine Diagnose_ScalarGradientEquation(is, vars)
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
    subroutine Diagnose_MoistThermodynamics(vars)
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
    subroutine Diagnose_Damkohler(vars)
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
