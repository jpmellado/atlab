#include "tlab_error.h"

module FDM
    use TLab_Constants, only: wp, wi, roundoff_wp, efile, wfile
    use TLab_Constants, only: BCS_DD, BCS_ND, BCS_DN, BCS_NN, BCS_MIN, BCS_MAX, BCS_NONE
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, stagger_on
    use TLab_Grid, only: grid_dt
    use FDM_Base
    use FDM_Derivative_1order
    use FDM_Derivative_2order
    use FDM_Integral_1
    implicit none
    private

    class(der_dt), allocatable, public, protected :: fdm_der1_X, fdm_der1_Y, fdm_der1_Z
    class(der2_dt), allocatable, public, protected :: fdm_der2_X, fdm_der2_Y, fdm_der2_Z

    type(fdm_integral_dt), public, protected :: fdm_Int0(2)    ! fdm integral plans along Oz (ode for lambda = 0)

    public :: FDM_Initialize
    public :: FDM_CreatePlan_Der1
    public :: FDM_CreatePlan_Der2

contains
    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Initialize(inifile)
        use TLab_Grid, only: x, y, z, grid
        character(len=*), optional, intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes
        integer ig, der1_type, der2_type

        !########################################################################
        ! Reading
        bakfile = trim(adjustl(inifile))//'.bak'

        ! -------------------------------------------------------------------
        block = 'Space'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#SchemeDerivative1=<CompactJacobian4/CompactJacobian6/CompactJacobian6Penta/CompactDirect4/CompactDirect6>')
call TLab_Write_ASCII(bakfile, '#SchemeDerivative2=<CompactJacobian4/CompactJacobian6/CompactJacobian6Hyper/CompactDirect4/CompactDirect6/CompactDirect6Hyper>')

        call ScanFile_Char(bakfile, inifile, block, 'SchemeDerivative1', 'compactjacobian6', sRes)
        if (trim(adjustl(sRes)) == 'compactjacobian4') then; der1_type = FDM_COM4_JACOBIAN; 
        elseif (trim(adjustl(sRes)) == 'compactjacobian6') then; der1_type = FDM_COM6_JACOBIAN; 
        elseif (trim(adjustl(sRes)) == 'compactjacobian6penta') then; der1_type = FDM_COM6_JACOBIAN_PENTA; 
        elseif (trim(adjustl(sRes)) == 'compactdirect4') then; der1_type = FDM_COM4_DIRECT; 
        elseif (trim(adjustl(sRes)) == 'compactdirect6') then; der1_type = FDM_COM6_DIRECT; 
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong SchemeDerivative1 option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'SchemeDerivative2', 'CompactJacobian6Hyper', sRes)
        if (trim(adjustl(sRes)) == 'compactjacobian4') then; der2_type = FDM_COM4_JACOBIAN; 
        elseif (trim(adjustl(sRes)) == 'compactjacobian6') then; der2_type = FDM_COM6_JACOBIAN; 
        elseif (trim(adjustl(sRes)) == 'compactjacobian6hyper') then; der2_type = FDM_COM6_JACOBIAN_HYPER; 
        elseif (trim(adjustl(sRes)) == 'compactdirect4') then; der2_type = FDM_COM4_DIRECT; 
        elseif (trim(adjustl(sRes)) == 'compactdirect6') then; der2_type = FDM_COM6_DIRECT; 
        elseif (trim(adjustl(sRes)) == 'compactdirect6hyper') then; der2_type = FDM_COM6_DIRECT_HYPER; 
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong SchemeDerivative2 option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        if (der1_type == FDM_COM6_JACOBIAN_PENTA) then     ! CFL_max depends on max[g(ig)%der1%mwn(:)]
           call TLab_Write_ASCII(wfile, trim(adjustl(eStr))//'CompactJacobian6Penta requires adjusted CFL-number depending on alpha and beta values.')
        end if

        ! -------------------------------------------------------------------
        block = 'Grid'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#XUniform=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#YUniform=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#ZUniform=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#XPeriodic=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#YPeriodic=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#ZPeriodic=<yes/no>')

        grid(1)%name = 'x'
        grid(2)%name = 'y'
        grid(3)%name = 'z'

        do ig = 1, 3
            call ScanFile_Char(bakfile, inifile, block, grid(ig)%name(1:1)//'Uniform', 'void', sRes)
            if (trim(adjustl(sRes)) == 'yes') then; grid(ig)%uniform = .true.
            else if (trim(adjustl(sRes)) == 'no') then; grid(ig)%uniform = .false.
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Error in Uniform '//grid(ig)%name(1:1)//' grid')
                call TLab_Stop(DNS_ERROR_UNIFORMX)
            end if

            call ScanFile_Char(bakfile, inifile, block, grid(ig)%name(1:1)//'Periodic', 'void', sRes)
            if (trim(adjustl(sRes)) == 'yes') then; grid(ig)%periodic = .true.
            else if (trim(adjustl(sRes)) == 'no') then; grid(ig)%periodic = .false.
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Error in Periodic '//grid(ig)%name(1:1)//' grid')
                call TLab_Stop(DNS_ERROR_IBC)
            end if

            ! consistency check
            if (grid(ig)%periodic .and. (.not. grid(ig)%uniform)) then
                call TLab_Write_ASCII(efile, __FILE__//'. Grid must be uniform in periodic direction.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

        end do

        !########################################################################
        ! Initializing fdm plan for derivatives
        call FDM_CreatePlan_Der1(x, fdm_der1_X, der1_type)
        call FDM_CreatePlan_Der1(y, fdm_der1_Y, der1_type)
        call FDM_CreatePlan_Der1(z, fdm_der1_Z, der1_type)

        call FDM_CreatePlan_Der2(x, fdm_der2_X, der2_type, fdm_der1_X)
        call FDM_CreatePlan_Der2(y, fdm_der2_Y, der2_type, fdm_der1_Y)
        call FDM_CreatePlan_Der2(z, fdm_der2_Z, der2_type, fdm_der1_Z)

        ! ###################################################################
        ! Initializing fdm plans for first-order integrals (cases lambda = 0.0_wp)
        call FDM_Int1_Initialize(fdm_der1_Z, 0.0_wp, BCS_MIN, fdm_Int0(BCS_MIN))
        call FDM_Int1_Initialize(fdm_der1_Z, 0.0_wp, BCS_MAX, fdm_Int0(BCS_MAX))

        return
    end subroutine FDM_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_CreatePlan_Der1(x, locDer, type)
        type(grid_dt), intent(in) :: x
        class(der_dt), allocatable, intent(out) :: locDer
        integer, intent(in) :: type

        ! -------------------------------------------------------------------
        real(wp), allocatable :: x_aux(:, :)        ! uniform grid to calculate Jacobians; shape to comply with %compute routines below
        real(wp), allocatable :: dx(:, :)           ! Jacobians

        integer i

        ! ###################################################################
        if (x%size == 1) return

        if (x%periodic) then
            allocate (der1_periodic :: locDer)
        else
            allocate (der1_biased :: locDer)
        end if

        ! -------------------------------------------------------------------
        ! Calculate Jacobian for the Jacobian formulations
        if (allocated(x_aux)) deallocate (x_aux)
        allocate (x_aux(1, x%size))
        if (allocated(dx)) deallocate (dx)
        allocate (dx(1, x%size))

        select type (locDer)
        type is (der1_biased)
            ! -------------------------------------------------------------------
            ! uniform grid to calculate Jacobian (used for the stencils below and also as grid spacing in the code).
            x_aux(1, :) = [(real(i - 1, wp), i=1, x%size)]
            dx(1, :) = 1.0_wp

            call locDer%initialize(x_aux(1, :), dx(1, :), type)

            ! Calculating derivative dxds
            x_aux(1, :) = x%nodes(:)                    ! I need shape (1,nx) in the procedure
            call locDer%compute(1, x_aux, dx)

        type is (der1_periodic)
            dx(1, :) = x%nodes(2) - x%nodes(1)

        end select

        ! -------------------------------------------------------------------
        ! Actual grid; possibly nonuniform
        call locDer%initialize(x%nodes, dx(1, :), type)

        return
    end subroutine

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_CreatePlan_Der2(x, locDer, type, fdm_der1)
        type(grid_dt), intent(in) :: x
        class(der2_dt), allocatable, intent(out) :: locDer
        integer, intent(in) :: type
        class(der_dt), intent(in) :: fdm_der1

        ! -------------------------------------------------------------------
        real(wp), allocatable :: x_aux(:, :)    ! uniform grid to calculate Jacobians; shape to comply with %compute routines below
        real(wp), allocatable :: dx(:, :)       ! Jacobians

        integer i

        ! ###################################################################
        if (x%size == 1) return

        if (x%periodic) then
            allocate (der2_periodic :: locDer)
        else
            allocate (der2_biased :: locDer)
        end if

        ! -------------------------------------------------------------------
        ! Calculate Jacobian for the Jacobian formulations
        if (allocated(x_aux)) deallocate (x_aux)
        allocate (x_aux(1, x%size))
        if (allocated(dx)) deallocate (dx)
        allocate (dx(2, x%size))

        select type (locDer)
        type is (der2_biased)
            ! -------------------------------------------------------------------
            ! uniform grid to calculate Jacobian (used for the stencils below and also as grid spacing in the code).
            x_aux(1, :) = [(real(i - 1, wp), i=1, x%size)]
            dx(1, :) = 1.0_wp
            dx(2, :) = 0.0_wp

            call locDer%initialize(x_aux(1, :), dx(:, :), type, uniform=.true.)

            ! Calculating derivative dxds; calculate dsdx and invert it
            call fdm_der1%compute(1, x_aux, dx(1:1, :))
            dx(1, :) = 1.0_wp/dx(1, :)

            ! Calculating derivative dx2ds2
            x_aux(1, :) = x%nodes(:)                    ! I need shape (1,nx) in the procedure
            call locDer%compute(1, x_aux, dx(2:2, :))

        type is (der2_periodic)
            dx(1, :) = x%nodes(2) - x%nodes(1)

        end select

        ! -------------------------------------------------------------------
        ! Actual grid; possibly nonuniform
        call locDer%initialize(x%nodes, dx(:, :), type, uniform=x%uniform)

        return
    end subroutine

end module FDM
