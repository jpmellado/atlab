#include "tlab_error.h"

module FDM
    use TLab_Constants, only: wp, wi, roundoff_wp, efile, wfile
    use TLab_Constants, only: BCS_DD, BCS_ND, BCS_DN, BCS_NN, BCS_MIN, BCS_MAX, BCS_NONE
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, stagger_on
    use TLab_Grid, only: grid_dt
    use FDM_Derivative
    use FDM_Derivative_1order_X
    use FDM_Derivative_2order_X
    use FDM_Interpolate
    use FDM_Integral
    implicit none
    private

    ! type, public :: fdm_dt
    !     sequence
    !     character*8 name
    !     integer(wi) size
    !     logical :: uniform = .false.
    !     logical :: periodic = .false.
    !     real(wp) scale
    !     !
    !     real(wp), allocatable :: nodes(:)
    !     real(wp), allocatable :: jac(:, :)      ! Jacobian of 1. order derivative; grid spacing
    !     !                                         Jacobian of 2. order derivative; grid stretching
    !     type(fdm_derivative_dt) :: der1
    !     type(fdm_derivative_dt) :: der2
    !     type(fdm_interpol_dt) :: intl

    ! end type fdm_dt

    ! type(fdm_dt), public, protected :: g(3)                    ! fdm derivative plans along 3 directions
    type(fdm_integral_dt), public, protected :: fdm_Int0(2)    ! fdm integral plans along Oz (ode for lambda = 0)

    class(der_dt), allocatable, public, protected :: fdm_der1_X, fdm_der1_Y, fdm_der1_Z
    class(der2_dt), allocatable, public, protected :: fdm_der2_X, fdm_der2_Y, fdm_der2_Z

    public :: FDM_Initialize
    ! public :: FDM_CreatePlan
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

        ! ! old
        ! g(1:3)%der1%mode_fdm = der1_type
        ! g(1:3)%der2%mode_fdm = der2_type
        ! g(1)%periodic = x%periodic
        ! g(2)%periodic = y%periodic
        ! g(3)%periodic = z%periodic
        ! g(1)%uniform = x%uniform
        ! g(2)%uniform = y%uniform
        ! g(3)%uniform = z%uniform

        ! call FDM_CreatePlan(x, g(1))
        ! call FDM_CreatePlan(y, g(2))
        ! call FDM_CreatePlan(z, g(3))

        ! ###################################################################
        ! Initializing fdm plans for first-order integrals (cases lambda = 0.0_wp)
        ! call FDM_Int1_Initialize(g(3)%der1, 0.0_wp, BCS_MIN, fdm_Int0(BCS_MIN))
        ! call FDM_Int1_Initialize(g(3)%der1, 0.0_wp, BCS_MAX, fdm_Int0(BCS_MAX))
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

        ! select type (locDer)
        ! type is (der1_periodic)
        !     call FDM_Der1_ModifyWavenumbers(size(locDer%lhs, 1), locDer%lhs(1, :), locDer%rhs(1, :), locDer%mwn)
        ! end select

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

        ! select type (locDer)
        ! type is (der2_periodic)
        !     call FDM_Der2_ModifyWavenumbers(size(locDer%lhs, 1), locDer%lhs(1, :), locDer%rhs(1, :), locDer%mwn)
        ! end select

        return
    end subroutine

!     ! ###################################################################
!     ! ###################################################################
!     subroutine FDM_CreatePlan(x, g, locScale)
!         type(grid_dt), intent(in) :: x
!         type(fdm_dt), intent(inout) :: g                            ! fdm plan for derivatives
!         real(wp), intent(in), optional :: locScale                  ! for consistency check

!         ! -------------------------------------------------------------------
!         integer(wi) i, nx

!         ! ###################################################################
!         ! Consistency check
!         ! ###################################################################
!         if (g%periodic) then            ! if periodic, the grid is uniform,
!             !                             if uniform, direct and Jacobian formulations are the same,
!             !                             and periodic case is only implemented in Jacobian.
!             !                             Maybe one could call the Jacobian routines inside the direct ones
!             if (g%der1%mode_fdm == FDM_COM4_DIRECT) g%der1%mode_fdm = FDM_COM4_JACOBIAN
!             if (g%der1%mode_fdm == FDM_COM6_DIRECT) g%der1%mode_fdm = FDM_COM6_JACOBIAN

!             if (g%der2%mode_fdm == FDM_COM4_DIRECT) g%der2%mode_fdm = FDM_COM4_JACOBIAN
!             if (g%der2%mode_fdm == FDM_COM6_DIRECT) g%der2%mode_fdm = FDM_COM6_JACOBIAN

!         end if
!         if (g%der2%mode_fdm == FDM_COM6_DIRECT_HYPER) g%der2%mode_fdm = FDM_COM6_JACOBIAN_HYPER

!         g%size = x%size

!         if (g%size > 1) then
!             g%scale = x%nodes(g%size) - x%nodes(1)
!             if (g%periodic) g%scale = g%scale*(1.0_wp + 1.0_wp/real(g%size - 1, wp))
!         else
!             g%scale = 1.0_wp            ! to avoid conditionals and NaN in some of the calculations below
!         end if

!         if (present(locScale)) then
! ! print *, abs((scale_loc - g%scale)/scale_loc)
!             if (abs((locScale - g%scale)/g%scale) > roundoff_wp) then
!                 call TLab_Write_ASCII(efile, __FILE__//'. Unmatched domain scale.')
!                 call TLab_Stop(DNS_ERROR_OPTION)
!             end if
!         end if

!         ! ###################################################################
!         nx = g%size                     ! # of nodes, for clarity below

!         if (allocated(g%nodes)) deallocate (g%nodes)
!         if (allocated(g%jac)) deallocate (g%jac)
!         allocate (g%nodes(nx))
!         allocate (g%jac(nx, 2 + 1))     ! I need 1 aux array to calculate the Jacobian for 2. order derivative; to be fixed

!         if (nx == 1) then
!             g%jac(:, :) = 1.0_wp
!             return
!         end if

!         ! ###################################################################
!         ! first-order derivative
!         ! ###################################################################
!         ! uniform grid to calculate Jacobian (used for the stencils below and also as grid spacing in the code).
!         g%nodes(:) = [(real(i - 1, wp), i=1, g%size)]
!         g%jac(:, 1) = 1.0_wp

!         call FDM_Der1_Initialize(g%nodes, g%jac(:, 1), g%der1, periodic=.false., bcs_cases=[BCS_DD])

!         ! Calculating derivative dxds into g%jac(:, 1)
!         g%der1%periodic = .false.
!         call FDM_Der1_Solve(1, g%der1, g%der1%lu, x%nodes, g%jac(:, 1), g%jac(:, 2)) !g%jac(:, 2) is used as aux array...

!         ! -------------------------------------------------------------------
!         ! Actual grid; possibly nonuniform
!         g%nodes(:) = x%nodes(1:nx)

!         call FDM_Der1_Initialize(g%nodes, g%jac(:, 1), g%der1, periodic=g%periodic, bcs_cases=[BCS_DD, BCS_ND, BCS_DN, BCS_NN])

!         if (g%periodic) g%der1%mwn(:) = g%der1%mwn(:)/g%jac(1, 1)           ! normalized by dx

!         ! ###################################################################
!         ! second-order derivative
!         ! ###################################################################
!         ! -------------------------------------------------------------------
!         ! uniform grid to calculate Jacobian (needed to set up the Jacobian formulations and also a grid stretching in the code).
!         g%nodes(:) = [(real(i - 1, wp), i=1, g%size)]
!         g%jac(:, 2) = 1.0_wp
!         g%jac(:, 3) = 0.0_wp

!         call FDM_Der2_Initialize(g%nodes, g%jac(:, 2:), g%der2, periodic=.false., uniform=.true.)

!         ! Calculating derivative d2xds2 into g%jac(:, 3)
!         g%der2%periodic = .false.
!         call FDM_Der2_Solve(1, g%der2, g%der2%lu, x%nodes, g%jac(:, 2), g%jac(:, 3), g%jac(:, 3)) !g%jac(:, 3) is used as aux array...

!         ! -------------------------------------------------------------------
!         ! Actual grid; possibly nonuniform
!         g%nodes(:) = x%nodes(1:nx)

!         call FDM_Der2_Initialize(g%nodes, g%jac, g%der2, g%periodic, g%uniform)

!         if (g%der2%periodic) g%der2%mwn(:) = g%der2%mwn(:)/(g%jac(1, 1)**2)      ! normalized by dx

!         ! ###################################################################
!         ! interpolation for staggered cases
!         ! ###################################################################
!         if (stagger_on) then
!             if (g%periodic) then

!                 call FDM_Interpol_Initialize(x%nodes(:), g%jac(:, 1), g%intl, g%der1%mwn(:))

!                 ! else
!                 !     call TLab_Write_ASCII(efile, 'Staggered grid only along periodic directions.')
!                 !     call TLab_Stop(DNS_ERROR_UNDEVELOP)

!             end if

!         end if

!         return
!     end subroutine FDM_CreatePlan

end module FDM
