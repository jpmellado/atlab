#include "tlab_error.h"

module BoundaryConditions
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_DD, BCS_DN, BCS_ND, BCS_NN, BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH, MAX_VARS
    use TLab_Constants, only: efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use FDM, only: g
    implicit none
    private

    public :: BCS_Initialize
    public :: BCS_Neumann_Z
    public :: BCS_SURFACE_Z

    type bcs_dt
        sequence
        integer type(MAX_VARS)                              ! dirichlet, neumann for incompressible
        integer SfcType(MAX_VARS)                           ! Type of Surface Model
        real(wp) cpl(MAX_VARS)                              ! Coupling parameter for surface model
        ! real(wp) cinf, cout, ctan                           ! characteristic formulation for compressible
        real(wp), allocatable, dimension(:, :, :) :: ref    ! reference fields
    end type bcs_dt

    type(bcs_dt), public :: BcsFlowImin, BcsFlowImax, BcsFlowJmin, BcsFlowJmax, BcsFlowKmin, BcsFlowKmax
    type(bcs_dt), public :: BcsScalImin, BcsScalImax, BcsScalJmin, BcsScalJmax, BcsScalKmin, BcsScalKmax

    logical, public :: BcsDrift

    ! -------------------------------------------------------------------
    integer(wi) nmin, nmax, nsize

    ! Boundary conditions
    integer, parameter, public :: DNS_BCS_NONE = 0
    integer, parameter, public :: DNS_BCS_DIRICHLET = 1
    integer, parameter, public :: DNS_BCS_Neumann = 2

    ! Surface Models
    integer, parameter, public :: DNS_SFC_STATIC = 0
    integer, parameter, public :: DNS_SFC_LINEAR = 1

contains
! ###################################################################
! ###################################################################
    subroutine BCS_Initialize(inifile)
        use TLab_Memory, only: imax, jmax, kmax, inb_flow_array, inb_scal_array, inb_scal
        use NavierStokes

        character*(*) inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr

        integer is

        ! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'BoundaryConditions'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#VelocityImin=<none/dirichlet/neumman>')
        call TLab_Write_ASCII(bakfile, '#VelocityJmin=<none/dirichlet/neumman>')
        call TLab_Write_ASCII(bakfile, '#VelocityKmin=<none/dirichlet/neumman>')
        call TLab_Write_ASCII(bakfile, '#ScalarImin=<none/dirichlet/neumman>')
        call TLab_Write_ASCII(bakfile, '#ScalarJmin=<none/dirichlet/neumman>')
        call TLab_Write_ASCII(bakfile, '#ScalarKmin=<none/dirichlet/neumman>')

        call TLab_Write_ASCII(bakfile, '#ScalarSfcTypeKmin=<static/linear>')
        call TLab_Write_ASCII(bakfile, '#ScalarSfcTypeKmax=<static/linear>')
        call TLab_Write_ASCII(bakfile, '#ScalarCouplingKmin=<value>')
        call TLab_Write_ASCII(bakfile, '#ScalarCouplingKmax=<value>')
        
        ! -------------------------------------------------------------------
        ! Scalar terms (including surface model at vertical boundaries)
        ! -------------------------------------------------------------------
        BcsScalImin%type(:) = DNS_BCS_NONE; BcsScalImax%type(:) = DNS_BCS_NONE
        if (.not. g(1)%periodic) then
            call BCS_Scal_ReadBlock(bakfile, inifile, block, 'Imin', BcsScalImin)
            call BCS_Scal_ReadBlock(bakfile, inifile, block, 'Imax', BcsScalImax)
        end if

        BcsScalJmin%type(:) = DNS_BCS_NONE; BcsScalJmax%type(:) = DNS_BCS_NONE
        if (.not. g(2)%periodic) then
            call BCS_Scal_ReadBlock(bakfile, inifile, block, 'Jmin', BcsScalJmin)
            call BCS_Scal_ReadBlock(bakfile, inifile, block, 'Jmax', BcsScalJmax)
        end if

        BcsScalKmin%type(:) = DNS_BCS_NONE; BcsScalKmax%type(:) = DNS_BCS_NONE
        if (.not. g(3)%periodic) then
            call BCS_Scal_ReadBlock(bakfile, inifile, block, 'Kmin', BcsScalKmin)
            call BCS_Scal_ReadBlock(bakfile, inifile, block, 'Kmax', BcsScalKmax)
        end if

        ! -------------------------------------------------------------------
        ! Velocity terms / Euler part in compressible mode
        ! -------------------------------------------------------------------
        call BCS_Flow_ReadBlock(bakfile, inifile, block, 'Imin', BcsFlowImin)
        call BCS_Flow_ReadBlock(bakfile, inifile, block, 'Imax', BcsFlowImax)
        call BCS_Flow_ReadBlock(bakfile, inifile, block, 'Jmin', BcsFlowJmin)
        call BCS_Flow_ReadBlock(bakfile, inifile, block, 'Jmax', BcsFlowJmax)
        call BCS_Flow_ReadBlock(bakfile, inifile, block, 'Kmin', BcsFlowKmin)
        call BCS_Flow_ReadBlock(bakfile, inifile, block, 'Kmax', BcsFlowKmax)

        ! -------------------------------------------------------------------
        ! Boundary conditions
        ! -------------------------------------------------------------------
        ! Make sure periodic BCs are not modified
        if (g(1)%periodic) then; 
            BcsFlowImin%type(:) = DNS_BCS_NONE; BcsFlowImax%type(:) = DNS_BCS_NONE
            BcsScalImin%type(:) = DNS_BCS_NONE; BcsScalImax%type(:) = DNS_BCS_NONE
        end if
        if (g(2)%periodic) then; 
            BcsFlowJmin%type(:) = DNS_BCS_NONE; BcsFlowJmax%type(:) = DNS_BCS_NONE
            BcsScalJmin%type(:) = DNS_BCS_NONE; BcsScalJmax%type(:) = DNS_BCS_NONE
        end if
        if (g(3)%periodic) then; 
            BcsFlowKmin%type(:) = DNS_BCS_NONE; BcsFlowKmax%type(:) = DNS_BCS_NONE
            BcsScalKmin%type(:) = DNS_BCS_NONE; BcsScalKmax%type(:) = DNS_BCS_NONE
        end if

        ! -------------------------------------------------------------------
        ! Interactive Boundary conditions
        ! -------------------------------------------------------------------
        do is = 1, inb_scal
            if (BcsScalJmin%type(is) /= DNS_BCS_DIRICHLET .and. &
                BcsScalJmin%SfcType(is) /= DNS_SFC_STATIC) then
                call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Interactive BC at kmin not implemented for non-Dirichlet BC')
                call TLab_Stop(DNS_ERROR_JBC)
            end if
            if (BcsScalJmax%type(is) /= DNS_BCS_DIRICHLET .and. &
                BcsScalJmax%SfcType(is) /= DNS_SFC_STATIC) then
                write (*, *) BcsScalJmax%type(is), BcsScalJmax%SfcType(is), BcsScalJmax%cpl(is)
                call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Interactive BC at kmax not implemented for non-Dirichlet BC')
                call TLab_Stop(DNS_ERROR_JBC)
            end if
        end do

        !########################################################################
        select case (nse_eqns)
        case (DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC)
            ! we use inb_flow_array and inb_scal_array to have space for aux fields
            ! we add space at the end for the shape factor
            allocate (BcsFlowImin%ref(jmax, kmax, inb_flow_array + 1))
            allocate (BcsFlowImax%ref(jmax, kmax, inb_flow_array + 1))
            allocate (BcsFlowJmin%ref(imax, kmax, inb_flow_array + 1))
            allocate (BcsFlowJmax%ref(imax, kmax, inb_flow_array + 1))
            allocate (BcsFlowKmin%ref(imax, jmax, inb_flow_array + 1))
            allocate (BcsFlowKmax%ref(imax, jmax, inb_flow_array + 1))

            allocate (BcsScalImin%ref(jmax, kmax, inb_scal_array + 1))
            allocate (BcsScalImax%ref(jmax, kmax, inb_scal_array + 1))
            allocate (BcsScalJmin%ref(imax, kmax, inb_scal_array + 1))
            allocate (BcsScalJmax%ref(imax, kmax, inb_scal_array + 1))
            allocate (BcsScalKmin%ref(imax, jmax, inb_scal_array + 1))
            allocate (BcsScalKmax%ref(imax, jmax, inb_scal_array + 1))

        case (DNS_EQNS_COMPRESSIBLE)

        end select

        return
    end subroutine BCS_Initialize

    !########################################################################
    !########################################################################
    !# Calculate the boundary values of a field s.t. the normal derivative is zero
    !# Routine format extracted from OPR_Partial_Y
    subroutine BCS_Neumann_Z(ibc, nx, ny, nz, u, bcs_hb, bcs_ht, tmp1)
        integer(wi), intent(in) :: ibc     ! BCs at jmin/jmax: 1, for Neumann/-
        !                                                   2, for -      /Neumann
        !                                                   3, for Neumann/Neumann
        integer(wi) nx, ny, nz
        real(wp), intent(in) :: u(nx*ny, nz)            ! they are transposed below
        real(wp), intent(inout) :: tmp1(nx*ny, nz)      ! they are transposed below
        real(wp), intent(out) :: bcs_hb(nx*ny), bcs_ht(nx*ny)

        ! -------------------------------------------------------------------
        integer(wi) ip, idl, ic

        ! ###################################################################
        if (g(3)%size == 1) then ! Set to zero in 2D case
            bcs_hb = 0.0_wp; bcs_ht = 0.0_wp
            return
        end if

        ! ###################################################################
        ip = ibc*5

        nmin = 1; nmax = g(3)%size
        if (any([BCS_ND, BCS_NN] == ibc)) then
            tmp1(:, 1) = 0.0_wp      ! homogeneous bcs
            nmin = nmin + 1
        end if
        if (any([BCS_DN, BCS_NN] == ibc)) then
            tmp1(:, nz) = 0.0_wp
            nmax = nmax - 1
        end if
        nsize = nmax - nmin + 1

        ! -------------------------------------------------------------------
        ! Calculate RHS in system of equations A u' = B u
        call g(3)%der1%matmul(g(3)%der1%rhs, u, tmp1, ibc, g(3)%der1%rhs_b, g(3)%der1%rhs_t, bcs_hb, bcs_ht)

        ! -------------------------------------------------------------------
        ! Solve for u' in system of equations A u' = B u
        select case (g(3)%der1%nb_diag(1))
        case (3)
            call TRIDSS(nsize, nx*ny, g(3)%der1%lu(nmin:nmax, ip + 1), g(3)%der1%lu(nmin:nmax, ip + 2), g(3)%der1%lu(nmin:nmax, ip + 3), tmp1(:, nmin:nmax))
        case (5)
                call PENTADSS2(nsize, nx*ny, g(3)%der1%lu(nmin:nmax, ip + 1), g(3)%der1%lu(nmin:nmax, ip + 2), g(3)%der1%lu(nmin:nmax, ip + 3), g(3)%der1%lu(nmin:nmax, ip + 4), g(3)%der1%lu(nmin:nmax, ip + 5), tmp1(:, nmin:nmax))
        end select

        idl = g(3)%der1%nb_diag(1)/2 + 1
        if (any([BCS_ND, BCS_NN] == ibc)) then
            do ic = 1, idl - 1
                bcs_hb(:) = bcs_hb(:) + g(3)%der1%lu(1, ip + idl + ic)*tmp1(:, 1 + ic)
            end do
        end if
        if (any([BCS_DN, BCS_NN] == ibc)) then
            do ic = 1, idl - 1
                bcs_ht(:) = bcs_ht(:) + g(3)%der1%lu(nz, ip + idl - ic)*tmp1(:, nz - ic)
            end do
        end if

        return
    end subroutine BCS_Neumann_Z

    !########################################################################
    !########################################################################
    ! Calculates and updates interactive surface boundary condition
    subroutine BCS_SURFACE_Z(is, s, hs, tmp1, aux)
#ifdef TRACE_ON
        use TLab_Constants, only: tfile
        use TLab_WorkFlow, only: TLab_Write_ASCII
#endif
        use TLab_Constants, only: lfile
        use TLab_Memory, only: imax, jmax, kmax
        use TLab_Memory, only: isize_field
        use NavierStokes, only: visc, schmidt
        use Averages, only: AVG1V2D
        use OPR_Partial

        integer(wi) is
        real(wp), dimension(isize_field, *) :: s, hs
        real(wp), dimension(isize_field) :: tmp1
        real(wp), dimension(imax, kmax, 6), target :: aux

        integer(wi) nxy, ip, k
        real(wp), dimension(:, :), pointer :: hfx, hfx_anom
        real(wp) :: diff, hfx_avg

#ifdef TRACE_ON
        call TLab_Write_ASCII(tfile, 'ENTERING SUBROUTINE BCS_SURFACE_Y')
#endif
        diff = visc/schmidt(is)
        nxy = imax*jmax

        ! vertical derivative of scalar for flux at the boundaries
        call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, s(:, is), tmp1)

        ! ------------------------------------------------------------
        ! Bottom Boundary
        ! ------------------------------------------------------------
        if (BcsScalJmin%SfcType(is) == DNS_SFC_LINEAR) then
            hfx => aux(:, :, 1)
            hfx_anom => aux(:, :, 2)
            ip = 1
            do k = 1, kmax    ! Calculate the surface flux
                hfx(:, k) = diff*tmp1(ip:ip + imax - 1); ip = ip + nxy
            end do
            hfx_avg = diff*AVG1V2D(imax, jmax, kmax, 1, 1, tmp1)
            hfx_anom = hfx - hfx_avg
            BcsScalJmin%ref(:, :, is) = BcsScalJmin%ref(:, :, is) + BcsScalJmin%cpl(is)*hfx_anom
        end if

        ! ------------------------------------------------------------
        ! Top Boundary
        ! ------------------------------------------------------------
        if (BcsScalJmax%SfcType(is) == DNS_SFC_LINEAR) then
            hfx => aux(:, :, 3)
            hfx_anom => aux(:, :, 4)
            ip = imax*(jmax - 1) + 1
            do k = 1, kmax; ! Calculate the surface flux
                hfx(:, k) = -diff*tmp1(ip:ip + imax - 1); ip = ip + nxy; 
            end do
            hfx_avg = diff*AVG1V2D(imax, jmax, kmax, 1, 1, tmp1)
            hfx_anom = hfx - hfx_avg
            BcsScalJmax%ref(:, :, is) = BcsScalJmax%ref(:, :, is) + BcsScalJmax%cpl(is)*hfx_anom
        end if

#ifdef TRACE_ON
        call TLab_Write_ASCII(TFILE, 'LEAVING SUBROUTINE BOUNDAR_SURFACE_J')
#endif

        return

    end subroutine BCS_SURFACE_Z

! ###################################################################
! ###################################################################
    subroutine BCS_Scal_ReadBlock(bakfile, inifile, block, tag, locVar)
        use TLab_Memory, only: inb_scal

        character(len=*) bakfile, inifile, block, tag
        type(bcs_dt), intent(out) :: locVar

        character(len=512) sRes
        character(len=128) eStr, lstr
        integer is

        ! ###################################################################
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        do is = 1, inb_scal
            write (lstr, *) is; lstr = 'Scalar'//trim(adjustl(lstr))

            call ScanFile_Char(bakfile, inifile, block, trim(adjustl(lstr))//trim(adjustl(tag)), 'void', sRes)
            if (trim(adjustl(sRes)) == 'none') then; locVar%type(is) = DNS_BCS_NONE
            else if (trim(adjustl(sRes)) == 'dirichlet') then; locVar%type(is) = DNS_BCS_DIRICHLET
            else if (trim(adjustl(sRes)) == 'neumann') then; locVar%type(is) = DNS_BCS_Neumann
            else
                call TLab_Write_ASCII(efile, __FILE__//'. BoundaryConditions.'//trim(adjustl(lstr)))
                call TLab_Stop(DNS_ERROR_JBC)
            end if

            call ScanFile_Char(bakfile, inifile, block, trim(adjustl(lstr))//'SfcType'//trim(adjustl(tag)), 'static', sRes)
            if (sRes == 'static') then
                locVar%SfcType(is) = DNS_SFC_STATIC
            elseif (sRes == 'linear') then
                locVar%SfcType(is) = DNS_SFC_LINEAR
            else
                call TLab_Write_ASCII(efile, __FILE__//'. BoundaryConditions.'//trim(adjustl(lstr))//'SfcType'//trim(adjustl(tag)))
                call TLab_Stop(DNS_ERROR_JBC)
            end if

            call ScanFile_Real(bakfile, inifile, block, trim(adjustl(lstr))//'Coupling'//trim(adjustl(tag)), '0.0', locVar%cpl(is))

        end do

        return
    end subroutine BCS_Scal_ReadBlock

! ###################################################################
! ###################################################################
    subroutine BCS_Flow_ReadBlock(bakfile, inifile, block, tag, locVar)
        character(len=*) bakfile, inifile, block, tag
        type(bcs_dt) locVar

        character(len=512) sRes
        character(len=128) eStr
        integer inormal, itangential(2)

        ! ###################################################################
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        select case (trim(adjustl(tag)))
        case ('Imin', 'Imax')
            inormal = 1
            itangential = [2, 3]
        case ('Jmin', 'Jmax')
            inormal = 2
            itangential = [1, 3]
        case ('Kmin', 'Kmax')
            inormal = 3
            itangential = [1, 2]
        end select

        call ScanFile_Char(bakfile, inifile, block, 'Velocity'//trim(adjustl(tag)), 'freeslip', sRes)
        if (trim(adjustl(sRes)) == 'none') then; locVar%type(1:3) = DNS_BCS_NONE
        else if (trim(adjustl(sRes)) == 'noslip') then; locVar%type(1:3) = DNS_BCS_DIRICHLET
        else if (trim(adjustl(sRes)) == 'freeslip') then; locVar%type(inormal) = DNS_BCS_DIRICHLET
            locVar%type(itangential) = DNS_BCS_Neumann
        else
            call TLab_Write_ASCII(efile, __FILE__//'. BoundaryConditions.Velocity'//trim(adjustl(tag)))
            call TLab_Stop(DNS_ERROR_IBC)
        end if

    end subroutine BCS_Flow_ReadBlock

end module BoundaryConditions
