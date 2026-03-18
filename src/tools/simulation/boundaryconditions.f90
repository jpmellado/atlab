#include "tlab_error.h"

module BoundaryConditions
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_DD, BCS_DN, BCS_ND, BCS_NN, BCS_NONE, BCS_MIN, BCS_MAX, BCS_BOTH
    use TLab_Constants, only: efile, lfile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Memory, only: imax, jmax, kmax
    implicit none
    private

    public :: BcsFlowImin, BcsFlowImax, BcsFlowJmin, BcsFlowJmax, BcsFlowKmin, BcsFlowKmax
    public :: BcsScalImin, BcsScalImax, BcsScalJmin, BcsScalJmax, BcsScalKmin, BcsScalKmax

    public :: BCS_Initialize
    public :: BCS_Neumann_Z, BCS_Neumann_Z_PerVolume
    ! public :: BCS_SURFACE_Z

    ! -------------------------------------------------------------------
    integer, parameter, public :: DNS_BCS_NONE = 0
    integer, parameter, public :: DNS_BCS_DIRICHLET = 1
    integer, parameter, public :: DNS_BCS_Neumann = 2

    type bcs_dt
        ! sequence
        integer :: type = DNS_BCS_NONE
        integer :: SfcType                      ! Type of Surface Model
        real(wp) cpl                            ! Coupling parameter for surface model
        real(wp), allocatable :: ref(:, :)      ! reference fields
    end type bcs_dt

    type, extends(bcs_dt) :: bcs_flow_dt
    contains
        procedure :: initialize => initialize_flow_dt
    end type
    type, extends(bcs_dt) :: bcs_scal_dt
    contains
        procedure :: initialize => initialize_scal_dt
    end type

    type(bcs_flow_dt), allocatable :: BcsFlowImin(:), BcsFlowImax(:), BcsFlowJmin(:), BcsFlowJmax(:), BcsFlowKmin(:), BcsFlowKmax(:)
    type(bcs_scal_dt), allocatable :: BcsScalImin(:), BcsScalImax(:), BcsScalJmin(:), BcsScalJmax(:), BcsScalKmin(:), BcsScalKmax(:)

    real(wp), allocatable :: c_b(:), c_t(:)
    integer :: k_bcs_b, k_bcs_t

    ! Surface Models
    integer, parameter :: DNS_SFC_STATIC = 0
    integer, parameter :: DNS_SFC_LINEAR = 1

contains
! ###################################################################
! ###################################################################
    subroutine BCS_Initialize(inifile)
        use TLab_Grid, only: x, y, z
        use TLab_Memory, only: inb_flow, inb_scal
        use TLab_Arrays, only: wrk1d
        use FDM_Derivative_1order, only: der1_biased_extended
        use FDM, only: fdm_der1_Z
        use FDM_derivative_Neumann
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        integer iq, is
        character(len=32) str

        ! ###################################################################
        ! Flow fields
        allocate (BcsFlowImin(1:inb_flow))
        allocate (BcsFlowImax(1:inb_flow))
        allocate (BcsFlowJmin(1:inb_flow))
        allocate (BcsFlowJmax(1:inb_flow))
        allocate (BcsFlowKmin(1:inb_flow))
        allocate (BcsFlowKmax(1:inb_flow))
        do iq = 1, inb_flow
            if (.not. x%periodic) call BcsFlowImin(iq)%initialize(inifile, 'Imin', iq)
            if (.not. x%periodic) call BcsFlowImax(iq)%initialize(inifile, 'Imax', iq)
            if (.not. y%periodic) call BcsFlowJmin(iq)%initialize(inifile, 'Jmin', iq)
            if (.not. y%periodic) call BcsFlowJmax(iq)%initialize(inifile, 'Jmax', iq)
            if (.not. z%periodic) call BcsFlowKmin(iq)%initialize(inifile, 'Kmin', iq)
            if (.not. z%periodic) call BcsFlowKmax(iq)%initialize(inifile, 'Kmax', iq)
        end do

        ! ###################################################################
        ! Scalar fields
        allocate (BcsScalImin(1:inb_scal))
        allocate (BcsScalImax(1:inb_scal))
        allocate (BcsScalJmin(1:inb_scal))
        allocate (BcsScalJmax(1:inb_scal))
        allocate (BcsScalKmin(1:inb_scal))
        allocate (BcsScalKmax(1:inb_scal))
        do is = 1, inb_scal
            if (.not. x%periodic) call BcsScalImin(is)%initialize(inifile, 'Imin', is)
            if (.not. x%periodic) call BcsScalImax(is)%initialize(inifile, 'Imax', is)
            if (.not. y%periodic) call BcsScalJmin(is)%initialize(inifile, 'Jmin', is)
            if (.not. y%periodic) call BcsScalJmax(is)%initialize(inifile, 'Jmax', is)
            if (.not. z%periodic) call BcsScalKmin(is)%initialize(inifile, 'Kmin', is)
            if (.not. z%periodic) call BcsScalKmax(is)%initialize(inifile, 'Kmax', is)
        end do

        ! ###################################################################
        ! Using only truncated versions because typical decay index is 30-40
        allocate (c_b(kmax), c_t(kmax))

        select type (fdm_der1_Z)
        type is (der1_biased_extended)

            call FDM_Der1_NeumannMin_Initialize(fdm_der1_Z, c_b(:), wrk1d(1, 1), wrk1d(1, 2), k_bcs_b)
            write (str, *) k_bcs_b
            call TLab_Write_ASCII(lfile, 'Decay to round-off in bottom Neumann condition in '//trim(adjustl(str))//' indexes.')

            call FDM_Der1_NeumannMax_Initialize(fdm_der1_Z, c_t(:), wrk1d(1, 1), wrk1d(1, 2), k_bcs_t)
            write (str, *) k_bcs_t
            call TLab_Write_ASCII(lfile, 'Decay to round-off in top Neumann condition in '//trim(adjustl(str))//' indexes.')

        end select

        return
    end subroutine

! ###################################################################
! ###################################################################
    subroutine initialize_flow_dt(self, inifile, tag, idir)
        class(bcs_flow_dt), intent(inout) :: self
        character(len=*), intent(in) :: inifile
        character(len=*), intent(in) :: tag
        integer, intent(in) :: idir

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        integer inormal

        ! ###################################################################
        ! read
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'BoundaryConditions'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Velocity'//trim(adjustl(tag))//'=<none/noslip/freeslip>')

        select case (trim(adjustl(tag)))
        case ('Imin', 'Imax')
            inormal = 1                 ! velocity component normal to the boundary
        case ('Jmin', 'Jmax')
            inormal = 2
        case ('Kmin', 'Kmax')
            inormal = 3
        end select

        call ScanFile_Char(bakfile, inifile, block, 'Velocity'//trim(adjustl(tag)), 'freeslip', sRes)
        if (trim(adjustl(sRes)) == 'none') then; self%type = DNS_BCS_NONE
        else if (trim(adjustl(sRes)) == 'noslip') then; self%type = DNS_BCS_DIRICHLET
        else if (trim(adjustl(sRes)) == 'freeslip') then; self%type = DNS_BCS_NEUMANN
            if (idir == inormal) self%type = DNS_BCS_DIRICHLET   ! normal components are always nopenetration
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong '//trim(adjustl(tag)))
            call TLab_Stop(DNS_ERROR_IBC)
        end if

        ! ###################################################################
        ! initialize
        select case (trim(adjustl(tag)))
        case ('Imin', 'Imax')
            allocate (self%ref(jmax, kmax))
        case ('Jmin', 'Jmax')
            allocate (self%ref(imax, kmax))
        case ('Kmin', 'Kmax')
            allocate (self%ref(imax, jmax))
        end select

        return
    end subroutine

! ###################################################################
! ###################################################################
    subroutine initialize_scal_dt(self, inifile, tag, is)
        class(bcs_scal_dt), intent(inout) :: self
        character(len=*), intent(in) :: inifile
        character(len=*), intent(in) :: tag
        integer, intent(in) :: is

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block, locScal
        character(len=128) eStr
        character(len=512) sRes

        ! ###################################################################
        ! read
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'BoundaryConditions'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        write (locScal, *) is; locScal = 'Scalar'//trim(adjustl(locScal))

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#'//trim(adjustl(locScal))//trim(adjustl(tag))//'=<none/dirichlet/neumman>')
        ! call TLab_Write_ASCII(bakfile, '#'//trim(adjustl(locScal))//'SfcType'//trim(adjustl(tag))//'=<static/linear>')
        ! call TLab_Write_ASCII(bakfile, '#'//trim(adjustl(locScal))//'Coupling'//trim(adjustl(tag))//'=<value>')

        call ScanFile_Char(bakfile, inifile, block, trim(adjustl(locScal))//trim(adjustl(tag)), 'void', sRes)
        if (trim(adjustl(sRes)) == 'none') then; self%type = DNS_BCS_NONE
        else if (trim(adjustl(sRes)) == 'dirichlet') then; self%type = DNS_BCS_DIRICHLET
        else if (trim(adjustl(sRes)) == 'neumann') then; self%type = DNS_BCS_Neumann
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong '//trim(adjustl(locScal))//trim(adjustl(tag)))
            call TLab_Stop(DNS_ERROR_IBC)
        end if

        ! call ScanFile_Char(bakfile, inifile, block, trim(adjustl(locScal))//'SfcType'//trim(adjustl(tag)), 'static', sRes)
        ! if (sRes == 'static') then
        !     self%%SfcType = DNS_SFC_STATIC
        ! elseif (sRes == 'linear') then
        !     self%%SfcType = DNS_SFC_LINEAR
        ! else
        !     call TLab_Write_ASCII(efile, __FILE__//'. BoundaryConditions.'//trim(adjustl(lstr))//'SfcType'//trim(adjustl(tag)))
        !     call TLab_Stop(DNS_ERROR_JBC)
        ! end if

        ! call ScanFile_Real(bakfile, inifile, block, trim(adjustl(locScal))//'Coupling'//trim(adjustl(tag)), '0.0', self%%cpl)

        ! ###################################################################
        ! initialize
        select case (trim(adjustl(tag)))
        case ('Imin', 'Imax')
            allocate (self%ref(jmax, kmax))
        case ('Jmin', 'Jmax')
            allocate (self%ref(imax, kmax))
        case ('Kmin', 'Kmax')
            allocate (self%ref(imax, jmax))
        end select

        return
    end subroutine

    !########################################################################
    !########################################################################
    !# Calculate the boundary values of a field s.t. the normal derivative is zero
    !# see valid/fdm
    subroutine BCS_Neumann_Z(ibc, nlines, nz, u, bcs_hb, bcs_ht)
        integer(wi), intent(in) :: ibc     ! BCs at Kmin/Kmax: 1, for Neumann/-
        !                                                      2, for -      /Neumann
        !                                                      3, for Neumann/Neumann
        integer(wi) nlines, nz
        real(wp), intent(in) :: u(nlines, nz)
        real(wp), intent(out) :: bcs_hb(nlines), bcs_ht(nlines)

        integer k

        ! ###################################################################
        if (any([BCS_ND, BCS_NN] == ibc)) then
            bcs_hb(:) = 0.0_wp                  ! zero derivative at the wall
            do k = 2, k_bcs_b
                bcs_hb(:) = bcs_hb(:) + c_b(k)*u(:, k)
            end do
        end if
        if (any([BCS_DN, BCS_NN] == ibc)) then
            bcs_ht(:) = 0.0_wp                  ! zero derivative at the wall
            do k = nz - 1, nz - k_bcs_t + 1, -1
                bcs_ht(:) = bcs_ht(:) + c_t(k)*u(:, k)
            end do
        end if

        return
    end subroutine BCS_Neumann_Z

    !########################################################################
    !########################################################################
    subroutine BCS_Neumann_Z_PerVolume(ibc, nlines, nz, u, bcs_hb, bcs_ht)
        use Thermo_Anelastic, only: ribackground

        integer(wi), intent(in) :: ibc     ! BCs at Kmin/Kmax: 1, for Neumann/-
        !                                                      2, for -      /Neumann
        !                                                      3, for Neumann/Neumann
        integer(wi) nlines, nz
        real(wp), intent(in) :: u(nlines, nz)
        real(wp), intent(out) :: bcs_hb(nlines), bcs_ht(nlines)

        integer k

        ! ###################################################################
        if (any([BCS_ND, BCS_NN] == ibc)) then
            bcs_hb(:) = 0.0_wp                  ! zero derivative at the wall
            do k = 2, k_bcs_b
                bcs_hb(:) = bcs_hb(:) + c_b(k)*u(:, k)*ribackground(k)
            end do
        end if
        if (any([BCS_DN, BCS_NN] == ibc)) then
            bcs_ht(:) = 0.0_wp                  ! zero derivative at the wall
            do k = nz - 1, nz - k_bcs_t + 1, -1
                bcs_ht(:) = bcs_ht(:) + c_t(k)*u(:, k)*ribackground(k)
            end do
        end if

        return
    end subroutine BCS_Neumann_Z_PerVolume

end module BoundaryConditions
