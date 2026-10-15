#include "tlab_error.h"

module Radiation
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: MAX_PARS, MAX_VARS
    use TLab_Constants, only: efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Memory, only: TLab_Allocate_Real
    use FDM, only: fdm_Int0
    use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC
    use Thermo_Base, only: imixture, MIXT_TYPE_AIRWATER
    use Thermo_AirWater, only: inb_scal_e, inb_scal_qt, inb_scal_ql, inb_scal_T
    use Thermo_Anelastic, only: rbackground, Thermo_Anelastic_Weight_InPlace
    use TLab_CacheBlock
    use OPR_Radiation
    implicit none
    private

    public :: Radiation_Initialize
    public :: Radiation_Infrared_Z
    public :: Radiation_Infrared_Z_Cache

    ! -------------------------------------------------------------------
    integer, parameter :: ncomps_max = 10               ! maximum number of radiatively active components
    integer, parameter :: nbands_max = 10               ! maximum number of spectral bands
    type :: infrared_dt
        sequence
        integer type
        logical active(MAX_VARS), lpadding(3)           ! evolution equations affected by this term
        integer ncomps                                  ! number of radiatively active components
        integer nbands                                  ! number of spectral bands
        real(wp) :: kappa(ncomps_max, nbands_max)       ! mass absorption coefficients for each radiatively active component
        real(wp) beta(3, nbands_max)                    ! polynomial coefficients for band emission functions; assuming second-order polynomial
        real(wp) bcs_t(nbands_max)                      ! boundary conditions at the top, typically downward fluxes
        real(wp) bcs_b(nbands_max)                      ! boundary conditions at the bottom, upward fluxes or surface emissivity
    end type infrared_dt

    type(infrared_dt), public, protected :: infraredProps    ! Radiation parameters
    ! type(visible_dt), public, protected :: visibleProps    ! Radiation parameters

    integer, parameter, public :: TYPE_RAD_NONE = 0
    integer, parameter :: TYPE_IR_GRAY_LIQUID = 1
    integer, parameter :: TYPE_IR_GRAY = 2
    integer, parameter :: TYPE_IR_BAND = 3

    real(wp), allocatable, target :: bcs_ht(:)          ! flux boundary condition at the top of the domain
    real(wp), allocatable, target :: bcs_hb(:)          ! flux boundary condition at the bottom of the domain
    real(wp), allocatable, target :: tmp_rad(:, :)      ! 3D temporary arrays for radiation routine

    real(wp), pointer :: p_source(:) => null()

    ! -----------------------------------------------------------------------
    ! Cache blocking information; tuned to levante in a subdomain 48x32x768
    ! integer, parameter :: groupSizeX = 1024, groupSizeY = 512, groupSizeZ = 32
    integer, parameter :: groupSizeZ = CacheBlockSizeZ
    real(wp), allocatable :: ql(:, :), qv(:, :), T(:, :), a_source_x(:, :)

contains
    !########################################################################
    !########################################################################
    subroutine Radiation_Initialize(inifile)
        use TLab_Memory, only: imax, jmax, kmax
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        integer(wi) idummy, ib, ic
        integer(wi) :: inb_tmp_rad = 0
        real(wp) :: dummy(MAX_PARS)

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Infrared'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<value>')
        call TLab_Write_ASCII(bakfile, '#AbsorptionComponent#=<values>')
        call TLab_Write_ASCII(bakfile, '#BoundaryConditionsKmax=<values>')
        call TLab_Write_ASCII(bakfile, '#BoundaryConditionsKmin=<values>')
        call TLab_Write_ASCII(bakfile, '#BetaCoefficient#=<values>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; infraredProps%type = TYPE_RAD_NONE
        else if (trim(adjustl(sRes)) == 'grayliquid') then; infraredProps%type = TYPE_IR_GRAY_LIQUID
        else if (trim(adjustl(sRes)) == 'gray') then; infraredProps%type = TYPE_IR_GRAY
        else if (trim(adjustl(sRes)) == 'band') then; infraredProps%type = TYPE_IR_BAND
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Error in entry Type.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        infraredProps%active = .false.
        if (infraredProps%type /= TYPE_RAD_NONE) then
            infraredProps%bcs_t(:) = 0.0_wp
            call ScanFile_Char(bakfile, inifile, block, 'BoundaryConditionsKmax', '1.0', sRes)
            idummy = MAX_PARS
            call LIST_REAL(sRes, idummy, infraredProps%bcs_t)
            infraredProps%nbands = idummy

            infraredProps%bcs_b(:) = 0.0_wp
            call ScanFile_Char(bakfile, inifile, block, 'BoundaryConditionsKmin', '1.0', sRes)
            idummy = MAX_PARS
            call LIST_REAL(sRes, idummy, infraredProps%bcs_b)
            if (idummy == 1) then
                ! epsilon = infraredProps%bcs_b(idummy)
            else if (idummy /= infraredProps%nbands) then
                call TLab_Write_ASCII(efile, __FILE__//'. Error in '//trim(adjustl(block))//'.AbsorptionComponent.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

            infraredProps%kappa(:, :) = 0.0_wp
            do ic = 1, ncomps_max
                write (sRes, *) ic
                call ScanFile_Char(bakfile, inifile, block, 'AbsorptionComponent'//trim(adjustl(sRes)), 'void', sRes)
                if (trim(adjustl(sRes)) /= 'void') then
                    idummy = nbands_max
                    call LIST_REAL(sRes, idummy, dummy)
                    if (idummy /= infraredProps%nbands) then
                        call TLab_Write_ASCII(efile, __FILE__//'. Error in '//trim(adjustl(block))//'.AbsorptionComponent.')
                        call TLab_Stop(DNS_ERROR_OPTION)
                    end if
                    infraredProps%kappa(ic, 1:idummy) = dummy(1:idummy)
                else
                    exit
                end if
            end do
            infraredProps%ncomps = ic - 1           ! correct for the increment in the loop

            infraredProps%beta(:, :) = 0.0_wp
            infraredProps%beta(1:3, 1) = [2.6774e-1_wp, -1.3344e-3_wp, 1.8017e-6_wp] ! default coefficients for band 1 according to Jeevanjee, 2023 for vapor bands
            infraredProps%beta(1:3, 2) = [-2.2993e-2_wp, 8.7439e-5_wp, 1.4744e-7_wp] ! default coefficients for band 2
            do ic = 1, 3                                                        ! so far, only 3 coefficients, 2. order polynomial
                write (sRes, *) ic
                call ScanFile_Char(bakfile, inifile, block, 'BetaCoefficient'//trim(adjustl(sRes)), 'void', sRes)
                if (trim(adjustl(sRes)) /= 'void') then
                    idummy = nbands_max
                    call LIST_REAL(sRes, idummy, dummy)
                    if (idummy /= infraredProps%nbands - 1) then                      ! The coefficients need to add to 1, that is why we only read nbands -1
                        call TLab_Write_ASCII(efile, __FILE__//'. Error in '//trim(adjustl(block))//'.BetaCoefficient.')
                        call TLab_Stop(DNS_ERROR_OPTION)
                    end if
                    infraredProps%beta(ic, 1:idummy) = dummy(1:idummy)
                end if
            end do
            infraredProps%beta(1:3, infraredProps%nbands) = [1.0_wp, 0.0_wp, 0.0_wp]                ! last band from equation sum beta_i = 1
            do ib = 1, infraredProps%nbands - 1
                infraredProps%beta(1:3, infraredProps%nbands) = infraredProps%beta(1:3, infraredProps%nbands) - infraredProps%beta(1:3, ib)
            end do

            ! -------------------------------------------------------------------
            select case (imixture)
            case (MIXT_TYPE_AIRWATER)
                infraredProps%active(inb_scal_e) = .true.               ! energy

                ! first radiatively active scalar is liquid
                ! second radiatively active scalar is vapor
                ! third radiatively active scalar is a homogeneous field, e.g., CO2

            case default
                call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Infrared only derived for airwater mixture.')
                call TLab_Stop(DNS_ERROR_OPTION)

            end select

        end if

        ! !########################################################################
        ! bakfile = trim(adjustl(inifile))//'.bak'
        ! block = 'Visible'

        !########################################################################
        ! Local allocation
        allocate (bcs_ht(imax*jmax), bcs_hb(imax*jmax))!, t_ht(imax*jmax))
        select case (infraredProps%type)
        case (TYPE_IR_BAND)
            inb_tmp_rad = 1                                         ! Additional memory space
            call TLab_Allocate_Real(__FILE__, tmp_rad, [imax*jmax*kmax, inb_tmp_rad], 'tmp-rad')
            p_source(1:imax*jmax*kmax) => tmp_rad(1:imax*jmax*kmax, 1)
        end select

        return
    end subroutine Radiation_Initialize

    !########################################################################
    !########################################################################
    subroutine Radiation_Infrared_Z(localProps, nx, ny, nz, s, source, b, tmp1, tmp2, flux_down, flux_up)
        use TLab_Arrays, only: wrk3d
        use TLab_Constants, only: BCS_MAX, BCS_MIN, BCS_BOTH
        type(infrared_dt), intent(in) :: localProps
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(out) :: source(nx*ny*nz)               ! heating rate, also used for absorption coefficient below
        real(wp), intent(inout) :: b(nx*ny*nz)                  ! emission function
        real(wp), intent(inout) :: tmp1(nx*ny*nz), tmp2(nx*ny*nz)
        real(wp), intent(inout), optional :: flux_down(nx*ny*nz), flux_up(nx*ny*nz)

        ! -----------------------------------------------------------------------
        integer iband, ibc, nlines

        !########################################################################
        nlines = nx*ny

        ! -----------------------------------------------------------------------
        select case (localProps%type)
        case (TYPE_IR_GRAY_LIQUID)
            ! absorption coefficient
            source(1:nx*ny*nz) = localProps%kappa(1, 1)*s(:, inb_scal_ql)
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call Thermo_Anelastic_Weight_InPlace(nx, ny, nz, rbackground, source)
            end if

            ! solve radiative transfer equation
            ibc = BCS_MAX                                   ! Flux given at the top
            bcs_ht = localProps%bcs_t(1)                    ! downward flux at domain top
            if (localProps%bcs_b(1) > 0.0_wp) then
                ibc = BCS_BOTH                              ! Flux also given at the bottom
                bcs_hb = localProps%bcs_b(1)                ! upward flux at domain bottom
            end if
            if (present(flux_up)) then
                call IR_RTE1_OnlyLiquid(nlines, nz, fdm_Int0, source, wrk3d, bcs_hb, bcs_ht, ibc, flux_down, flux_up)
                ! call IR_RTE1_OnlyLiquid_UpwardFirst(nlines, nz, fdm_Int0, source,  wrk3d, bcs_hb, bcs_ht, ibc, flux_down, flux_up)
            else
                call IR_RTE1_OnlyLiquid(nlines, nz, fdm_Int0, source, wrk3d, bcs_hb, bcs_ht, ibc)
                ! call IR_RTE1_OnlyLiquid_UpwardFirst(nlines, nz, fdm_Int0, source, ibc)
            end if

            ! -----------------------------------------------------------------------
        case (TYPE_IR_GRAY)
            ! emission function, Stefan-Boltzmann law
            b(1:nx*ny*nz) = sigma*s(1:nx*ny*nz, inb_scal_T)**4

            ! absorption coefficient
            source(1:nx*ny*nz) = localProps%kappa(1, 1)*s(:, inb_scal_ql) + &
                                 localProps%kappa(2, 1)*(s(:, inb_scal_qt) - s(:, inb_scal_ql)) + &
                                 localProps%kappa(3, 1)
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call Thermo_Anelastic_Weight_InPlace(nx, ny, nz, rbackground, source)
            end if

            ! solve radiative transfer equation
            bcs_ht = localProps%bcs_t(1)                ! downward flux at domain top
            epsilon = infraredProps%bcs_b(1)            ! surface emmisivity at the bottom of the domain
            if (present(flux_up)) then
                call IR_RTE1_Global(nlines, nz, fdm_Int0, source, b, wrk3d, bcs_hb, bcs_ht, flux_down, flux_up)
            else
                call IR_RTE1_Global(nlines, nz, fdm_Int0, source, b, wrk3d, bcs_hb, bcs_ht, tmp1, tmp2)
            end if
            ! call IR_RTE1_Local(localProps, nlines, nz, fdm_Int0, source, b, flux_down, tmp2, flux_up)
            ! call IR_RTE1_Incremental(localProps, nlines, nz, fdm_Int0, source, b, flux_down, flux_up)

            ! -----------------------------------------------------------------------
        case (TYPE_IR_BAND)
            ! initialize accumulation of conttributions from each band
            source = 0.0_wp
            if (present(flux_up)) then
                flux_down = 0.0_wp
                flux_up = 0.0_wp
            end if

            do iband = 1, localProps%nbands
                ! emission function, Stefan-Boltzmann law
                b = sigma*s(:, inb_scal_T)**4* &
                    (localProps%beta(1, iband) + s(:, inb_scal_T)*(localProps%beta(2, iband) + s(:, inb_scal_T)*localProps%beta(3, iband)))

                ! absorption coefficient
                p_source(1:nx*ny*nz) = localProps%kappa(1, iband)*s(:, inb_scal_ql) + &
                                       localProps%kappa(2, iband)*(s(:, inb_scal_qt) - s(:, inb_scal_ql)) + &
                                       localProps%kappa(3, iband)
                if (nse_eqns == DNS_EQNS_ANELASTIC) then
                    call Thermo_Anelastic_Weight_InPlace(nx, ny, nz, rbackground, p_source)
                end if

                ! solve radiative transfer equation
                bcs_ht(1:nlines) = localProps%bcs_t(iband)   ! downward flux at domain top
                epsilon = infraredProps%bcs_b(1)            ! surface emmisivity at the bottom of the domain
                call IR_RTE1_Global(nlines, nz, fdm_Int0, p_source, b, wrk3d, bcs_hb, bcs_ht, tmp1, tmp2)
                ! call IR_RTE1_Local()
                ! call IR_RTE1_Incremental()

                ! accumulate
                source = source + p_source
                if (present(flux_up)) then
                    flux_down = flux_down + tmp1
                    flux_up = flux_up + tmp2
                end if

            end do

        end select

    end subroutine Radiation_Infrared_Z

    !########################################################################
    !########################################################################
    subroutine Radiation_Infrared_Z_Cache(localProps, nx, ny, nz, s, source, b, tmp1, tmp2, result, flux_down, flux_up)
        use TLab_Arrays, only: wrk3d
        use TLab_Constants, only: BCS_MAX, BCS_MIN, BCS_BOTH
        type(infrared_dt), intent(in) :: localProps
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz, *)
        real(wp), intent(inout) :: source(nx*ny*nz)               ! heating rate, also used for absorption coefficient below
        real(wp), intent(inout) :: b(nx*ny*nz)                  ! emission function
        real(wp), intent(inout) :: tmp1(nx*ny*nz), tmp2(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: flux_down(nx*ny*nz), flux_up(nx*ny*nz)

        ! -----------------------------------------------------------------------
        integer iband, ibc, nlines
        integer(wi) ib, ip

        !########################################################################
        nlines = nx*ny

        do ib = 1, nx*ny, groupSizeZ
            ip = (ib - 1)*nz + 1
            nlines = min(groupSizeZ, nx*ny - ib + 1)

            select case (localProps%type)
            case (TYPE_IR_GRAY_LIQUID)
                ! memory arrangement
                call tlab_cache_reduce_z(s(ib, inb_scal_ql), nlines, nx*ny, nz, b)

                ! absorption coefficient
                b(1:nlines*nz) = localProps%kappa(1, 1)*b(1:nlines*nz)
                if (nse_eqns == DNS_EQNS_ANELASTIC) then
                    call Thermo_Anelastic_Weight_InPlace(nlines, 1, nz, rbackground, b)
                end if

                ! solve radiative transfer equation
                ibc = BCS_MAX                                   ! Flux given at the top
                bcs_ht = localProps%bcs_t(1)                    ! downward flux at domain top
                if (localProps%bcs_b(1) > 0.0_wp) then
                    ibc = BCS_BOTH                              ! Flux also given at the bottom
                    bcs_hb = localProps%bcs_b(1)                ! upward flux at domain bottom
                end if
                if (present(flux_up)) then
                    call IR_RTE1_OnlyLiquid(nlines, nz, fdm_Int0, b, wrk3d, bcs_hb, bcs_ht, ibc, tmp1, tmp2)
                else
                    call IR_RTE1_OnlyLiquid(nlines, nz, fdm_Int0, b, wrk3d, bcs_hb, bcs_ht, ibc)
                end if

            case (TYPE_IR_GRAY)
            case (TYPE_IR_BAND)

            end select

            ! memory arrangement
            call tlab_cache_spread_add_z(b, nlines, nx*ny, nz, result(ib))
            if (present(flux_up)) then
                call tlab_cache_spread_z(tmp1, nlines, nx*ny, nz, flux_down(ib))
                call tlab_cache_spread_z(tmp2, nlines, nx*ny, nz, flux_up(ib))
            end if

        end do

    end subroutine Radiation_Infrared_Z_Cache

end module
