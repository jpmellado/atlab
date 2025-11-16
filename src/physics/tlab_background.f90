#include "tlab_error.h"

! information to set up bcs, ics, and reference background profiles

module Tlab_Background
    use TLab_Constants, only: wp, wi, efile, lfile, wfile, MAX_VARS
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use Profiles, only: profiles_dt, Profiles_Calculate
    implicit none
    private

    public :: TLab_Initialize_Background

    type(profiles_dt), public :: qbg(3)                     ! Velocity
    type(profiles_dt), public :: sbg(MAX_VARS)              ! Scalars
    type(profiles_dt), public :: pbg, rbg, tbg, hbg         ! Pressure, density, temperature, enthalpy

    ! background, reference profiles
    real(wp), allocatable :: sbackground(:, :)              ! Scalars

contains
    !########################################################################
    !# Initialize data of reference profiles
    !########################################################################
    subroutine TLab_Initialize_Background(inifile)
        use TLab_Arrays, only: wrk1d
        use TLab_Memory, only: inb_scal, inb_scal_array
        use TLab_Grid, only: z
        use FDM, only: fdm_Int0
        use NavierStokes, only: schmidt
        use Thermodynamics, only: imode_thermo, THERMO_TYPE_ANELASTIC
        use Thermo_Base, only: imixture
        use Thermo_Base, only: MIXT_TYPE_AIR, MIXT_TYPE_AIRVAPOR, MIXT_TYPE_AIRWATER
        use Thermo_Anelastic
        use Profiles, only: Profiles_ReadBlock
        use Gravity

        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=64) lstr

        integer(wi) is, k

        ! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'

        ! -----------------------------------------------------------------------
        block = 'Scalar'

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        do is = 1, MAX_VARS
            write (lstr, *) is
            call Profiles_ReadBlock(bakfile, inifile, block, 'Scalar'//trim(adjustl(lstr)), sbg(is))
        end do

        ! -----------------------------------------------------------------------
        block = 'Flow'

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'VelocityX', qbg(1))
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'VelocityY', qbg(2))
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'VelocityZ', qbg(3))

        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Pressure', pbg)
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Density', rbg)
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Temperature', tbg)
        call Profiles_ReadBlock(bakfile, inifile, 'Flow', 'Enthalpy', hbg)

        ! -----------------------------------------------------------------------
        do is = 1, size(qbg)
            if (qbg(is)%relative) qbg(is)%zmean = z%nodes(1) + z%scale*qbg(is)%zmean_rel
        end do
        if (pbg%relative) pbg%zmean = z%nodes(1) + z%scale*pbg%zmean_rel
        if (rbg%relative) rbg%zmean = z%nodes(1) + z%scale*rbg%zmean_rel
        if (tbg%relative) tbg%zmean = z%nodes(1) + z%scale*tbg%zmean_rel
        if (hbg%relative) hbg%zmean = z%nodes(1) + z%scale*hbg%zmean_rel
        do is = 1, size(sbg)
            if (sbg(is)%relative) sbg(is)%zmean = z%nodes(1) + z%scale*sbg(is)%zmean_rel
        end do

        ! #######################################################################
        ! -----------------------------------------------------------------------
        ! Construct reference  profiles
        allocate (sbackground(z%size, inb_scal_array))   ! scalar profiles

        do is = 1, inb_scal
            do k = 1, z%size
                sbackground(k, is) = Profiles_Calculate(sbg(is), z%nodes(k))
            end do
        end do

        if (imode_thermo == THERMO_TYPE_ANELASTIC) then     ! thermodynamic profiles
            allocate (epbackground(z%size))
            allocate (tbackground(z%size))
            allocate (pbackground(z%size))
            allocate (rbackground(z%size))
            allocate (ribackground(z%size))

            select case (imixture)
            case (MIXT_TYPE_AIRWATER)
           call Gravity_Hydrostatic_Enthalpy(fdm_Int0, z%nodes(:), sbackground, epbackground, tbackground, pbackground, pbg%zmean, pbg%mean, equilibrium=.true.)

            case default
                call Gravity_Hydrostatic_Enthalpy(fdm_Int0, z%nodes(:), sbackground, epbackground, tbackground, pbackground, pbg%zmean, pbg%mean)

            end select

            call Thermo_Anelastic_Rho(1, 1, z%size, sbackground, rbackground, wrk1d)
            ribackground = 1.0_wp/rbackground

        end if

        if (any(gravityProps%active(1:3))) then
            allocate (bbackground(z%size))                   ! buoyancy profiles
            bbackground(:) = 0.0_wp

            if (gravityProps%active(3)) then
                wrk1d(1:z%size, 1) = 0.0_wp
                call Gravity_AddSource(gravityProps, 1, 1, z%size, sbackground(:, 1), wrk1d, 1.0_wp)
                bbackground(1:z%size) = wrk1d(1:z%size, 1)
            end if

        end if

        ! -----------------------------------------------------------------------
        ! Add diagnostic fields to reference profile data, if any
        do is = inb_scal + 1, inb_scal_array ! Add diagnostic fields, if any
            sbg(is) = sbg(1)
            schmidt(is) = schmidt(1)
        end do

        return
    end subroutine TLab_Initialize_Background

end module Tlab_Background
