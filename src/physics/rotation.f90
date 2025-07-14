#include "tlab_error.h"

module Rotation
    use TLab_Constants, only: wp, wi, pi_wp, efile, lfile, MAX_VARS, MAX_PARS
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private

    public :: Rotation_Initialize
    public :: Rotation_AddCoriolis

    real(wp), public, protected :: rossby

    type coriolis_dt
        sequence
        integer type
        ! logical active(MAX_VARS), lpadding(3)   ! fields affected by this term
        real(wp) parameters(MAX_PARS)
        real(wp) vector(3)
    end type coriolis_dt
    type(coriolis_dt), public :: coriolisProps

    ! -------------------------------------------------------------------
    integer, parameter :: TYPE_COR_NONE = 0
    integer, parameter :: TYPE_COR_AGEOSTROPHIC_Z = 1
    integer, parameter :: TYPE_COR_EXPLICIT_3D = 2

contains
    !########################################################################
    !########################################################################
    subroutine Rotation_Initialize(inifile)
        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes
        integer(wi) idummy
        real(wp) dummy

        !########################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Coriolis'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<none/ageostrophic/explicit>')
        call TLab_Write_ASCII(bakfile, '#Rossby=<value>')
        call TLab_Write_ASCII(bakfile, '#Vector=<Fx,Fy,Fz>')
        call TLab_Write_ASCII(bakfile, '#Parameters=<value>')

        ! Coriolis
        call ScanFile_Char(bakfile, inifile, block, 'Type', 'None', sRes)
        if (trim(adjustl(sRes)) == 'none') then; coriolisProps%type = TYPE_COR_NONE
        else if (trim(adjustl(sRes)) == 'ageostrophic') then; coriolisProps%type = TYPE_COR_AGEOSTROPHIC_Z
        else if (trim(adjustl(sRes)) == 'explicit') then; coriolisProps%type = TYPE_COR_EXPLICIT_3D
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Error in entry Type.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        coriolisProps%vector(:) = 0.0_wp
        call ScanFile_Char(bakfile, inifile, block, 'Vector', '0.0,0.0,0.0', sRes)
        idummy = 3
        call LIST_REAL(sRes, idummy, coriolisProps%vector)

        call ScanFile_Real(bakfile, inifile, block, 'Rossby', '-1.0', rossby)
        if (rossby <= 0.0) then
            call ScanFile_Real(bakfile, inifile, block, 'Coriolis', '1.0', dummy)   ! default value
            rossby = 1.0_wp/dummy
        end if

        ! coriolisProps%active = .false.
        ! if (abs(coriolisProps%vector(1)) > 0.0_wp) then; coriolisProps%active(2) = .true.; coriolisProps%active(3) = .true.; call TLab_Write_ASCII(lfile, 'Angular velocity along Ox.'); end if
        ! if (abs(coriolisProps%vector(2)) > 0.0_wp) then; coriolisProps%active(3) = .true.; coriolisProps%active(1) = .true.; call TLab_Write_ASCII(lfile, 'Angular velocity along Oy.'); end if
        ! if (abs(coriolisProps%vector(3)) > 0.0_wp) then; coriolisProps%active(1) = .true.; coriolisProps%active(2) = .true.; call TLab_Write_ASCII(lfile, 'Angular velocity along Oz.'); end if

        if (rossby > 0.0_wp) then
            coriolisProps%vector(:) = coriolisProps%vector(:)/rossby ! adding the rossby number into the vector
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Rossby number must be nonzero.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        if (coriolisProps%type == TYPE_COR_AGEOSTROPHIC_Z) then
            coriolisProps%parameters(:) = 0.0_wp
            call ScanFile_Char(bakfile, inifile, block, 'Parameters', '1.0, 0.0', sRes) ! Magnitude and direction of geostrophic wind
            idummy = MAX_PARS
            call LIST_REAL(sRes, idummy, coriolisProps%parameters)

            ! if (coriolisProps%active(3)) then       ! Consistency check
            !     call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Ageostrophic formulation only for angular velocity along Oz.')
            !     call TLab_Stop(DNS_ERROR_OPTION)
            ! end if

        end if

        return
    end subroutine Rotation_Initialize

    !########################################################################
    !########################################################################
    ! Remember that coriolisProps%vector already contains the Rossby #.
    subroutine Rotation_AddCoriolis(locProps, nx, ny, nz, u, r)
        type(coriolis_dt), intent(in) :: locProps
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz, *)
        real(wp), intent(out) :: r(nx*ny*nz, *)

        ! -----------------------------------------------------------------------
        integer(wi) i
        real(wp) fx, fy, fz, geo_u, geo_v

        ! #######################################################################
        select case (locProps%type)
        case (TYPE_COR_EXPLICIT_3D)
            fx = locProps%vector(1)
            fy = locProps%vector(2)
            fz = locProps%vector(3)
            do i = 1, nx*ny*nz
                r(i, 1) = r(i, 1) + fz*u(i, 2) - fy*u(i, 3)
                r(i, 2) = r(i, 2) + fx*u(i, 3) - fz*u(i, 1)
                r(i, 3) = r(i, 3) + fy*u(i, 1) - fx*u(i, 2)
            end do

        case (TYPE_COR_AGEOSTROPHIC_Z)
            geo_u = cos(locProps%parameters(2)*pi_wp/180.0_wp)*locProps%parameters(1)
            geo_v = sin(locProps%parameters(2)*pi_wp/180.0_wp)*locProps%parameters(1)
            fy = locProps%vector(3)
            do i = 1, nx*ny*nz
                r(i, 1) = r(i, 1) + fy*(u(i, 2) - geo_v)
                r(i, 2) = r(i, 2) + fy*(geo_u - u(i, 1))
            end do
        end select

    end subroutine Rotation_AddCoriolis

end module Rotation
