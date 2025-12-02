#include "tlab_error.h"

module DNS_Arrays
    use TLab_Constants, only: wp
    implicit none

    real(wp), allocatable, target :: hq(:, :)       ! Right-hand sides Eulerian fields
    real(wp), allocatable, target :: hs(:, :)       ! Right-hand sides Eulerian fields
    ! real(wp), allocatable, target :: l_hq(:, :)     ! Right-hand sides Lagrangian fields

    real(wp), pointer :: p_hq(:, :, :, :) => null()
    real(wp), pointer :: p_hs(:, :, :, :) => null()
    real(wp), pointer :: pxy_hq(:, :, :) => null()
    real(wp), pointer :: pxy_hs(:, :, :) => null()

end module DNS_Arrays

! ###################################################################
! ###################################################################
module DNS_LOCAL
    use TLab_Constants, only: wp, wi, sp
    ! use TLab_Constants, only: MAX_PATH_LENGTH
    implicit none

    integer :: nitera_first     ! First iteration in current run
    integer :: nitera_last      ! Last iteration in current run
    integer :: nitera_save      ! Iteration step to check-point: save restart files
    integer :: nitera_stats     ! Iteration step to check-point: save statistical data
    integer :: nitera_pln       ! Iteration step to save planes
    integer :: nitera_filter    ! Iteration step for domain filter, if any
    integer :: nitera_log       ! Iteration step for data logger with simulation information

    real(wp) :: nruntime_sec    ! Maximum runtime of the simulation in seconds
    real(wp) :: wall_time       ! Actual elapsed time during the simulation in seconds
    integer :: start_clock      ! Starting time of the simulation on the system

    ! Variable viscosity
    logical :: flag_viscosity
    real(wp) :: visc_stop, visc_time, visc_rate

contains
    ! ###################################################################
    ! ###################################################################
    subroutine DNS_Initialize_Parameters(inifile)
        use TLab_Constants, only: wp, wi, big_wp, efile, lfile, wfile
        use TLab_Memory, only: inb_txc
        use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
        ! use PLANES
        ! use OPR_Filters, only: FilterDomain, PressureFilter, DNS_FILTER_NONE

        character*(*) inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block
        character(len=128) eStr
        ! character(len=512) sRes

        ! ###################################################################
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Time'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Start=<integral start time>')
        call TLab_Write_ASCII(bakfile, '#End=<integral stop time>')
        call TLab_Write_ASCII(bakfile, '#Restart=<restart time step>')
        call TLab_Write_ASCII(bakfile, '#Statistics=<statistics time step>')
        call TLab_Write_ASCII(bakfile, '#Saveplanes=<value>')
        call TLab_Write_ASCII(bakfile, '#RunAvera=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#Runtime=<seconds>')
        call TLab_Write_ASCII(bakfile, '#Logs=<value>')

        call ScanFile_Int(bakfile, inifile, block, 'Start', '0', nitera_first)
        call ScanFile_Int(bakfile, inifile, block, 'End', '0', nitera_last)
        call ScanFile_Int(bakfile, inifile, block, 'Restart', '50', nitera_save)
        call ScanFile_Int(bakfile, inifile, block, 'Statistics', '50', nitera_stats)
        call ScanFile_Int(bakfile, inifile, block, 'Saveplanes', '-1', nitera_pln)
        call ScanFile_Int(bakfile, inifile, block, 'Logs', '10', nitera_log)
        call ScanFile_Real(bakfile, inifile, block, 'Runtime', '10000000', nruntime_sec)

        ! ! Domain Filter (Should we move it to Iteration?)
        ! call ScanFile_Int(bakfile, inifile, 'Filter', 'Step', '0', nitera_filter)
        ! if (nitera_filter == 0) FilterDomain(:)%type = DNS_FILTER_NONE

        if (nitera_first > nitera_last) then
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Not started because nitera_first > nitera_last.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        ! ###################################################################
        ! Viscosity Control
        ! ###################################################################
        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#[ViscChange]')
        call TLab_Write_ASCII(bakfile, '#Time=<time>')

        call ScanFile_Real(bakfile, inifile, 'ViscChange', 'Time', '0.0', visc_time)

        ! ! ###################################################################
        ! ! Save planes to disk
        ! ! ###################################################################
        !         call TLab_Write_ASCII(bakfile, '#')
        !         call TLab_Write_ASCII(bakfile, '#[SavePlanes]')

        !         call PLANES_READBLOCK(bakfile, inifile, 'SavePlanes', 'PlanesI', iplanes)
        !         call PLANES_READBLOCK(bakfile, inifile, 'SavePlanes', 'PlanesJ', jplanes)
        !         call PLANES_READBLOCK(bakfile, inifile, 'SavePlanes', 'PlanesK', kplanes)

        ! ###################################################################
        ! Final initialization and consistency check
        ! ###################################################################
        ! Avoid dividing by zero in time_integration routine
        if (nitera_save <= 0) nitera_save = nitera_last - nitera_first + 1
        if (nitera_stats <= 0) nitera_stats = nitera_last - nitera_first + 1
        if (nitera_log <= 0) nitera_log = nitera_last - nitera_first + 1
        if (nitera_pln <= 0) nitera_pln = nitera_last - nitera_first + 1
        if (nitera_filter <= 0) nitera_filter = nitera_last - nitera_first + 1

        ! -------------------------------------------------------------------
        ! Array sizes
        ! -------------------------------------------------------------------
        inb_txc = 7

        return
    end subroutine DNS_Initialize_Parameters

end module DNS_LOCAL
