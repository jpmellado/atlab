#include "tlab_error.h"

module TLab_WorkFlow
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: lfile, efile, fmt_r
! #ifdef USE_OPENMP
!     use OMP_LIB
! #endif
#ifdef USE_MPI
    use mpi_f08
#endif
    implicit none
    private

    public :: TLab_Start
    public :: TLab_Stop
    public :: TLab_Write_ASCII
    public :: TLab_Runtime

#ifdef USE_MPI
    integer, protected, public :: num_processors
    integer, protected, public :: rank
    real(wp), public :: runtime_mpi             ! runtime accumulated in mpi part
#endif

    integer, public :: imode_verbosity = 1      ! level of verbosity used in log files

    logical, public :: flow_on = .true.         ! calculate flow parts of the code
    logical, public :: scal_on = .true.         ! calculate scal parts of the code
    ! logical, public :: stagger_on = .false.     ! horizontal staggering of pressure

    ! -------------------------------------------------------------------
#ifdef USE_MPI
    real(wp) :: walltime_ref                    ! reference time in seconds
    integer :: ierr
#else
    integer :: clock_ref                        ! reference time in clock cycles
#endif

    character(len=128) :: line

contains

    ! ###################################################################
    ! ###################################################################
    subroutine TLab_Start()
        character(len=10) clock(2)

        !#####################################################################
        ! Initialize MPI parallel mode
#ifdef USE_MPI
        call MPI_INIT(ierr)

        call MPI_COMM_SIZE(MPI_COMM_WORLD, num_processors, ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

        write (line, *) num_processors
        line = 'Number of MPI tasks '//trim(adjustl(line))
        call TLab_Write_ASCII(lfile, line)

        if (num_processors == 0) then
            call TLab_Write_ASCII(efile, __FILE__//'. Number of processors is zero.')
            call TLab_Stop(DNS_ERROR_MINPROC)
        end if

#endif

        !#####################################################################
        ! Initialize calculation of running time
#ifdef USE_MPI
        walltime_ref = MPI_WTIME()
#ifdef PROFILE_ON
        runtime_mpi = 0.0_wp
#endif
#else
        call system_clock(clock_ref)
#endif

        !########################################################################
        ! First output
        call date_and_time(clock(1), clock(2))
        line = 'Starting on '//trim(adjustl(clock(1) (1:8)))//' at '//trim(adjustl(clock(2)))
        call TLab_Write_ASCII(lfile, line)

        line = 'Git-hash '//GITHASH
        call TLab_Write_ASCII(lfile, line)

        line = 'Git-branch '//GITBRANCH
        call TLab_Write_ASCII(lfile, line)

        return
    end subroutine TLab_Start

    ! ###################################################################
    ! ###################################################################
    subroutine TLab_Stop(error_code)
        integer, intent(in) :: error_code

        real(wp) runtime

        ! ###################################################################
        if (error_code /= 0) then
            write (line, *) error_code
            line = 'Error code '//trim(adjustl(line))//'.'
            call TLab_Write_ASCII(efile, line)
        end if

        call GETARG(0, line)
        write (line, *) 'Finalizing program '//trim(adjustl(line))
        if (error_code == 0) then
            line = trim(adjustl(line))//' normally.'
        else
            line = trim(adjustl(line))//' abnormally. Check '//trim(adjustl(efile))
        end if
        call TLab_Write_ASCII(lfile, line)

        !#####################################################################
        runtime = TLab_Runtime()
        write (line, fmt_r) runtime
        line = 'Time elapse ....................: '//trim(adjustl(line))
        call TLab_Write_ASCII(lfile, line)

#ifdef USE_MPI
#ifdef PROFILE_ON
        write (line, fmt_r) runtime_mpi
        line = 'Time in mpi transposition ......: '//trim(ADJUST(line))
        call TLab_Write_ASCII(lfile, line)
#endif
#endif
        !#####################################################################
        call TLab_Write_ASCII(lfile, '########################################')

#ifdef USE_MPI
        if (error_code == 0) then
            call MPI_FINALIZE(ierr)
        else
            call MPI_Abort(MPI_COMM_WORLD, error_code, ierr)
        end if
#else
        stop
#endif

        return
    end subroutine TLab_Stop

    ! ###################################################################
    ! ###################################################################
    ! Execution elapsed time in seconds
    function TLab_Runtime() result(runtime)
        real(wp) :: runtime

#ifdef USE_MPI
#else
        integer clock_loc, int_dummy
#endif

#ifdef USE_MPI
        runtime = MPI_WTIME() - walltime_ref
        call MPI_BCast(runtime, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
#else
        ! call ETIME(tdummy, clock_loc)
        call system_clock(clock_loc, int_dummy)
        runtime = real(clock_loc - clock_ref)/int_dummy
#endif

        return
    end function

    ! ###################################################################
    ! ###################################################################
    subroutine TLab_Write_ASCII(file, locLine, flag_all)
        character*(*), intent(in) :: file, locLine
        logical, intent(in), optional :: flag_all

        ! -----------------------------------------------------------------------
        character*10 clock(2)

        ! #######################################################################
#ifdef USE_MPI
        if (rank == 0 .or. present(flag_all)) then
#endif

            if (imode_verbosity > 0) then

                open (UNIT=22, FILE=file, STATUS='unknown', POSITION='APPEND')
                if (imode_verbosity == 1) then
                    write (22, '(a)') trim(adjustl(locLine))
                else if (imode_verbosity == 2) then
                    call date_and_time(clock(1), clock(2))
                    write (22, '(a)') '['//trim(adjustr(clock(2)))//'] '//trim(adjustl(locLine))
                end if
                close (22)

            end if

#ifndef PARALLEL
            if (file == efile) then
                write (*, *) trim(adjustl(locLine))
            end if
#endif

#ifdef USE_MPI
        end if
#endif

        return
    end subroutine TLab_Write_ASCII

end module TLab_WorkFlow
