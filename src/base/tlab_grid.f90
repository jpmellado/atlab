#include "tlab_error.h"

module TLab_Grid
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none
    private

    public :: TLab_Grid_Initialize
    public :: TLab_Grid_Read
    public :: TLab_Grid_Write

    public :: axis_dt
    public :: globalGrid, x, y, z
    public :: subaxis_dt
    public :: xSubgrid, ySubgrid, zSubgrid

    ! -----------------------------------------------------------------------
    type :: axis_dt
        ! sequence
        character*8 name
        integer(wi) size
        logical :: uniform = .false.
        logical :: periodic = .false.
        real(wp) scale
        real(wp), allocatable :: nodes(:)
    end type

    type :: grid_dt
        type(axis_dt) axes(3)
    end type

    type(grid_dt), target :: globalGrid
    type(axis_dt), pointer :: x => globalGrid%axes(1)
    type(axis_dt), pointer :: y => globalGrid%axes(2)
    type(axis_dt), pointer :: z => globalGrid%axes(3)

    type :: subaxis_dt
        type(axis_dt), pointer :: parent
        integer :: offset
        integer :: size
    contains
        procedure :: initialize => subaxis_initialize
    end type

    type :: subgrid_dt
        type(subaxis_dt) axes(3)
    end type

    type(subgrid_dt), target :: subgrid
    type(subaxis_dt), pointer :: xSubgrid => subgrid%axes(1)
    type(subaxis_dt), pointer :: ySubgrid => subgrid%axes(2)
    type(subaxis_dt), pointer :: zSubgrid => subgrid%axes(3)

contains
    !########################################################################
    !########################################################################
    subroutine subaxis_initialize(self, referenceAxis)
        class(subaxis_dt), intent(inout) :: self
        type(axis_dt), intent(in), target :: referenceAxis

        self%parent => referenceAxis
        self%offset = 0
        self%size = referenceAxis%size

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine TLab_Grid_Initialize(locInifile)
        use TLab_Constants, only: ifile, gfile
        use TLab_Memory, only: imax, jmax, kmax
        character(len=*), optional, intent(in) :: locInifile

        ! -------------------------------------------------------------------
        character(len=32) inifile
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        integer ig

        ! ###################################################################
        if (present(locInifile)) then
            inifile = locInifile
        else
            inifile = ifile
        end if

        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'Grid'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Imax=<value>')
        call TLab_Write_ASCII(bakfile, '#XUniform=<yes/no>')
        call TLab_Write_ASCII(bakfile, '#XPeriodic=<yes/no>')
        ! and same fo Y and Z

        call ScanFile_Int(bakfile, inifile, block, 'Imax', '0', imax)
        call ScanFile_Int(bakfile, inifile, block, 'Jmax', '0', jmax)
        call ScanFile_Int(bakfile, inifile, block, 'Kmax', '0', kmax)

        globalGrid%axes(1)%name = 'x'
        globalGrid%axes(2)%name = 'y'
        globalGrid%axes(3)%name = 'z'

        do ig = 1, 3
            call ScanFile_Char(bakfile, inifile, block, globalGrid%axes(ig)%name(1:1)//'Uniform', 'void', sRes)
            if (trim(adjustl(sRes)) == 'yes') then
                globalGrid%axes(ig)%uniform = .true.
            else if (trim(adjustl(sRes)) == 'no') then
                globalGrid%axes(ig)%uniform = .false.
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Error in Uniform '//globalGrid%axes(ig)%name(1:1)//' grid')
                call TLab_Stop(DNS_ERROR_UNIFORMX)
            end if

            call ScanFile_Char(bakfile, inifile, block, globalGrid%axes(ig)%name(1:1)//'Periodic', 'void', sRes)
            if (trim(adjustl(sRes)) == 'yes') then
                globalGrid%axes(ig)%periodic = .true.
            else if (trim(adjustl(sRes)) == 'no') then
                globalGrid%axes(ig)%periodic = .false.
            else
                call TLab_Write_ASCII(efile, __FILE__//'. Error in Periodic '//globalGrid%axes(ig)%name(1:1)//' grid')
                call TLab_Stop(DNS_ERROR_IBC)
            end if

            ! consistency check
            if (globalGrid%axes(ig)%periodic .and. (.not. globalGrid%axes(ig)%uniform)) then
                call TLab_Write_ASCII(efile, __FILE__//'. Grid must be uniform in periodic direction.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

        end do

        ! ###################################################################
        call TLab_Grid_Read(gfile, x, y, z, [imax, jmax, kmax])

        ! ###################################################################
        ! Default local subgrid is equal to global grid
        do ig = 1, 3
            call subgrid%axes(ig)%initialize(globalGrid%axes(ig))
        end do

#ifdef USE_MPI
        ! MPI domain decomposition in X and Y directions
        call TLab_Write_ASCII(bakfile, '#Imax(*)=<value>')
        call TLab_Write_ASCII(bakfile, '#Jmax(*)=<value>')

        write (sRes, *) x%size
        call ScanFile_Int(bakfile, inifile, block, 'Imax(*)', trim(adjustl(sRes)), subgrid%axes(1)%size)
        write (sRes, *) y%size
        call ScanFile_Int(bakfile, inifile, block, 'Jmax(*)', trim(adjustl(sRes)), subgrid%axes(2)%size)

#endif
        ! ###################################################################
        ! For readability in the code
        imax = xSubgrid%size
        jmax = ySubgrid%size
        kmax = zSubgrid%size

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine TLab_Grid_Read(name, x, y, z, sizes)
        character*(*) name
        type(axis_dt), intent(inout) :: x, y, z
        integer(wi), intent(in), optional :: sizes(3)

        ! -----------------------------------------------------------------------
        character*(32) line

        ! #######################################################################
        open (50, file=name, status='old', form='unformatted')
        rewind (50)

        ! -----------------------------------------------------------------------
        read (50) x%size, y%size, z%size

        if (present(sizes)) then        ! check
            if (any([x%size, y%size, z%size] /= sizes)) then
                close (50)
                write (line, 100) x%size, y%size, z%size
                call TLab_Write_ASCII(efile, __FILE__//'. Dimensions ('//trim(line)//') unmatched.')
                call TLab_Stop(DNS_ERROR_DIMGRID)
            end if
        end if

        read (50) x%scale, y%scale, z%scale

        if (allocated(x%nodes)) deallocate (x%nodes)
        if (allocated(y%nodes)) deallocate (y%nodes)
        if (allocated(z%nodes)) deallocate (z%nodes)
        allocate (x%nodes(x%size), y%nodes(y%size), z%nodes(z%size))

        read (50) x%nodes(:)
        read (50) y%nodes(:)
        read (50) z%nodes(:)

        ! -----------------------------------------------------------------------
        close (50)

        return

100     format(I5, ',', I5, ',', I5)

    end subroutine TLab_Grid_Read

!########################################################################
!########################################################################
    subroutine TLab_Grid_Write(name, x, y, z)
        character*(*) name
        type(axis_dt), intent(in) :: x, y, z

        !########################################################################
        open (unit=51, file=name, form='unformatted', status='unknown')

        write (51) x%size, y%size, z%size
        write (51) x%scale, y%scale, z%scale

        write (51) x%nodes(1:x%size)
        write (51) y%nodes(1:y%size)
        write (51) z%nodes(1:z%size)

        close (51)

        return
    end subroutine TLab_Grid_Write

end module TLab_Grid
