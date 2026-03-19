#include "tlab_error.h"

module Buffer
    use TLab_Constants, only: wp, wi, MAX_PARS
    use TLab_Constants, only: efile, lfile
    use TLab_Constants, only: tag_flow, tag_scal, fmt_r
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Memory, only: inb_flow, inb_scal
    implicit none
    private

    public :: Buffer_Initialize
    public :: Buffer_Nudge_Flow, Buffer_Nudge_Scal

    ! -------------------------------------------------------------------
    integer, parameter :: BUFFER_TYPE_NONE = 0
    integer, parameter :: BUFFER_TYPE_NUDGE = 1
    ! integer, parameter :: BUFFER_TYPE_FILTER = 2

    integer(wi) :: bufferType = BUFFER_TYPE_NONE

    type :: buffer_dt
        ! sequence
        integer type
        integer(wi) size, offset
        real(wp), allocatable :: tau(:)         ! relaxation timescale
        real(wp), allocatable :: ref(:, :, :)   ! reference field
        integer(wi) form                        ! form of function of relaxation term
        real(wp) :: parameters(MAX_PARS)
    contains
        procedure :: readblock => readblock_dt
    end type buffer_dt

    type, extends(buffer_dt) :: bufferZ_dt
    contains
        procedure :: initialize => bufferZ_initialize_dt
        procedure :: nudge => nudgeZ_dt
    end type

    type(bufferZ_dt), allocatable :: bufferFlowKmin(:), bufferFlowKmax(:)
    type(bufferZ_dt), allocatable :: bufferScalKmin(:), bufferScalKmax(:)

    integer(wi), parameter :: FORM_POWER_MIN = 1
    integer(wi), parameter :: FORM_POWER_MAX = 2

    logical :: bufferLoad

contains
    !########################################################################
    !########################################################################
    subroutine Buffer_Initialize(inifile)
        use TLab_Arrays, only: q, s
        use TLab_Grid, only: z

        character(len=*), intent(in) :: inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block, name
        character(len=128) eStr
        character(len=512) sRes
        integer is, iq

        !########################################################################
        ! Read
        bakfile = trim(adjustl(inifile))//'.bak'
        block = 'BufferZone'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Type=<none/relaxation/filter/both>')
        call TLab_Write_ASCII(bakfile, '#LoadBuffer=<yes/no>')

        call ScanFile_Char(bakfile, inifile, block, 'Type', 'none', sRes)
        if (trim(adjustl(sRes)) == 'none') then; bufferType = BUFFER_TYPE_NONE
        else if (trim(adjustl(sRes)) == 'relaxation') then; bufferType = BUFFER_TYPE_NUDGE
            ! else if (trim(adjustl(sRes)) == 'filter') then; bufferType = BUFFER_TYPE_FILTER
            ! else if (trim(adjustl(sRes)) == 'both') then; bufferType = BUFFER_TYPE_BOTH
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong Type option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'LoadBuffer', 'no', sRes)
        if (trim(adjustl(sRes)) == 'yes') then
            bufferLoad = .true.
        else
            bufferLoad = .false.
        end if

        !########################################################################
        ! Read & Initialize
        if (bufferType /= BUFFER_TYPE_NONE) then
            allocate (bufferFlowKmin(1:inb_flow))
            bufferFlowKmin(:)%type = bufferType            ! So far, all the same
            bufferFlowKmin(:)%form = FORM_POWER_MIN
            do iq = 1, inb_flow
                call bufferFlowKmin(iq)%readblock(bakfile, inifile, block, 'UKmin')
                bufferFlowKmin(iq)%offset = 0
                write (name, *) iq
                call bufferFlowKmin(iq)%initialize(trim(adjustl(tag_flow))//'bcs.kmin.'//trim(adjustl(name)), &
                                                   q(:, iq))
            end do

            allocate (bufferFlowKmax(1:inb_flow))
            bufferFlowKmax(:)%type = bufferType            ! So far, all the same
            bufferFlowKmax(:)%form = FORM_POWER_MAX
            do iq = 1, inb_flow
                call bufferFlowKmax(iq)%readblock(bakfile, inifile, block, 'UKmax')
                bufferFlowKmax(iq)%offset = z%size - bufferFlowKmax(iq)%size
                write (name, *) iq
                call bufferFlowKmax(iq)%initialize(trim(adjustl(tag_flow))//'bcs.kmax.'//trim(adjustl(name)), &
                                                   q(:, iq))
            end do

            allocate (bufferScalKmin(1:inb_scal))
            bufferScalKmin(:)%type = bufferType             ! So far, all the same
            bufferScalKmin(:)%form = FORM_POWER_MIN
            do is = 1, inb_scal
                call bufferScalKmin(is)%readblock(bakfile, inifile, block, 'SKmin')
                bufferScalKmin(is)%offset = 0
                write (name, *) is
                call bufferScalKmin(is)%initialize(trim(adjustl(tag_scal))//'bcs.kmin.'//trim(adjustl(name)), &
                                                   s(:, is))
            end do

            allocate (bufferScalKmax(1:inb_scal))
            bufferScalKmax(:)%type = bufferType; ! So far, all the same
            bufferScalKmax(:)%form = FORM_POWER_MAX
            do is = 1, inb_scal
                call bufferScalKmax(is)%readblock(bakfile, inifile, block, 'SKmax')
                bufferScalKmax(is)%offset = z%size - bufferScalKmax(is)%size
                write (name, *) is
                call bufferScalKmax(is)%initialize(trim(adjustl(tag_scal))//'bcs.kmax.'//trim(adjustl(name)), &
                                                   s(:, is))
            end do

        end if

        return
    end subroutine Buffer_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine readblock_dt(self, bakfile, inifile, block, tag)
        class(buffer_dt), intent(inout) :: self
        character(len=*) bakfile, inifile, block, tag

        character(len=512) sRes
        character(len=128) eStr
        integer idummy

        ! ###################################################################
        call TLab_Write_ASCII(bakfile, '#Points'//trim(adjustl(tag))//'=<value>')
        call TLab_Write_ASCII(bakfile, '#Parameters'//trim(adjustl(tag))//'=<values>')

        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call ScanFile_Int(bakfile, inifile, block, 'Points'//trim(adjustl(tag)), '0', self%size)

        if (self%size > 0) then
            call ScanFile_Char(bakfile, inifile, block, 'Parameters'//trim(adjustl(tag)), '1.0, 2.0', sRes)

            idummy = 2; call LIST_REAL(sRes, idummy, self%parameters)
            if (idummy /= 2) then
                call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong number of values in buffer parameters.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

            if (self%parameters(1) == 0.0_wp) self%type = BUFFER_TYPE_NONE

        else
            self%type = BUFFER_TYPE_NONE

        end if

    end subroutine

    ! ###################################################################
    ! ###################################################################
    subroutine bufferZ_initialize_dt(self, tag, field)
        use TLab_Memory, only: imax, jmax, kmax
        use TLab_Grid, only: z
        use IO_Fields
        use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC
        use Thermo_Anelastic, only: rbackground
        use Averages, only: AVG1V2D
        use TLab_Time, only: itime
#ifdef USE_MPI
        use mpi_f08, only: MPI_COMM_WORLD, MPI_REAL8
#endif
        class(bufferZ_dt), intent(inout) :: self
        character(len=*), intent(in) :: tag         ! File name information
        real(wp), intent(in) :: field(imax, jmax, kmax)

        character(len=32) str, name
        character(len=128) line
        integer(wi) k, kglobal
        real(wp) dummy
        type(io_subarray_dt) :: io_subarrays
        integer(wi) io_sizes(5), idummy

        ! ###################################################################
        if (self%size <= 0) return

        ! ###################################################################
        ! Reference fields
        allocate (self%ref(imax, jmax, self%size))

        io_subarrays%offset = 0
        io_subarrays%precision = wp
#ifdef USE_MPI
        io_subarrays%active = .true.
        io_subarrays%communicator = MPI_COMM_WORLD
        io_subarrays%subarray = IO_Create_Subarray_XOY(imax, jmax, self%size, MPI_REAL8)
#endif
        idummy = imax*jmax*self%size; io_sizes = [idummy, 1, idummy, 1, 1]

        name = trim(adjustl(tag))

        if (bufferLoad) then
            call IO_Read_Subarray(io_subarrays, name, [' '], self%ref, io_sizes)

        else
            do k = 1, self%size
                kglobal = k + self%offset
                self%ref(:, :, k) = AVG1V2D(imax, jmax, kmax, kglobal, 1, field)
            end do
            idummy = index(name, '.', back=.true.)
            write (str, *) itime; name = trim(name(1:idummy - 1))//'.'//trim(adjustl(str))//trim(name(idummy:))
            call IO_Write_Subarray(io_subarrays, name, [' '], self%ref, io_sizes)

        end if

        ! Control
        line = 'Checking bounds of field '//trim(adjustl(tag))//':'
        write (str, fmt_r) minval(self%ref(:, :, :)); line = trim(adjustl(line))//' '//trim(adjustl(str))//','
        write (str, fmt_r) maxval(self%ref(:, :, :)); line = trim(adjustl(line))//' '//trim(adjustl(str))//'.'
        call TLab_Write_ASCII(lfile, line)

        ! ###################################################################
        ! Inverse of relaxation time
        allocate (self%tau(self%size))

#define strength parameters(1)
#define sigma parameters(2)

        dummy = 1.0_wp/(z%nodes(self%offset + self%size) - z%nodes(self%offset + 1))    ! Inverse of segment length
        do k = 1, self%size
            kglobal = self%offset + k
            if (self%form == FORM_POWER_MAX) &
                self%tau(k) = self%strength*((z%nodes(kglobal) - z%nodes(self%offset + 1))*dummy)**self%sigma

            if (self%form == FORM_POWER_MIN) &
                self%tau(k) = self%strength*((z%nodes(self%offset + self%size) - z%nodes(kglobal))*dummy)**self%sigma

            if (nse_eqns == DNS_EQNS_ANELASTIC) &                                       ! formulation per unit volume
                self%tau(k) = self%tau(k)*rbackground(kglobal)

        end do

#undef strength
#undef sigma

        return
    end subroutine bufferZ_initialize_dt

    ! ###################################################################
    ! ###################################################################
    subroutine nudgeZ_dt(self, s, hs)
        class(bufferZ_dt), intent(in) :: self
        real(wp), intent(in) :: s(:, :, :)
        real(wp), intent(out) :: hs(:, :, :)

        integer(wi) k, kglobal

        do k = 1, self%size
            kglobal = k + self%offset
            hs(:, :, kglobal) = hs(:, :, kglobal) - self%tau(k)*(s(:, :, kglobal) - self%ref(:, :, k))
        end do

        return
    end subroutine nudgeZ_dt

    ! ###################################################################
    ! ###################################################################
    subroutine Buffer_Nudge_Flow(q, hq)
        use TLab_Memory, only: imax, jmax, kmax, inb_flow
        real(wp), intent(inout) :: q(imax, jmax, kmax, inb_flow)
        real(wp), intent(inout) :: hq(imax, jmax, kmax, inb_flow)

        integer iq

        if (bufferType /= BUFFER_TYPE_NUDGE) return

        do iq = 1, inb_flow
            if (bufferFlowKmin(iq)%type == BUFFER_TYPE_NUDGE) call bufferFlowKmin(iq)%nudge(q(:, :, :, iq), hq(:, :, :, iq))
            if (bufferFlowKmax(iq)%type == BUFFER_TYPE_NUDGE) call bufferFlowKmax(iq)%nudge(q(:, :, :, iq), hq(:, :, :, iq))
        end do

        return
    end subroutine

    ! ###################################################################
    ! ###################################################################
    subroutine Buffer_Nudge_Scal(s, hs)
        use TLab_Memory, only: imax, jmax, kmax, inb_scal
        real(wp), intent(inout) :: s(imax, jmax, kmax, inb_scal)
        real(wp), intent(inout) :: hs(imax, jmax, kmax, inb_scal)

        integer is

        if (bufferType /= BUFFER_TYPE_NUDGE) return

        do is = 1, inb_scal
            if (bufferScalKmin(is)%type == BUFFER_TYPE_NUDGE) call bufferScalKmin(is)%nudge(s(:, :, :, is), hs(:, :, :, is))
            if (bufferScalKmax(is)%type == BUFFER_TYPE_NUDGE) call bufferScalKmax(is)%nudge(s(:, :, :, is), hs(:, :, :, is))
        end do

        return
    end subroutine

end module Buffer
