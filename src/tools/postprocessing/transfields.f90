#include "tlab_error.h"

program TransFields
    use TLab_Constants, only: wp, wi, MAX_PARS
    use TLab_Constants, only: ifile, lfile, efile, tag_flow, tag_scal
    use TLab_Memory, only: imax, jmax, kmax, inb_txc, inb_flow, inb_scal
    use TLab_Time, only: itime, rtime
    use TLab_Arrays
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start, flow_on, scal_on
    use TLab_Memory, only: TLab_Initialize_Memory
#ifdef USE_MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use IO_Fields
    use TLab_Grid
    use NavierStokes
    use Thermodynamics, only: Thermo_Initialize
    use Gravity, only: Gravity_Initialize
    use Rotation, only: Rotation_Initialize
    implicit none

    integer(wi), parameter :: itime_size_max = 3000 ! iterations to be processed
    integer(wi) itime_size, it
    integer(wi) itime_vec(itime_size_max)

    integer(wi), parameter :: iopt_size_max = 30    ! options to be processed
    integer(wi) iopt_size, iv
    integer(wi) opt_vec(iopt_size_max)
    character(len=64) opt_name(iopt_size_max)

    character*32 fname
    character(len=32) time_str                      ! Time stamp

    integer is, iq
    logical iread_flow, iread_scal
    real(wp) params(MAX_PARS)

    integer opt_main

    !########################################################################
    !########################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
    call IO_Initialize()

    call TLab_Grid_Initialize()

#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Trp_Initialize(ifile)
#endif

    call NavierStokes_Initialize_Parameters(ifile)
    call Thermo_Initialize(ifile)

    call Gravity_Initialize(ifile)
    call Rotation_Initialize(ifile)

    call TLab_Consistency_Check()

    call TransFields_Initialize()

    ! #######################################################################
    call TLab_Initialize_Memory(__FILE__)

    ! ###################################################################
    ! Postprocess given list of files
    ! ###################################################################
    do it = 1, itime_size
        itime = itime_vec(it)

        write (time_str, *) itime

        call TLab_Write_ASCII(lfile, 'Processing iteration It'//trim(adjustl(time_str)))

        if (iread_scal) then
            fname = trim(adjustl(tag_scal))//trim(adjustl(time_str))
            call IO_Read_Fields(fname, imax, jmax, kmax, itime, inb_scal, 0, s, params(1:1))
            rtime = params(1)
        end if

        if (iread_flow) then
            fname = trim(adjustl(tag_flow))//trim(adjustl(time_str))
            call IO_Read_Fields(fname, imax, jmax, kmax, itime, inb_flow, 0, q, params(1:1))
            rtime = params(1)
        end if

        iv = 1
        select case (trim(adjustl(opt_name(opt_vec(iv)))))
        case ('Convert double to single precision')
            call DoubleToSingle()

        case ('Swap Y-Z coordinates')
            do iq = 1, inb_flow
                call SwapYZ(q(:, iq))
            end do
            do is = 1, inb_scal
                call SwapYZ(s(:, is))
            end do

        end select

        if (flow_on) then
            write (fname, *) itime; fname = trim(adjustl(tag_flow))//'trn.'//trim(adjustl(fname))
            io_header_q(1)%params(1) = rtime
            call IO_Write_Fields(fname, imax, jmax, kmax, itime, inb_flow, q, io_header_q(1:1))
        end if

        if (scal_on) then
            write (fname, *) itime; fname = trim(adjustl(tag_scal))//'trn.'//trim(adjustl(fname))
            io_header_s(:)%params(1) = rtime
            call IO_Write_Fields(fname, imax, jmax, kmax, itime, inb_scal, s, io_header_s(1:inb_scal))
        end if

        select case (trim(adjustl(opt_name(opt_vec(iv)))))
        case ('Convert double to single precision')
            call DoubleToSingle(datatype="original")

        end select

    end do

    call TLab_Stop(0)

contains
    !########################################################################
    !########################################################################
    subroutine TransFields_Initialize()

        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        ! #######################################################################
        ! Read from tlab.ini
        bakfile = trim(adjustl(ifile))//'.bak'
        block = 'PostProcessing'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        ! -------------------------------------------------------------------
        call ScanFile_Char(bakfile, ifile, block, 'Files', '-1', sRes)
        if (sRes == '-1') then
#ifdef USE_MPI
#else
            write (*, *) 'Iteration numbers ?'
            read (*, '(A512)') sRes
#endif
        end if
        itime_size = itime_size_max
        call LIST_INTEGER(sRes, itime_size, itime_vec)

        if (itime_vec(1) < 0) then ! Check
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Missing Files input.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
        end if

        ! -------------------------------------------------------------------
        iv = 0
        iv = iv + 1; opt_name(iv) = 'Convert double to single precision'
        iv = iv + 1; opt_name(iv) = 'Swap Y-Z coordinates' ! old to new file format

        if (iv > iopt_size_max) then ! Check
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Increase number of options.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
        end if

        call ScanFile_Char(bakfile, ifile, block, 'ParamTransFields', '-1', sRes)
        if (sRes == '-1') then
#ifdef USE_MPI
#else
            write (*, '(A)') 'Option?'
            do is = 1, iv
                write (*, '(I2,A)') is, '. '//trim(adjustl(opt_name(is)))
            end do
            read (*, '(A512)') sRes
#endif
        end if
        opt_vec(:) = -1
        iopt_size = iopt_size_max
        call LIST_INTEGER(sRes, iopt_size, opt_vec)

        if (opt_vec(1) < 0) then ! Check
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Missing input [ParamAverages] in tlab.ini.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
        end if

        opt_main = opt_vec(1)

        ! ###################################################################
        ! Initialization of array sizes
        ! ###################################################################
        iread_flow = .false.
        iread_scal = .false.
        inb_txc = 0

        iv = 1
        select case (trim(adjustl(opt_name(opt_vec(iv)))))
        case ('Convert double to single precision')
            iread_flow = flow_on; iread_scal = scal_on
        case ('Swap Y-Z coordinates')
            iread_flow = flow_on; iread_scal = scal_on
        end select

        return
    end subroutine TransFields_Initialize

    !########################################################################
    !########################################################################
    subroutine DoubleToSingle(datatype)
        character(len=*), optional :: datatype
        integer, save :: original_datatype

        if (present(datatype)) then
            if (trim(adjustl(datatype)) == "original") then
                io_datatype = original_datatype
                io_subarray_main%precision = original_datatype
            end if
        else
            original_datatype = io_datatype
            io_datatype = IO_TYPE_SINGLE
            io_subarray_main%precision = io_datatype
        end if

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine SwapYZ(field)
        real(wp), intent(inout) :: field(:)

        integer j, k
        target field
        real(wp), pointer :: p_org(:, :, :) => null(), p_dst(:, :, :) => null()

        p_org(1:imax, 1:kmax, 1:jmax) => field(1:imax*jmax*kmax)
        p_dst(1:imax, 1:jmax, 1:kmax) => wrk3d(1:imax*jmax*kmax)

        do k = 1, kmax
            do j = 1, jmax
                p_dst(1:imax, j, k) = p_org(1:imax, k, j)
            end do
        end do
        field(1:imax*jmax*kmax) = wrk3d(1:imax*jmax*kmax)

        nullify (p_org, p_dst)
        
        return
    end subroutine

end program TransFields
