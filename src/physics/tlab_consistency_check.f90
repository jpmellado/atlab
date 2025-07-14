#include "tlab_error.h"

! Everything has been read from input file.
! Check for cross dependencies and undeveloped options.
subroutine TLab_Consistency_Check()
    use TLab_Constants, only: wi, wp, efile, lfile, MAX_VARS
    use TLab_Memory, only: inb_flow, inb_scal
    use TLab_Memory, only: imax, jmax, kmax
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_j, ims_npro_k
#endif
    use IO_Fields
    use TLab_Time, only: rtime
    use TLab_Grid, only: x, y, z
    use NavierStokes
    use Thermodynamics
    use Thermo_Base, only: gamma0
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    implicit none

    ! -------------------------------------------------------------------
    integer is, ig
    integer(wi) grid_sizes(3), grid_sizes_aux(3)
    character(len=32) lstr

    ! ###################################################################
#ifdef USE_MPI
    grid_sizes = [imax*ims_npro_i, jmax*ims_npro_j, kmax*ims_npro_k]
#else
    grid_sizes = [imax, jmax, kmax]
#endif
    grid_sizes_aux(1:3) = [x%size, y%size, z%size]
    do ig = 1, 3
        if (grid_sizes_aux(ig) /= grid_sizes(ig)) then
            write (lstr, *) ig
            call TLab_Write_ASCII(efile, __FILE__//'. Grid size mismatch along direction'//trim(adjustl(lstr))//'.')
            call TLab_Stop(DNS_ERROR_DIMGRID)
        end if
    end do

    ! ###################################################################
    if (max(inb_flow, inb_scal) > MAX_VARS) then
        call TLab_Write_ASCII(efile, __FILE__//'. Error MAX_VARS should be larger than or equal to inb_flow and inb_scal')
        call TLab_Stop(DNS_ERROR_TOTALVARS)
    end if

    ! ! ###################################################################
    ! if (any([EQNS_TRANS_SUTHERLAND, EQNS_TRANS_POWERLAW] == itransport)) inb_flow_array = inb_flow_array + 1    ! space for viscosity

    ! ###################################################################
    if (imode_thermo == THERMO_TYPE_COMPRESSIBLE .and. nse_eqns /= DNS_EQNS_COMPRESSIBLE) then
        call TLab_Write_ASCII(efile, __FILE__//'. Incorrect combination of compressible thermodynamics and type of evolution equations.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

    if (imode_thermo == THERMO_TYPE_ANELASTIC .and. nse_eqns /= DNS_EQNS_ANELASTIC) then
        call TLab_Write_ASCII(efile, __FILE__//'. Incorrect combination of anelastic thermodynamics and type of evolution equations.')
        call TLab_Stop(DNS_ERROR_OPTION)
    end if

    ! ###################################################################
    ! preparing headers of restart files
    io_header_q(1)%size = 0
    io_header_q(1)%size = io_header_q(1)%size + 1; io_header_q(1)%params(io_header_q(1)%size) = rtime
    io_header_q(1)%size = io_header_q(1)%size + 1; io_header_q(1)%params(io_header_q(1)%size) = visc
    io_header_q(1)%size = io_header_q(1)%size + 1; io_header_q(1)%params(io_header_q(1)%size) = froude
    io_header_q(1)%size = io_header_q(1)%size + 1; io_header_q(1)%params(io_header_q(1)%size) = rossby
    if (nse_eqns == DNS_EQNS_COMPRESSIBLE) then
        io_header_q(1)%size = io_header_q(1)%size + 1; io_header_q(1)%params(io_header_q(1)%size) = gamma0
        io_header_q(1)%size = io_header_q(1)%size + 1; io_header_q(1)%params(io_header_q(1)%size) = prandtl
        io_header_q(1)%size = io_header_q(1)%size + 1; io_header_q(1)%params(io_header_q(1)%size) = mach
    end if

    do is = 1, inb_scal
        io_header_s(is)%size = 0
        io_header_s(is)%size = io_header_s(is)%size + 1; io_header_s(is)%params(io_header_s(is)%size) = rtime
        io_header_s(is)%size = io_header_s(is)%size + 1; io_header_s(is)%params(io_header_s(is)%size) = visc
        io_header_s(is)%size = io_header_s(is)%size + 1; io_header_s(is)%params(io_header_s(is)%size) = schmidt(is)
    end do

    return
end subroutine TLab_Consistency_Check
