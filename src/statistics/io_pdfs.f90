!########################################################################
!#
!# Write N PDFs in netCDF or ASCII format
!#
!########################################################################
module IO_PDFS
contains
    subroutine IO_Write_PDFs(fname, itime, rtime, z, pdf, varnames)
        use TLab_Constants, only: wp, sp, wi
        use TLab_Constants, only: lfile
        use TLab_WorkFlow, only: TLab_Write_ASCII
#ifdef USE_NETCDF
        use NETCDF
#else
        use IO_Fields
#endif
#ifdef USE_MPI
        use mpi_f08, only: MPI_COMM_WORLD, MPI_COMM_RANK
#endif
        implicit none

        character(len=*), intent(in) :: fname
        integer(wi), intent(in) :: itime
        real(wp), intent(in) :: rtime
        real(wp), intent(in) :: z(:)
        real(wp), intent(in) :: pdf(:, :, :)
        character(len=*), intent(in) :: varnames(size(pdf, 3))

        ! -------------------------------------------------------------------
        integer nbins, nv, iv, nz

#ifdef USE_NETCDF
        integer fid, dtid, dbid, drid, dzid, tid, zid, itid
        integer, allocatable :: vid(:), vid_range(:)
        integer, allocatable :: vid_global(:), vid_global_range(:)
#else
        character(len=64) name
#endif

#ifdef USE_MPI
        integer locRank, ims_err
        call MPI_COMM_RANK(MPI_COMM_WORLD, locRank, ims_err)
#endif

        ! ###################################################################
        call TLab_Write_ASCII(lfile, 'Writing '//trim(adjustl(fname))//'...')

#ifdef USE_MPI
        if (locRank == 0) then
#endif

            nbins = size(pdf, 1) - 2    ! First index is random variable; last 2 values are min and max
            nz = size(pdf, 2) - 1       ! Second index is height; last one is global PDF
            nv = size(pdf, 3)           ! Third index is variable

            if (allocated(vid)) deallocate (vid)
            allocate (vid(1:nv))
            if (allocated(vid_range)) deallocate (vid_range)
            allocate (vid_range(1:nv))
            if (allocated(vid_global)) deallocate (vid_global)
            allocate (vid_global(1:nv))
            if (allocated(vid_global_range)) deallocate (vid_global_range)
            allocate (vid_global_range(1:nv))

            ! -----------------------------------------------------------------------
            ! Using NetCDF format
#ifdef USE_NETCDF
            call NC_CHECK(NF90_CREATE(trim(adjustl(fname))//'.nc', NF90_NETCDF4, fid))

            call NC_CHECK(NF90_DEF_DIM(fid, "t", NF90_UNLIMITED, dtid))
            call NC_CHECK(NF90_DEF_DIM(fid, "bins", nbins, dbid))
            call NC_CHECK(NF90_DEF_DIM(fid, "z", nz, dzid))
            call NC_CHECK(NF90_DEF_DIM(fid, "range", 2, drid))

            call NC_CHECK(Nf90_DEF_VAR(fid, "t", NF90_FLOAT, (/dtid/), tid))
            call NC_CHECK(Nf90_DEF_VAR(fid, "z", NF90_FLOAT, (/dzid/), zid))
            call NC_CHECK(Nf90_DEF_VAR(fid, "it", NF90_INT, (/dtid/), itid))
            do iv = 1, nv
                call NC_CHECK(Nf90_DEF_VAR(fid, trim(adjustl(varnames(iv))), NF90_FLOAT, (/dbid, dzid, dtid/), vid(iv)))
                call NC_CHECK(Nf90_DEF_VAR(fid, trim(adjustl(varnames(iv)))//'_range', NF90_FLOAT, (/drid, dzid, dtid/), vid_range(iv)))
                call NC_CHECK(Nf90_DEF_VAR(fid, trim(adjustl(varnames(iv)))//'_global', NF90_FLOAT, (/dbid, dtid/), vid_global(iv)))
                call NC_CHECK(Nf90_DEF_VAR(fid, trim(adjustl(varnames(iv)))//'_global_range', NF90_FLOAT, (/drid, dtid/), vid_global_range(iv)))
            end do

            call NC_CHECK(NF90_ENDDEF(fid))

            call NC_CHECK(NF90_PUT_VAR(fid, tid, real(rtime, sp)))
            call NC_CHECK(NF90_PUT_VAR(fid, itid, itime))
            call NC_CHECK(NF90_PUT_VAR(fid, zid, real(z, sp)))
            do iv = 1, nv
                call NC_CHECK(NF90_PUT_VAR(fid, vid(iv), real(pdf(1:nbins, 1:nz, iv), sp)))
                call NC_CHECK(NF90_PUT_VAR(fid, vid_range(iv), real(pdf(nbins + 1:nbins + 2, 1:nz, iv), sp)))
                call NC_CHECK(NF90_PUT_VAR(fid, vid_global(iv), real(pdf(1:nbins, nz + 1, iv), sp)))
                call NC_CHECK(NF90_PUT_VAR(fid, vid_global_range(iv), real(pdf(nbins + 1:nbins + 2, nz + 1, iv), sp)))
            end do

            call NC_CHECK(NF90_CLOSE(fid))

            ! -----------------------------------------------------------------------
            ! Using ASCII format; to be checked
#else

#define LOC_UNIT_ID 21
#define LOC_STATUS 'unknown'
            do iv = 1, nv
                name = trim(adjustl(fname))
                if (varnames(iv) /= '') name = trim(adjustl(fname))//'.'//trim(adjustl(varnames(iv)))
                call TLab_Write_ASCII(lfile, 'Writing field '//trim(adjustl(name))//'...')
                call IO_Open_File(name, LOC_STATUS, LOC_UNIT_ID)
                if (nz > 1) then
                    write (LOC_UNIT_ID) SNGL(rtime), nz, nbins, SNGL(z(:)), SNGL(pdf(:, :, iv))
                else
                    write (LOC_UNIT_ID) SNGL(rtime), nz, nbins, SNGL(z(:)), SNGL(pdf(:, 1, iv))
                end if
                close (LOC_UNIT_ID)
            end do

#endif

#ifdef USE_MPI
        end if
#endif

        return
    contains
! ###################################################################
#include "tlab_error.h"

#ifdef USE_NETCDF
        subroutine NC_CHECK(status)
            use TLab_Constants, only: efile
            use TLab_WorkFlow, only: TLab_Stop
            integer, intent(in) :: status

            if (status /= nf90_noerr) then
                call TLab_Write_ASCII(efile, __FILE__//'. NETCDF error signal '//trim(adjustl(NF90_STRERROR(status))))
                call TLab_Stop(DNS_ERROR_UNDEVELOP)
            end if

            return
        end subroutine NC_CHECK
#endif

    end subroutine IO_Write_PDFs
end module
