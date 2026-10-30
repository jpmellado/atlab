#include "tlab_error.h"

module IO_Fields_NetCDF
    use TLab_Constants, only: wp, wi
    use TLab_WorkFlow, only: TLab_Stop, TLab_Log, TLab_Error
    implicit none
    private

    public :: IO_Write_NETCDF

contains
    subroutine IO_WRITE_NETCDF(fname, itime, rtime, x, y, z, vars)
        use TLab_Pointers_3D, only: pointers3d_dt
        use NETCDF

        character(len=*), intent(in) :: fname
        integer(wi), intent(in) :: itime
        real(wp), intent(in) :: rtime
        real(wp), intent(in) :: x(:), y(:), z(:)
        type(pointers3d_dt), intent(in) :: vars(:)

        ! -------------------------------------------------------------------
        integer fid, dtid, dxid, dyid, dzid, tid, xid, yid, zid, itid
        integer, allocatable :: vid(:)
        integer iv

        ! ###################################################################
#ifdef USE_MPI
        call TLab_Error(__FILE__//'NETCDF IO not implemented in parallel mode.')
        call TLab_Stop(DNS_ERROR_UNDEVELOP)
#endif

        call TLab_Log('Writing '//trim(adjustl(fname))//'...')

        call NC_CHECK(NF90_CREATE(trim(adjustl(fname))//'.nc', NF90_NETCDF4, fid))

        call NC_CHECK(NF90_DEF_DIM(fid, "t", NF90_UNLIMITED, dtid))
        call NC_CHECK(NF90_DEF_DIM(fid, "x", size(x), dxid))
        call NC_CHECK(NF90_DEF_DIM(fid, "y", size(y), dyid))
        call NC_CHECK(NF90_DEF_DIM(fid, "z", size(z), dzid))

        call NC_CHECK(Nf90_DEF_VAR(fid, "t", NF90_FLOAT, (/dtid/), tid))
        call NC_CHECK(Nf90_DEF_VAR(fid, "x", NF90_FLOAT, (/dxid/), xid))
        call NC_CHECK(Nf90_DEF_VAR(fid, "y", NF90_FLOAT, (/dyid/), yid))
        call NC_CHECK(Nf90_DEF_VAR(fid, "z", NF90_FLOAT, (/dzid/), zid))
        call NC_CHECK(Nf90_DEF_VAR(fid, "it", NF90_INT, (/dtid/), itid))

        if (allocated(vid)) deallocate (vid)
        allocate (vid(size(vars)))
        do iv = 1, size(vars)
            call NC_CHECK(Nf90_DEF_VAR(fid, trim(adjustl(vars(iv)%tag)), NF90_FLOAT, (/dxid, dyid, dzid, dtid/), vid(iv)))
        end do

        call NC_CHECK(NF90_ENDDEF(fid))

        call NC_CHECK(NF90_PUT_VAR(fid, tid, SNGL(rtime)))
        call NC_CHECK(NF90_PUT_VAR(fid, itid, itime))
        call NC_CHECK(NF90_PUT_VAR(fid, xid, SNGL(x)))
        call NC_CHECK(NF90_PUT_VAR(fid, yid, SNGL(y)))
        call NC_CHECK(NF90_PUT_VAR(fid, zid, SNGL(z)))
        do iv = 1, size(vars)
            call NC_CHECK(NF90_PUT_VAR(fid, vid(iv), SNGL(vars(iv)%field(1:size(x), 1:size(y), 1:size(z)))))
        end do
        call NC_CHECK(NF90_CLOSE(fid))

        return
    contains
! ###################################################################
        subroutine NC_CHECK(status)
            integer, intent(in) :: status

            if (status /= nf90_noerr) then
                call TLab_Error('NETCDF error signal '//trim(adjustl(NF90_STRERROR(status))))
                call TLab_Stop(DNS_ERROR_UNDEVELOP)
            end if

            return
        end subroutine NC_CHECK

    end subroutine IO_WRITE_NETCDF

end module
