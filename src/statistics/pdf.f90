!########################################################################
!#
!# Calcualte {pdf(u_i), i=1,...,nv] in ny planes. (Histograms, not normalized.)
!# A last j-plane is added with the PDF in all the volume.
!#
!# ibc       In    BCs: 0 homogeneous interval
!#                      1 local interval
!#                      2 local interval, analysis and drop no point
!#                      3 local interval, analysis and drop left point
!#                      4 local interval, analysis and drop right point
!#                      5 local interval, analysis and drop both points
!#
!########################################################################
subroutine PDF1V_N(fname, nx, ny, nz, nv, nbins, ibc, umin, umax, u, igate, gate, z, pdf)
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: lfile
    use TLab_Time, only: itime, rtime
    use TLab_Pointers, only: pointers_dt
    use TLab_Arrays, only: wrk1d
    use TLab_WorkFlow, only: TLab_Write_ASCII
    use IO_Fields
    use PDFS
    use IO_PDFS
    implicit none

    character(len=*), intent(IN) :: fname
    integer(wi), intent(IN) :: nx, ny, nz, nv, nbins, ibc(nv)
    real(wp), intent(IN) :: umin(nv), umax(nv)              ! Random variables
    type(pointers_dt), intent(IN) :: u(nv)
    integer(1), intent(IN) :: gate(*), igate                ! discrete conditioning criteria
    real(wp), intent(IN) :: z(nz)                           ! heights of each plane
    real(wp), intent(inout) :: pdf(nbins + 2, nz + 1, nv)     ! last 2 bins contain the interval bounds

    ! -------------------------------------------------------------------
    integer(wi) iv, k, nplim, ibc_loc
    real(wp) plim, umin_loc, umax_loc
 
    ! ###################################################################
    call TLab_Write_ASCII(lfile, 'Calculating '//trim(adjustl(fname))//'...')

    plim = 1.0e-4_wp                    ! relative threshold in PDF analysis; adapt to sample size

    do iv = 1, nv

        do k = 1, nz                    ! calculation in planes
            if (igate == 0) then
                call PDF1V2D(ibc(iv), nx, ny, nz, k, umin(iv), umax(iv), u(iv)%field, nbins, pdf(1, k, iv), wrk1d)
            else
                ! call PDF1V2D1G(ibc(iv), nx, ny, nz, k, igate, gate, umin(iv), umax(iv), u(iv)%field, nbins, pdf(1, k, iv), wrk1d)
            end if

            if (ibc(iv) > 1) then
                ibc_loc = ibc(iv) - 2; umin_loc = umin(iv); umax_loc = umax(iv)
                call PDF_ANALYZE(ibc_loc, nbins, pdf(1, k, iv), umin_loc, umax_loc, plim, nplim)
                if (igate == 0) then
                    call PDF1V2D(0, nx, ny, nz, k, umin_loc, umax_loc, u(iv)%field, nbins, pdf(1, k, iv), wrk1d)
                else
                    ! call PDF1V2D1G(0, nx, ny, nz, k, igate, gate, umin_loc, umax_loc, u(iv)%field, nbins, pdf(1, k, iv), wrk1d)
                end if
            end if

        end do

        if (nz > 1) then                ! calculation in whole volume, saved as plane k=nz+1
            if (igate == 0) then
                call PDF1V2D(ibc(iv), nx, ny*nz, 1, 1, umin(iv), umax(iv), u(iv)%field, nbins, pdf(1, k, iv), wrk1d)
            else
                ! call PDF1V2D1G(ibc(iv), nx, ny*nz, 1, 1, igate, gate, umin(iv), umax(iv), u(iv)%field, nbins, pdf(1, k, iv), wrk1d)
            end if

            if (ibc(iv) > 1) then
                ibc_loc = ibc(iv) - 2; umin_loc = umin(iv); umax_loc = umax(iv)
                call PDF_ANALYZE(ibc_loc, nbins, pdf(1, k, iv), umin_loc, umax_loc, plim, nplim)
                if (igate == 0) then
                    call PDF1V2D(0, nx, ny*nz, 1, 1, umin_loc, umax_loc, u(iv)%field, nbins, pdf(1, k, iv), wrk1d)
                else
                    ! call PDF1V2D1G(0, nx, ny*nz, 1, 1, igate, gate, umin_loc, umax_loc, u(iv)%field, nbins, pdf(1, k, iv), wrk1d)
                end if
            end if

        end if

    end do

    ! ###################################################################
    call IO_Write_PDFs(fname, itime, rtime, z, pdf, u(1:nv)%tag)

    return

end subroutine PDF1V_N

!########################################################################
!########################################################################
subroutine PDF2V(fname, time, nx, ny, nz, nbins, u, v, z, pdf)
    use TLab_Constants, only: lfile, wp, wi
    use TLab_Arrays, only: wrk2d
    use TLab_WorkFlow, only: TLab_Write_ASCII
    use IO_Fields
    use PDFS
#ifdef USE_MPI
    use mpi_f08
#endif

    implicit none

    character(len=*), intent(IN) :: fname
    real(wp), intent(IN) :: time
    integer(wi), intent(IN) :: nx, ny, nz, nbins(2)
    real(wp), intent(IN) :: u(nx*ny*nz), v(nx*ny*nz)
    real(wp), intent(IN) :: z(nz)
    real(wp), intent(OUT) :: pdf(nbins(1)*nbins(2) + 2 + 2*nbins(1), nz + 1)

    ! -------------------------------------------------------------------
    integer(wi) k
    character*64 name

#ifdef USE_MPI
    integer ims_pro, ims_err
    call MPI_COMM_RANK(MPI_COMM_WORLD, ims_pro, ims_err)
#endif

    ! ###################################################################
    call TLab_Write_ASCII(lfile, 'Calculating '//trim(adjustl(fname))//'...')

    do k = 1, nz                ! calculation in planes
        call PDF2V2D(nx, ny, nz, k, u, v, nbins, pdf(1, k), wrk2d)
    end do

    if (nz > 1) then            ! calculation in whole volume, saved as plane nz+1
        call PDF2V2D(nx, ny*nz, 1, 1, u, v, nbins, pdf(1, k), wrk2d)
    end if

    ! ###################################################################
#ifdef USE_MPI
    if (ims_pro == 0) then
#endif

#define LOC_UNIT_ID 21
#define LOC_STATUS 'unknown'
        name = trim(adjustl(fname))
        call TLab_Write_ASCII(lfile, 'Writing field '//trim(adjustl(name))//'...')
        call IO_Open_File(name, LOC_STATUS, LOC_UNIT_ID)
        if (nz > 1) then
            write (LOC_UNIT_ID) SNGL(time), nz, nbins, SNGL(z(:)), SNGL(pdf(:, :))
        else
            write (LOC_UNIT_ID) SNGL(time), nz, nbins, SNGL(z(:)), SNGL(pdf(:, 1))
        end if
        close (LOC_UNIT_ID)
#ifdef USE_MPI
    end if
#endif

    return

end subroutine PDF2V
