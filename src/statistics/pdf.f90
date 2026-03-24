module StatsPDFs
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: lfile
    use TLab_WorkFlow, only: TLab_Write_ASCII
    use TLab_Time, only: itime, rtime
    use PDFS
    use IO_PDFS
    implicit none
    private

    public :: PDF1V_N
    public :: PDF2V

contains
!########################################################################
!#
!# Calculate {pdf(u_i), i=1,...,nv] in nz planes. (Histograms, not normalized.)
!# A last k-plane is added with the PDF in all the volume.
!#
!# ibc       In    BCs: 0 homogeneous interval
!#                      1 local interval
!#                      2 local interval, analysis and drop no point
!#                      3 local interval, analysis and drop left point
!#                      4 local interval, analysis and drop right point
!#                      5 local interval, analysis and drop both points
!#
!########################################################################
    subroutine PDF1V_N(fname, nx, ny, nz, nv, nbins, ibc, umin, umax, u, z, pdf, maskGate)
        use TLab_Pointers, only: pointers_dt
        use TLab_Arrays, only: wrk1d

        character(len=*), intent(in) :: fname
        integer(wi), intent(in) :: nx, ny, nz, nv, nbins, ibc(nv)
        real(wp), intent(in) :: umin(nv), umax(nv)                  ! Random variables
        type(pointers_dt), intent(in) :: u(nv)
        real(wp), intent(in) :: z(nz)                               ! heights of each plane
        real(wp), intent(inout) :: pdf(nbins + 2, nz + 1, nv)       ! last 2 bins contain the interval bounds
        logical, intent(in), optional :: maskGate(:)                ! discrete conditioning criteria

        ! -------------------------------------------------------------------
        integer iv, k
        integer nplim, ibc_loc
        real(wp) plim, umin_loc, umax_loc

        ! ###################################################################
        call TLab_Write_ASCII(lfile, 'Calculating '//trim(adjustl(fname))//'...')

        plim = 1.0e-4_wp                    ! relative threshold in PDF analysis; adapt to sample size

        if (present(maskGate)) then
            ! same but using conditional routines like this one:
            ! call PDF1V2D1G(0, nx, ny, nz, k, igate, gate, umin_loc, umax_loc, u(iv)%field, nbins, pdf(1, k, iv), wrk1d)

        else
            do iv = 1, nv
                do k = 1, nz
                    call PDF1V2D(ibc(iv), nx, ny, nz, k, umin(iv), umax(iv), u(iv)%field, nbins, pdf(1, k, iv), wrk1d)
                    if (ibc(iv) > 1) then
                        ibc_loc = ibc(iv) - 2; umin_loc = umin(iv); umax_loc = umax(iv)
                        call PDF_ANALYZE(ibc_loc, nbins, pdf(1, k, iv), umin_loc, umax_loc, plim, nplim)
                        call PDF1V2D(0, nx, ny, nz, k, umin_loc, umax_loc, u(iv)%field, nbins, pdf(1, k, iv), wrk1d)
                    end if

                end do

                ! PDF over the whole domain
                call PDF1V2D(ibc(iv), nx, ny*nz, 1, 1, umin(iv), umax(iv), u(iv)%field, nbins, pdf(1, k, iv), wrk1d)
                if (ibc(iv) > 1) then
                    ibc_loc = ibc(iv) - 2; umin_loc = umin(iv); umax_loc = umax(iv)
                    call PDF_ANALYZE(ibc_loc, nbins, pdf(1, k, iv), umin_loc, umax_loc, plim, nplim)
                    call PDF1V2D(0, nx, ny*nz, 1, 1, umin_loc, umax_loc, u(iv)%field, nbins, pdf(1, k, iv), wrk1d)
                end if

            end do

        end if

        ! ###################################################################
        call IO_Write_PDFs_1D(fname, itime, rtime, z, pdf, u(1:nv)%tag)

        return
    end subroutine PDF1V_N

    !########################################################################
    !########################################################################
    subroutine PDF2V(fname, nx, ny, nz, nbins, u, v, z, pdfMem)
        use TLab_Arrays, only: wrk2d
        character(len=*), intent(in) :: fname
        integer(wi), intent(in) :: nx, ny, nz, nbins(2)
        real(wp), intent(in) :: u(nx*ny*nz), v(nx*ny*nz)
        real(wp), intent(in) :: z(nz)
        real(wp), intent(out) :: pdfMem((nbins(1)*nbins(2) + 2 + 2*nbins(1))*(nz + 1))
        target pdfMem

        ! -------------------------------------------------------------------
        integer(wi) k, offset

        real(wp), pointer :: pdf(:, :, :) => null()
        real(wp), pointer :: ranges(:, :, :) => null()

        ! ###################################################################
        call TLab_Write_ASCII(lfile, 'Calculating '//trim(adjustl(fname))//'...')

        offset = nbins(1)*nbins(2)*(nz + 1)

        pdf(1:nbins(1), 1:nbins(2), 1:nz + 1) => pdfMem(1:offset)
        ranges(1:2, 0:nbins(1), 1:nz + 1) => pdfMem(1 + offset:)

        do k = 1, nz                ! calculation in planes
            call PDF2V2D(nx, ny, nz, k, u, v, nbins, pdf(:, :, k), ranges(:, :, k), wrk2d)
        end do

        if (nz > 1) then            ! calculation in whole volume, saved as plane nz+1
            call PDF2V2D(nx, ny*nz, 1, 1, u, v, nbins, pdf(:, :, k), ranges(:, :, k), wrk2d)
        end if

        ! ###################################################################
        call IO_Write_PDFs_2D(fname, itime, rtime, z, &
                              pdf(:, :, :), &
                              ranges(:, :, :))

        return
    end subroutine PDF2V

end module
