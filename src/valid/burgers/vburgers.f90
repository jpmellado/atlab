program VBURGERS
    use TLab_Constants, only: wp, wi, big_wp
    use TLab_Constants, only: gfile, ifile
    use TLab_Time, only: itime
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Memory, only: imax, jmax, kmax, inb_txc
    use TLab_Memory, only: TLab_Initialize_Memory
    use TLab_Arrays
    use TLab_Pointers_3D, only: tmp1, tmp2
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_VARS
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use FDM, only: FDM_Initialize
    use NavierStokes !, only: NavierStokes_Initialize_Parameters, visc
    use Thermodynamics, only: Thermo_Initialize
    use Thermo_Anelastic, only: ribackground
    use Gravity, only: Gravity_Initialize
    use LargeScaleForcing, only: LargeScaleForcing_Initialize
    use TLab_Grid
    use IO_Fields
    use OPR_Partial
    use NSE_Burgers
    use TLab_Background, only: TLab_Initialize_Background

    implicit none

    real(wp), dimension(:, :, :), pointer :: a, b, c

    integer(wi) i, j, k
    real(wp) params(0)

! ###################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Trp_Initialize(ifile)
#endif

    call TLab_Grid_Read(gfile, x, y, z)
    call FDM_Initialize(ifile)

    call NavierStokes_Initialize_Parameters(ifile)
    ! visc = 0.0_wp                           ! inviscid
    call Thermo_Initialize(ifile)
    call Gravity_Initialize(ifile)
    call LargeScaleForcing_Initialize(ifile)

    inb_txc = 5
    call TLab_Initialize_Memory(__FILE__)

    call OPR_Partial_Initialize(ifile)

    call TLab_Initialize_Background(ifile)
    call NSE_Burgers_Initialize(ifile)

    a(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 3)
    b(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 4)
    c(1:imax, 1:jmax, 1:kmax) => txc(1:imax*jmax*kmax, 5)

    ! ###################################################################
    ! Define forcing term
    ! ###################################################################
    call IO_Read_Fields('field.inp', imax, jmax, kmax, itime, 1, 0, a, params)

    ! ###################################################################
    print *, new_line('a'), 'Derivative along x.'
    call OPR_Partial_X(OPR_P2_P1, imax, jmax, kmax, a, b, c)
    if (nse_eqns == DNS_EQNS_ANELASTIC) then
        do k = 1, kmax
            do j = 1, jmax
                do i = 1, imax
                    b(i, j, k) = b(i, j, k)*visc*ribackground(k) - a(i, j, k)*c(i, j, k)
                end do
            end do
        end do
    else
        b = b*visc - a*c
    end if
    ! call IO_Write_Fields('fieldXdirect.out', imax, jmax, kmax, itime, 1, b, io_header_s(1:1))

    c = 0.0_wp
    call NSE_AddBurgers_PerVolume_X(0, imax, jmax, kmax, a, c, tmp1, tmp2)
    ! call IO_Write_Fields('fieldXburgers.out', imax, jmax, kmax, itime, 1, c, io_header_s(1:1))

    call check(b, c, tmp1)!, 'fieldX.dif')

    ! ###################################################################
    if (y%size > 1) then

        print *, new_line('a'), 'Derivative along y.'
        call OPR_Partial_Y(OPR_P2_P1, imax, jmax, kmax, a, b, c)
        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            do k = 1, kmax
                do j = 1, jmax
                    do i = 1, imax
                        b(i, j, k) = b(i, j, k)*visc*ribackground(k) - a(i, j, k)*c(i, j, k)
                    end do
                end do
            end do
        else
            b = b*visc - a*c
        end if
        ! call IO_Write_Fields('fieldYdirect.out', imax, jmax, kmax, itime, 1, b, io_header_s(1:1))

        c = 0.0_wp
        call NSE_AddBurgers_PerVolume_Y(0, imax, jmax, kmax, a, c, tmp1, tmp2)
        ! call IO_Write_Fields('fieldYburgers.out', imax, jmax, kmax, itime, 1, c, io_header_s(1:1))

        call check(b, c, tmp1)!, 'fieldY.dif')

    end if

    ! ###################################################################
    ! Careful if you have subsidence activated
    print *, new_line('a'), 'Derivative along z.'
    call OPR_Partial_Z(OPR_P2_P1, imax, jmax, kmax, a, b, c)
    if (nse_eqns == DNS_EQNS_ANELASTIC) then
        do k = 1, kmax
            do j = 1, jmax
                do i = 1, imax
                    b(i, j, k) = b(i, j, k)*visc*ribackground(k) - a(i, j, k)*c(i, j, k)
                end do
            end do
        end do
    else
        b = b*visc - a*c
    end if
    ! call IO_Write_Fields('fieldZdirect.out', imax, jmax, kmax, itime, 1, b, io_header_s(1:1))

    c = 0.0_wp
    call NSE_AddBurgers_PerVolume_Z(0, imax, jmax, kmax, a, c, tmp1, rhou_in=a)
    ! call IO_Write_Fields('fieldZburgers.out', imax, jmax, kmax, itime, 1, c, io_header_s(1:1))

    call check(b, c, tmp1)!, 'fieldZ.dif')

    call TLab_Stop(0)

    ! ###################################################################
contains
    subroutine check(a1, a2, dif, name)
        real(wp), intent(in) :: a1(:, :, :), a2(:, :, :)
        real(wp), intent(inout) :: dif(:, :, :)
        character(len=*), optional :: name

        real(wp) dummy, error
#ifdef USE_MPI
        real(wp) sum_mpi
#endif

        dif = a2 - a1
        error = sum(dif**2)/real(size(a1), wp)
        dummy = sum(a1**2)/real(size(a1), wp)
#ifdef USE_MPI
        sum_mpi = error/real(ims_npro, wp)
        call MPI_ALLREDUCE(sum_mpi, error, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)
        sum_mpi = dummy/real(ims_npro, wp)
        call MPI_ALLREDUCE(sum_mpi, dummy, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ims_err)

        if (ims_pro == 0) then
#endif
            write (*, *) 'Solution L2-norm ...........:', sqrt(dummy)
            write (*, *) 'Relative error .............: ', sqrt(error)/sqrt(dummy)
#ifdef USE_MPI
        end if
#endif

        if (present(name)) then
            call IO_Write_Fields(name, imax, jmax, kmax, itime, 1, dif)
        end if

        return
    end subroutine check

end program VBURGERS
