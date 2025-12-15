! Calculate the non-linear operator N(u)(s) = dyn_visc* d^2/dx^2 s - rho u d/dx s

module NSE_Burgers
    use TLab_Constants, only: wp, wi, efile, lfile, BCS_NONE, MAX_VARS
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Arrays, only: wrk2d, wrk3d
    use TLab_Transpose
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_j
    use TLabMPI_Transpose
#endif
    use TLab_Grid, only: x, y, z
    use FDM, only: fdm_dt, g
    use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC
    use Thermo_Anelastic, only: ribackground, rbackground
    use NavierStokes, only: visc, schmidt
    use FDM_Derivative, only: FDM_Der1_Solve, FDM_Der2_Solve
#ifdef USE_MPI
    use OPR_Partial
#endif
    use LargeScaleForcing, only: subsidenceProps, TYPE_SUB_CONSTANT, wbackground
    implicit none
    private

    public :: NSE_Burgers_Initialize
    public :: NSE_AddBurgers_PerVolume_X
    public :: NSE_AddBurgers_PerVolume_Y
    public :: NSE_AddBurgers_PerVolume_Z

    ! -----------------------------------------------------------------------
    procedure(nse_burgers_ice) :: NSE_AddBurgers_PerVolume_dt
    abstract interface
        subroutine nse_burgers_ice(is, nx, ny, nz, s, rhs, result, tmp1, rhou_in)
            use TLab_Constants, only: wi, wp
            integer, intent(in) :: is                           ! scalar index; if 0, then velocity
            integer(wi), intent(in) :: nx, ny, nz
            real(wp), intent(in) :: s(nx*ny*nz)
            real(wp), intent(out) :: result(nx*ny*nz)
            real(wp), intent(inout) :: rhs(nx*ny*nz)
            real(wp), intent(out) :: tmp1(nx*ny*nz)             ! transposed field s times density
            real(wp), intent(in), optional :: rhou_in(nx*ny*nz) ! transposed field u times density
        end subroutine
    end interface
    procedure(NSE_AddBurgers_PerVolume_dt), pointer :: NSE_AddBurgers_PerVolume_X, NSE_AddBurgers_PerVolume_Y

    type :: fdm_diffusion_dt
        sequence
        real(wp), allocatable :: lu(:, :, :)
    end type fdm_diffusion_dt
    type(fdm_diffusion_dt) :: fdmDiffusion(3)

    real(wp), allocatable :: rho_yz(:), rho_xz(:)   ! rho in anelastic formulation, x and y directions
    real(wp) :: diffusivity(0:MAX_VARS)

    real(wp), allocatable :: rho_wbackground(:)     ! subsidence velocity (times density)

contains
    !########################################################################
    !########################################################################
    subroutine NSE_Burgers_Initialize(inifile)
        use TLab_Memory, only: imax, jmax, kmax
        use TLab_Memory, only: TLab_Allocate_Real
        use TLab_Memory, only: inb_scal
#ifdef USE_MPI
        use TLabMPI_VARS, only: ims_pro_i, ims_npro_i, ims_pro_j, ims_npro_j
#endif

        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        character(len=32) bakfile

        integer(wi) ig, is, ip, j, ndl, idl
        integer(wi) nlines, offset

        ! ###################################################################
        ! Read input data
        bakfile = trim(adjustl(inifile))//'.bak'

        ! ###################################################################
        ! Initialize LU factorization of the second-order derivative times the diffusivity
        do ig = 1, 3
            if (g(ig)%size == 1) cycle

            ndl = g(ig)%der2%nb_diag(1)
            if (g(ig)%periodic) then
                allocate (fdmDiffusion(ig)%lu(g(ig)%size, ndl + 2, 0:inb_scal))
            else
                allocate (fdmDiffusion(ig)%lu(g(ig)%size, ndl, 0:inb_scal))
            end if

            idl = ndl/2 + 1
            do is = 0, inb_scal ! case 0 for the reynolds number
                if (is == 0) then
                    diffusivity(is) = visc
                else
                    diffusivity(is) = visc/schmidt(is)
                end if

                fdmDiffusion(ig)%lu(:, :, is) = g(ig)%der2%lu(:, :)                 ! Check routines Thomas3C_LU and Thomas3C_Solve
                fdmDiffusion(ig)%lu(:, idl, is) = g(ig)%der2%lu(:, idl)*diffusivity(is)

            end do
        end do

        ! ###################################################################
        ! Initialize anelastic density correction
        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call TLab_Write_ASCII(lfile, 'Initialize anelastic density correction in burgers operator.')

            ! -----------------------------------------------------------------------
            ! Density correction term in the burgers operator along X
#ifdef USE_MPI
            if (ims_npro_i > 1 .and. der_mode_i == TYPE_TRANSPOSE) then
                nlines = tmpi_plan_dx%nlines
                offset = nlines*ims_pro_i
            else
#endif
                nlines = jmax*kmax
                offset = 0
#ifdef USE_MPI
            end if
#endif
            allocate (rho_yz(nlines))
            do j = 1, nlines
                ip = (offset + j - 1)/jmax + 1
                rho_yz(j) = rbackground(ip)
            end do

            ! -----------------------------------------------------------------------
            ! Density correction term in the burgers operator along Y
#ifdef USE_MPI
            if (ims_npro_j > 1 .and. der_mode_j == TYPE_TRANSPOSE) then
                nlines = tmpi_plan_dy%nlines
                offset = nlines*ims_pro_j
            else
#endif
                nlines = imax*kmax
                offset = 0
#ifdef USE_MPI
            end if
#endif
            allocate (rho_xz(nlines))
            do j = 1, nlines
                ip = mod(offset + j - 1, z%size) + 1
                rho_xz(j) = rbackground(ip)
            end do

        end if

        ! ###################################################################
        ! Initialize subsidence velocity (times density) to handle both Boussinesq and anelastic
        allocate (rho_wbackground(z%size))
        if (nse_eqns == DNS_EQNS_ANELASTIC) then    ! evolution equations per unit volume
            rho_wbackground(:) = wbackground(:)*rbackground(:)
        else
            rho_wbackground(:) = wbackground(:)
        end if

        ! ###################################################################
        ! Setting procedure pointers
#ifdef USE_MPI
        if (ims_npro_i > 1) then
            select case (der_mode_i)
            case (TYPE_TRANSPOSE)
                NSE_AddBurgers_PerVolume_X => NSE_AddBurgers_PerVolume_X_MPITranspose
            case (TYPE_SPLIT)
                NSE_AddBurgers_PerVolume_X => NSE_AddBurgers_PerVolume_X_MPISplit
            end select
        else
#endif
            NSE_AddBurgers_PerVolume_X => NSE_AddBurgers_PerVolume_X_Serial
#ifdef USE_MPI
        end if
#endif

#ifdef USE_MPI
        if (ims_npro_j > 1) then
            select case (der_mode_j)
            case (TYPE_TRANSPOSE)
                NSE_AddBurgers_PerVolume_Y => NSE_AddBurgers_PerVolume_Y_MPITranspose
            case (TYPE_SPLIT)
                NSE_AddBurgers_PerVolume_Y => NSE_AddBurgers_PerVolume_Y_MPISplit
            end select
        else
#endif
            NSE_AddBurgers_PerVolume_Y => NSE_AddBurgers_PerVolume_Y_Serial
#ifdef USE_MPI
        end if
#endif
        return
    end subroutine NSE_Burgers_Initialize

    !########################################################################
    !########################################################################
    subroutine NSE_AddBurgers_PerVolume_X_Serial(is, nx, ny, nz, s, rhs, result, tmp1, rhou_in)
        integer, intent(in) :: is                           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout) :: rhs(nx*ny*nz)
        real(wp), intent(out) :: tmp1(nx*ny*nz)             ! transposed field s times density
        real(wp), intent(in), optional :: rhou_in(nx*ny*nz) ! transposed field u times density

        ! -------------------------------------------------------------------
        integer(wi) nlines

        ! ###################################################################
        if (x%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            return
        end if

        ! Transposition: make x-direction the last one
#ifdef USE_ESSL
        call DGETMO(s, nx, nx, ny*nz, tmp1, ny*nz)
#else
        call TLab_Transpose_Real(s, nx, ny*nz, nx, tmp1, ny*nz)
#endif

        nlines = ny*nz

        call FDM_Der1_Solve(nlines, g(1)%der1, g(1)%der1%lu, tmp1, result, wrk2d)
        call FDM_Der2_Solve(nlines, g(1)%der2, fdmDiffusion(1)%lu(:, :, is), tmp1, wrk3d, result, wrk2d)

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            wrk3d(1:nx*ny*nz) = wrk3d(1:nx*ny*nz) - rhou_in(:)*result(:)
        else
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call NSE_Burgers_1D(nlines, nx, &
                                    der1=result, &
                                    der2=wrk3d, &
                                    rhou=tmp1, &
                                    rho_xy=rho_yz(:))
            else
                wrk3d(1:nx*ny*nz) = wrk3d(1:nx*ny*nz) - tmp1(:)*result(:)
            end if
        end if

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(wrk3d, ny*nz, ny*nz, nx, result, nx)
        rhs = rhs + result
#else
        call TLab_AddTranspose(wrk3d, ny*nz, nx, ny*nz, rhs, nx)
#endif

        return
    end subroutine NSE_AddBurgers_PerVolume_X_Serial

    !########################################################################
    !########################################################################
#ifdef USE_MPI
    subroutine NSE_AddBurgers_PerVolume_X_MPITranspose(is, nx, ny, nz, s, rhs, result, tmp1, rhou_in)
        integer, intent(in) :: is                           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout) :: rhs(nx*ny*nz)
        real(wp), intent(out) :: tmp1(nx*ny*nz)             ! transposed field s times density
        real(wp), intent(in), optional :: rhou_in(nx*ny*nz) ! transposed field u times density

        ! -------------------------------------------------------------------
        integer(wi) nlines

        ! ###################################################################
        if (x%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            return
        end if

        nlines = tmpi_plan_dx%nlines

        ! Transposition: make x-direction the last one
        call TLabMPI_Trp_ExecI_Forward(s, result, tmpi_plan_dx)
#ifdef USE_ESSL
        call DGETMO(result, g(1)%size, g(1)%size, nlines, tmp1, nlines)
#else
        call TLab_Transpose_Real(result, g(1)%size, nlines, g(1)%size, tmp1, nlines)
#endif

        call FDM_Der1_Solve(nlines, g(1)%der1, g(1)%der1%lu, tmp1, wrk3d, wrk2d)
        call FDM_Der2_Solve(nlines, g(1)%der2, fdmDiffusion(1)%lu(:, :, is), tmp1, result, wrk3d, wrk2d)

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            result(:) = result(:) - rhou_in(:)*wrk3d(1:nx*ny*nz)
        else
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call NSE_Burgers_1D(nlines, nx*ims_npro_i, &
                                    der1=wrk3d, &
                                    der2=result, &
                                    rhou=tmp1, &
                                    rho_xy=rho_yz(:))
            else
                result(:) = result(:) - tmp1(:)*wrk3d(1:nx*ny*nz)
            end if
        end if

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(result, nlines, nlines, g(1)%size, wrk3d, g(1)%size)
#else
        call TLab_Transpose_Real(result, nlines, g(1)%size, nlines, wrk3d, g(1)%size)
#endif
        call TLabMPI_Trp_ExecI_Backward(wrk3d, result, tmpi_plan_dx)
        rhs = rhs + result

        return
    end subroutine NSE_AddBurgers_PerVolume_X_MPITranspose

    !########################################################################
    !########################################################################
    subroutine NSE_AddBurgers_PerVolume_X_MPISplit(is, nx, ny, nz, s, rhs, result, tmp1, rhou_in)
        use TLabMPI_PROCS, only: TLabMPI_Halos_X
        use FDM_Derivative_MPISplit
        integer, intent(in) :: is                           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout) :: rhs(nx*ny*nz)
        real(wp), intent(out) :: tmp1(nx*ny*nz)             ! transposed field s times density
        real(wp), intent(in), optional :: rhou_in(nx*ny*nz) ! transposed field u times density

        ! -------------------------------------------------------------------
        integer(wi) nlines
        integer np, np1, np2

        ! ###################################################################
        if (x%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            return
        end if

        ! Transposition: make x-direction the last one
#ifdef USE_ESSL
        call DGETMO(s, nx, nx, ny*nz, tmp1, ny*nz)
#else
        call TLab_Transpose_Real(s, nx, ny*nz, nx, tmp1, ny*nz)
#endif

        nlines = ny*nz

        np1 = size(der1_split_x%rhs)/2
        np2 = size(der2_split_x%rhs)/2
        np = max(np1, np2)
        call TLabMPI_Halos_X(tmp1, nlines, np, pyz_halo_m(:, 1), pyz_halo_p(:, 1))

        call FDM_MPISplit_Solve(nlines, nx, der1_split_x, tmp1, &
                                pyz_halo_m(:, np - np1 + 1:np), pyz_halo_p, result, wrk2d)
        call FDM_MPISplit_Solve(nlines, nx, der2_split_x, tmp1, &
                                pyz_halo_m(:, np - np2 + 1:np), pyz_halo_p, wrk3d, wrk2d)

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            wrk3d(1:nx*ny*nz) = wrk3d(1:nx*ny*nz)*diffusivity(is) - rhou_in(:)*result(:)
        else
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call NSE_Burgers_1D_Split(nlines, nx, &
                                          der1=result, &
                                          der2=wrk3d, &
                                          rhou=tmp1, &
                                          rho_xy=rho_yz(:), &
                                          diff=diffusivity(is))
            else
                wrk3d(1:nx*ny*nz) = wrk3d(1:nx*ny*nz)*diffusivity(is) - tmp1(:)*result(:)
            end if
        end if

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(wrk3d, ny*nz, ny*nz, nx, result, nx)
        rhs = rhs + result
#else
        call TLab_AddTranspose(wrk3d, ny*nz, nx, ny*nz, rhs, nx)
#endif

        return
    end subroutine NSE_AddBurgers_PerVolume_X_MPISplit

#endif

    !########################################################################
    !########################################################################
    subroutine NSE_AddBurgers_PerVolume_Y_Serial(is, nx, ny, nz, s, rhs, result, tmp1, rhou_in)
        integer, intent(in) :: is                           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout) :: rhs(nx*ny*nz)
        real(wp), intent(out) :: tmp1(nx*ny*nz)             ! transposed field s times density
        real(wp), intent(in), optional :: rhou_in(nx*ny*nz) ! transposed field u times density

        ! -------------------------------------------------------------------
        integer(wi) nlines

        ! ###################################################################
        if (y%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            return
        end if

        ! Transposition: make y-direction the last one
#ifdef USE_ESSL
        call DGETMO(s, nx*ny, nx*ny, nz, tmp1, nz)
#else
        call TLab_Transpose_Real(s, nx*ny, nz, nx*ny, tmp1, nz, locBlock=trans_y)
#endif

        nlines = nx*nz

        call FDM_Der1_Solve(nlines, g(2)%der1, g(2)%der1%lu, tmp1, result, wrk2d)
        call FDM_Der2_Solve(nlines, g(2)%der2, fdmDiffusion(2)%lu(:, :, is), tmp1, wrk3d, result, wrk2d)

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            wrk3d(1:nx*ny*nz) = wrk3d(1:nx*ny*nz) - rhou_in(:)*result(:)
        else
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call NSE_Burgers_1D(nlines, ny, &
                                    der1=result, &
                                    der2=wrk3d, &
                                    rhou=tmp1, &
                                    rho_xy=rho_xz(:))
            else
                wrk3d(1:nx*ny*nz) = wrk3d(1:nx*ny*nz) - tmp1(:)*result(:)
            end if

        end if

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(wrk3d, nz, nz, nx*ny, result, nx*ny)
        rhs = rhs + result
#else
        call TLab_AddTranspose(wrk3d, nz, nx*ny, nz, rhs, nx*ny, locBlock=trans_y)
#endif

        return
    end subroutine NSE_AddBurgers_PerVolume_Y_Serial

    !########################################################################
    !########################################################################
#ifdef USE_MPI
    subroutine NSE_AddBurgers_PerVolume_Y_MPITranspose(is, nx, ny, nz, s, rhs, result, tmp1, rhou_in)
        integer, intent(in) :: is                           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout) :: rhs(nx*ny*nz)
        real(wp), intent(out) :: tmp1(nx*ny*nz)             ! transposed field s times density
        real(wp), intent(in), optional :: rhou_in(nx*ny*nz) ! transposed field u times density

        ! -------------------------------------------------------------------
        integer(wi) nlines

        ! ###################################################################
        if (y%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            return
        end if

        ! Transposition: make y-direction the last one
#ifdef USE_ESSL
        call DGETMO(s, nx*ny, nx*ny, nz, wrk3d, nz)
#else
        call TLab_Transpose_Real(s, nx*ny, nz, nx*ny, wrk3d, nz)
#endif
        call TLabMPI_Trp_ExecJ_Forward(wrk3d, tmp1, tmpi_plan_dy)
        nlines = tmpi_plan_dy%nlines

        call FDM_Der1_Solve(nlines, g(2)%der1, g(2)%der1%lu, tmp1, wrk3d, wrk2d)
        call FDM_Der2_Solve(nlines, g(2)%der2, fdmDiffusion(2)%lu(:, :, is), tmp1, result, wrk3d, wrk2d)

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            result(:) = result(:) - rhou_in(:)*wrk3d(1:nx*ny*nz)
        else
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call NSE_Burgers_1D(nlines, ny*ims_npro_j, &
                                    der1=wrk3d, &
                                    der2=result, &
                                    rhou=tmp1, &
                                    rho_xy=rho_xz(:))
            else
                result(:) = result(:) - tmp1(:)*wrk3d(1:nx*ny*nz)
            end if

        end if

        ! Put arrays back in the order in which they came in
        call TLabMPI_Trp_ExecJ_Backward(result, wrk3d, tmpi_plan_dy)
#ifdef USE_ESSL
        call DGETMO(wrk3d, nz, nz, nx*ny, result, nx*ny)
        rhs = rhs + result
#else
        call TLab_AddTranspose(wrk3d, nz, nx*ny, nz, rhs, nx*ny)
#endif

        return
    end subroutine NSE_AddBurgers_PerVolume_Y_MPITranspose

    !########################################################################
    !########################################################################
    subroutine NSE_AddBurgers_PerVolume_Y_MPISplit(is, nx, ny, nz, s, rhs, result, tmp1, rhou_in)
        use TLabMPI_PROCS, only: TLabMPI_Halos_Y
        use FDM_Derivative_MPISplit
        integer, intent(in) :: is                           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout) :: rhs(nx*ny*nz)
        real(wp), intent(out) :: tmp1(nx*ny*nz)             ! transposed field s times density
        real(wp), intent(in), optional :: rhou_in(nx*ny*nz) ! transposed field u times density

        ! -------------------------------------------------------------------
        integer(wi) nlines
        integer np, np1, np2

        ! ###################################################################
        if (y%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            return
        end if

        ! Transposition: make y-direction the last one
#ifdef USE_ESSL
        call DGETMO(s, nx*ny, nx*ny, nz, tmp1, nz)
#else
        call TLab_Transpose_Real(s, nx*ny, nz, nx*ny, tmp1, nz, locBlock=trans_y)
#endif

        nlines = nx*nz

        np1 = size(der1_split_y%rhs)/2
        np2 = size(der2_split_y%rhs)/2
        np = max(np1, np2)
        call TLabMPI_Halos_Y(tmp1, nlines, np, pxz_halo_m(:, 1), pxz_halo_p(:, 1))

        call FDM_MPISplit_Solve(nlines, ny, der1_split_y, tmp1, &
                                pxz_halo_m(:, np - np1 + 1:np), pxz_halo_p, result, wrk2d)
        call FDM_MPISplit_Solve(nlines, ny, der2_split_y, tmp1, &
                                pxz_halo_m(:, np - np2 + 1:np), pxz_halo_p, wrk3d, wrk2d)

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            wrk3d(1:nx*ny*nz) = wrk3d(1:nx*ny*nz)*diffusivity(is) - rhou_in(:)*result(:)
        else
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call NSE_Burgers_1D_Split(nlines, ny, &
                                          der1=result, &
                                          der2=wrk3d, &
                                          rhou=tmp1, &
                                          rho_xy=rho_xz(:), &
                                          diff=diffusivity(is))
            else
                wrk3d(1:nx*ny*nz) = wrk3d(1:nx*ny*nz)*diffusivity(is) - tmp1(:)*result(:)
            end if

        end if

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(wrk3d, nz, nz, nx*ny, result, nx*ny)
        rhs = rhs + result
#else
        call TLab_AddTranspose(wrk3d, nz, nx*ny, nz, rhs, nx*ny, locBlock=trans_y)
#endif

        return
    end subroutine NSE_AddBurgers_PerVolume_Y_MPISplit

#endif

    !########################################################################
    !########################################################################
    subroutine NSE_AddBurgers_PerVolume_Z(is, nx, ny, nz, s, result, tmp1, rhou_in, rhou_out)
        use TLab_Pointers_2D, only: pxy_wrk3d
        integer, intent(in) :: is                       ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny, nz)
        real(wp), intent(inout) :: result(nx*ny, nz)
        real(wp), intent(inout) :: tmp1(nx*ny, nz)
        real(wp), intent(in), optional :: rhou_in(nx*ny, nz)
        real(wp), intent(out), optional :: rhou_out(nx*ny, nz)

        ! -------------------------------------------------------------------
        integer(wi) k, nlines

        ! ###################################################################
        if (z%size == 1) then ! Set to zero in 2D case nx*ny
            return
        end if

        nlines = nx*ny

        call FDM_Der1_Solve(nlines, g(3)%der1, g(3)%der1%lu, s, wrk3d, wrk2d)
        call FDM_Der2_Solve(nlines, g(3)%der2, fdmDiffusion(3)%lu(:, :, is), s, tmp1, wrk3d, wrk2d)

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            if (subsidenceProps%type == TYPE_SUB_CONSTANT) then
                do k = 1, nz
                    result(:, k) = result(:, k) + tmp1(:, k) + (rho_wbackground(k) - rhou_in(:, k))*pxy_wrk3d(:, k)
                end do
            else
                result(:, :) = result(:, :) + tmp1(:, :) - rhou_in(:, :)*pxy_wrk3d(:, :)
            end if

        else                            ! Only used in anelastic formulation
            if (subsidenceProps%type == TYPE_SUB_CONSTANT) then
                do k = 1, nz
                    rhou_out(:, k) = s(:, k)*rbackground(k)
                    result(:, k) = result(:, k) + tmp1(:, k) + (rho_wbackground(k) - rhou_out(:, k))*pxy_wrk3d(:, k)
                end do
            else
                do k = 1, nz
                    rhou_out(:, k) = s(:, k)*rbackground(k)
                    result(:, k) = result(:, k) + tmp1(:, k) - rhou_out(:, k)*pxy_wrk3d(:, k)
                end do
            end if

        end if

        return
    end subroutine NSE_AddBurgers_PerVolume_Z

    !########################################################################
    !########################################################################
    subroutine NSE_Burgers_1D(nlines, nsize, der1, der2, rhou, rho_xy)
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(inout) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)
        real(wp), intent(in) :: rho_xy(nlines)

        integer(wi) ij

        do ij = 1, nsize
            rhou(:, ij) = rhou(:, ij)*rho_xy(:)
            der2(:, ij) = der2(:, ij) - rhou(:, ij)*der1(:, ij)
        end do

        return
    end subroutine NSE_Burgers_1D

    !########################################################################
    !########################################################################
#ifdef USE_MPI
    subroutine NSE_Burgers_1D_Split(nlines, nsize, der1, der2, rhou, rho_xy, diff)
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(inout) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)
        real(wp), intent(in) :: rho_xy(nlines)
        real(wp), intent(in) :: diff

        integer(wi) ij

        do ij = 1, nsize
            rhou(:, ij) = rhou(:, ij)*rho_xy(:)
            der2(:, ij) = der2(:, ij)*diff - rhou(:, ij)*der1(:, ij)
        end do

        return
    end subroutine NSE_Burgers_1D_Split
#endif

end module NSE_Burgers
