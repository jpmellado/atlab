! Calculate the non-linear operator N(u)(s) = dyn_visc* d^2/dx^2 s - rho u d/dx s

module NSE_Burgers
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: MAX_VARS
    use TLab_Constants, only: efile, lfile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Arrays, only: wrk3d
    use TLab_Transpose
#ifdef USE_MPI
    use TLabMPI_VARS, only: xMpi, yMpi
    use TLabMPI_Transpose_DerivedTypes, only: tmpi_trp_X, tmpi_trp_Y
    use FDM_Derivative_MPISplit
#endif
    use TLab_Grid, only: x, y, z
    use FDM, only: fdm_der1_X, fdm_der1_Y, fdm_der1_Z
    use FDM, only: fdm_der2_X, fdm_der2_Y, fdm_der2_Z
    use FDM_Derivative_1order, only: der1_periodic
    use FDM_Derivative_2order, only: der2_extended_periodic
    use FDM_Derivative_Burgers
    use OPR_Burgers
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
    procedure(NSE_AddBurgers_PerVolume_dt), pointer :: NSE_AddBurgers_PerVolume_X
    procedure(NSE_AddBurgers_PerVolume_dt), pointer :: NSE_AddBurgers_PerVolume_Y

    ! -----------------------------------------------------------------------
#ifdef USE_MPI
    type(der_burgers_mpisplit) :: fdm_burgersX_split, fdm_burgersY_split
#endif
    type(der_burgers) :: fdm_burgersX, fdm_burgersY

    ! -----------------------------------------------------------------------
    real(wp) :: diffusivity(0:MAX_VARS)

    real(wp), allocatable :: rho_wbackground(:)     ! subsidence velocity (times density)

    class(burgers1d), allocatable :: burgers1d_X(:), burgers1d_Y(:), burgers1d_Z(:)

contains
    !########################################################################
    !########################################################################
    subroutine NSE_Burgers_Initialize(inifile)
        use TLab_Memory, only: TLab_Allocate_Real
        use TLab_Memory, only: inb_scal
        use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC, DNS_EQNS_BOUSSINESQ
        use NavierStokes, only: visc, schmidt
        use Thermo_Anelastic, only: rbackground

        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        character(len=32) bakfile

        integer(wi) is

        ! ###################################################################
        ! Read input data
        bakfile = trim(adjustl(inifile))//'.bak'

        ! ###################################################################
        select case (nse_eqns)
        case (DNS_EQNS_ANELASTIC)
            allocate (burgers1d_anelastic :: burgers1d_X(0:inb_scal))
            allocate (burgers1d_anelastic :: burgers1d_Y(0:inb_scal))
            if (subsidenceProps%type == TYPE_SUB_CONSTANT) then
                allocate (burgers1d_subsidence_anelastic :: burgers1d_Z(0:inb_scal))
            else
                allocate (burgers1d_anelastic :: burgers1d_Z(0:inb_scal))
            end if

            do is = 0, inb_scal     ! is = 0 corresponds to velocity fields
                if (is == 0) then
                    diffusivity(is) = visc
                else
                    diffusivity(is) = visc/schmidt(is)
                end if
                call burgers1d_X(is)%initialize(diffusivity(is), 'x', rbackground)
                call burgers1d_Y(is)%initialize(diffusivity(is), 'y', rbackground)
                if (subsidenceProps%type == TYPE_SUB_CONSTANT) then
                    call burgers1d_Z(is)%initialize(diffusivity(is), 'z', rbackground, wbackground=wbackground)
                else
                    call burgers1d_Z(is)%initialize(diffusivity(is), 'z', rbackground)
                end if

            end do

        case (DNS_EQNS_BOUSSINESQ)
            allocate (burgers1d_boussinesq :: burgers1d_X(0:inb_scal))
            allocate (burgers1d_boussinesq :: burgers1d_Y(0:inb_scal))
            if (subsidenceProps%type == TYPE_SUB_CONSTANT) then
                allocate (burgers1d_subsidence_boussinesq :: burgers1d_Z(0:inb_scal))
            else
                allocate (burgers1d_boussinesq :: burgers1d_Z(0:inb_scal))
            end if

            do is = 0, inb_scal     ! is = 0 corresponds to velocity fields
                if (is == 0) then
                    diffusivity(is) = visc
                else
                    diffusivity(is) = visc/schmidt(is)
                end if
                call burgers1d_X(is)%initialize(diffusivity(is))
                call burgers1d_Y(is)%initialize(diffusivity(is))
                if (subsidenceProps%type == TYPE_SUB_CONSTANT) then
                    call burgers1d_Z(is)%initialize(diffusivity(is), wbackground=wbackground)
                else
                    call burgers1d_Z(is)%initialize(diffusivity(is))
                end if

            end do

        end select

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
        if (xMpi%num_processors > 1) then
            select case (der_mode_i)
            case (TYPE_TRANSPOSE)
                NSE_AddBurgers_PerVolume_X => NSE_AddBurgers_PerVolume_X_MPITranspose
            case (TYPE_SPLIT)
                NSE_AddBurgers_PerVolume_X => NSE_AddBurgers_PerVolume_X_MPISplit
                call fdm_burgersX_split%initialize(fdm_der1_X_split, fdm_der2_X_split)
            end select
        else
#endif
            NSE_AddBurgers_PerVolume_X => NSE_AddBurgers_PerVolume_X_Serial
            select type (fdm_der1_X)
            type is (der1_periodic)
                select type (fdm_der2_X)
                type is (der2_extended_periodic)
                    call fdm_burgersX%initialize(fdm_der1_X, fdm_der2_X%der2)
                end select
            end select

#ifdef USE_MPI
        end if
#endif

#ifdef USE_MPI
        if (yMpi%num_processors > 1) then
            select case (der_mode_j)
            case (TYPE_TRANSPOSE)
                NSE_AddBurgers_PerVolume_Y => NSE_AddBurgers_PerVolume_Y_MPITranspose
            case (TYPE_SPLIT)
                NSE_AddBurgers_PerVolume_Y => NSE_AddBurgers_PerVolume_Y_MPISplit
                call fdm_burgersY_split%initialize(fdm_der1_Y_split, fdm_der2_Y_split)
            end select
        else
#endif
            NSE_AddBurgers_PerVolume_Y => NSE_AddBurgers_PerVolume_Y_Serial
            select type (fdm_der1_Y)
            type is (der1_periodic)
                select type (fdm_der2_Y)
                type is (der2_extended_periodic)
                    call fdm_burgersY%initialize(fdm_der1_Y, fdm_der2_Y%der2)
                end select
            end select

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
        call TLab_Transpose_Real(s, nx, ny*nz, nx, tmp1, ny*nz, locBlock=trans_x_forward)
#endif

        nlines = ny*nz

        call fdm_der1_X%compute(nlines, tmp1, result)
        call fdm_der2_X%compute(nlines, tmp1, wrk3d, result)
        ! call fdm_burgersX%compute(nlines, tmp1, result, wrk3d)

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            call burgers1d_X(is)%compute(nlines, nx, der1=result, der2=wrk3d, rhou=rhou_in)
        else
            call burgers1d_X(is)%compute_setrhou(nlines, nx, der1=result, der2=wrk3d, rhou=tmp1)
        end if

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(wrk3d, ny*nz, ny*nz, nx, result, nx)
        rhs = rhs + result
#else
        call TLab_AddTranspose(wrk3d, ny*nz, nx, ny*nz, rhs, nx, locBlock=trans_x_backward)
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

        nlines = tmpi_trp_X%nlines

        ! Transposition: make x-direction the last one
        call tmpi_trp_X%forward(s, result)
#ifdef USE_ESSL
        call DGETMO(result, x%size, x%size, nlines, tmp1, nlines)
#else
        call TLab_Transpose_Real(result, x%size, nlines, x%size, tmp1, nlines)
#endif

        call fdm_der1_X%compute(nlines, tmp1, wrk3d)
        call fdm_der2_X%compute(nlines, tmp1, result, wrk3d)

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            call burgers1d_X(is)%compute(nlines, nx*xMpi%num_processors, der1=wrk3d, der2=result, rhou=rhou_in)
        else
            call burgers1d_X(is)%compute_setrhou(nlines, nx*xMpi%num_processors, der1=wrk3d, der2=result, rhou=tmp1)
        end if

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(result, nlines, nlines, x%size, wrk3d, x%size)
#else
        call TLab_Transpose_Real(result, nlines, x%size, nlines, wrk3d, x%size)
#endif
        call tmpi_trp_X%backward(wrk3d, result)
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
        call TLab_Transpose_Real(s, nx, ny*nz, nx, tmp1, ny*nz, locBlock=trans_x_forward)
#endif

        nlines = ny*nz

        np1 = size(fdm_der1_X_split%rhs, 2)/2
        np2 = size(fdm_der2_X_split%rhs, 2)/2
        np = max(np1, np2)
        call TLabMPI_Halos_X(tmp1, nlines, np, pyz_halo_m(:, 1), pyz_halo_p(:, 1))

        call fdm_der1_X_split%compute(nlines, tmp1, pyz_halo_m(:, np - np1 + 1:np), pyz_halo_p, result)
        call fdm_der2_X_split%compute(nlines, tmp1, pyz_halo_m(:, np - np2 + 1:np), pyz_halo_p, wrk3d)
        ! call fdm_burgersX_split%compute(nlines, tmp1, pyz_halo_m(:, 1:np), pyz_halo_p(:, 1:np), result, wrk3d)

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            call burgers1d_X(is)%compute(nlines, nx, der1=result, der2=wrk3d, rhou=rhou_in)
        else
            call burgers1d_X(is)%compute_setrhou(nlines, nx, der1=result, der2=wrk3d, rhou=tmp1)
        end if

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(wrk3d, ny*nz, ny*nz, nx, result, nx)
        rhs = rhs + result
#else
        call TLab_AddTranspose(wrk3d, ny*nz, nx, ny*nz, rhs, nx, locBlock=trans_x_backward)
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
        call TLab_Transpose_Real(s, nx*ny, nz, nx*ny, tmp1, nz, locBlock=trans_y_forward)
#endif

        nlines = nx*nz

        call fdm_der1_Y%compute(nlines, tmp1, result)
        call fdm_der2_Y%compute(nlines, tmp1, wrk3d, result)
        ! call fdm_burgersY%compute(nlines, tmp1, result, wrk3d)

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            call burgers1d_Y(is)%compute(nlines, ny, der1=result, der2=wrk3d, rhou=rhou_in)
        else
            call burgers1d_Y(is)%compute_setrhou(nlines, ny, der1=result, der2=wrk3d, rhou=tmp1)
        end if

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(wrk3d, nz, nz, nx*ny, result, nx*ny)
        rhs = rhs + result
#else
        call TLab_AddTranspose(wrk3d, nz, nx*ny, nz, rhs, nx*ny, locBlock=trans_y_backward)
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
        call tmpi_trp_Y%forward(wrk3d, tmp1)
        nlines = tmpi_trp_Y%nlines

        call fdm_der1_Y%compute(nlines, tmp1, wrk3d)
        call fdm_der2_Y%compute(nlines, tmp1, result, wrk3d)

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            call burgers1d_Y(is)%compute(nlines, ny*yMpi%num_processors, der1=wrk3d, der2=result, rhou=rhou_in)
        else
            call burgers1d_Y(is)%compute_setrhou(nlines, ny*yMpi%num_processors, der1=wrk3d, der2=result, rhou=tmp1)
        end if

        ! Put arrays back in the order in which they came in
        call tmpi_trp_Y%backward(result, wrk3d)
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
        call TLab_Transpose_Real(s, nx*ny, nz, nx*ny, tmp1, nz, locBlock=trans_y_forward)
#endif

        nlines = nx*nz

        np1 = size(fdm_der1_Y_split%rhs, 2)/2
        np2 = size(fdm_der2_Y_split%rhs, 2)/2
        np = max(np1, np2)
        call TLabMPI_Halos_Y(tmp1, nlines, np, pxz_halo_m(:, 1), pxz_halo_p(:, 1))

        call fdm_der1_Y_split%compute(nlines, tmp1, pxz_halo_m(:, np - np1 + 1:np), pxz_halo_p, result)
        call fdm_der2_Y_split%compute(nlines, tmp1, pxz_halo_m(:, np - np2 + 1:np), pxz_halo_p, wrk3d)
        ! call fdm_burgersY_split%compute(nlines, tmp1, pxz_halo_m(:, 1:np), pxz_halo_p(:, 1:np), result, wrk3d)

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            call burgers1d_Y(is)%compute(nlines, ny, der1=result, der2=wrk3d, rhou=rhou_in)
        else
            call burgers1d_Y(is)%compute_setrhou(nlines, ny, der1=result, der2=wrk3d, rhou=tmp1)
        end if

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(wrk3d, nz, nz, nx*ny, result, nx*ny)
        rhs = rhs + result
#else
        call TLab_AddTranspose(wrk3d, nz, nx*ny, nz, rhs, nx*ny, locBlock=trans_y_backward)
#endif

        return
    end subroutine NSE_AddBurgers_PerVolume_Y_MPISplit

#endif

    !########################################################################
    !########################################################################
    subroutine NSE_AddBurgers_PerVolume_Z(is, nx, ny, nz, s, result, tmp1, rhou_in, rhou_out)
        use TLab_Pointers_2D, only: pxy_wrk3d
        use Thermo_Anelastic, only: rbackground
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

        call fdm_der1_Z%compute(nlines, s, wrk3d)
        call fdm_der2_Z%compute(nlines, s, tmp1, wrk3d)

        if (present(rhou_in)) then      ! velocity (times density) is passed as argument
            call burgers1d_Z(is)%add(nlines, nz, der1=wrk3d, der2=tmp1, rhou=rhou_in, result=result)
            ! if (subsidenceProps%type == TYPE_SUB_CONSTANT) then
            !     do k = 1, nz
            !         result(:, k) = result(:, k) + tmp1(:, k)*diffusivity(is) + (rho_wbackground(k) - rhou_in(:, k))*pxy_wrk3d(:, k)
            !     end do
            ! else
            !     result(:, :) = result(:, :) + tmp1(:, :)*diffusivity(is) - rhou_in(:, :)*pxy_wrk3d(:, :)
            ! end if

        else                            ! Only used in anelastic formulation
            ! call burgers1d_Z(is)%add_setrhou(nlines, nz, der1=wrk3d, der2=tmp1, u=s, rhou=rhou_out, result)
            if (subsidenceProps%type == TYPE_SUB_CONSTANT) then
                do k = 1, nz
                    rhou_out(:, k) = s(:, k)*rbackground(k)
                    result(:, k) = result(:, k) + tmp1(:, k)*diffusivity(is) + (rho_wbackground(k) - rhou_out(:, k))*pxy_wrk3d(:, k)
                end do
            else
                do k = 1, nz
                    rhou_out(:, k) = s(:, k)*rbackground(k)
                    result(:, k) = result(:, k) + tmp1(:, k)*diffusivity(is) - rhou_out(:, k)*pxy_wrk3d(:, k)
                end do
            end if

        end if

        return
    end subroutine NSE_AddBurgers_PerVolume_Z

end module NSE_Burgers
