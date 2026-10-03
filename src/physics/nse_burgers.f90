! Calculate the non-linear operator N(u)(s) = dyn_visc* d^2/dx^2 s - rho u d/dx s
! The handling of the density in the anelastic case requires special treatment (see below in the code).
!
! Calculate the non-linear operator N(u)(s) = dyn_visc* d^2/dx^2 s + ((rho u)_background- rho u) d/dx s (add subsidence term)
!

#include "tlab_error.h"

module NSE_Burgers
    use TLab_Constants, only: wp, wi
    use TLab_Arrays, only: wrk3d
    use TLab_Transpose
#ifdef USE_MPI
    use TLabMPI_VARS, only: xMpi, yMpi
#endif
    use TLab_Grid, only: x, y, z
    use FDM, only: fdm_der1_X, fdm_der1_Y, fdm_der1_Z
    use FDM, only: fdm_der2_X, fdm_der2_Y, fdm_der2_Z
    ! use FDM_Derivative_1order, only: der1_periodic
    ! use FDM_Derivative_2order, only: der2_extended_periodic
!     use FDM_Derivative_Burgers
#ifdef USE_MPI
!     use FDM_Derivative_MPISplit, only: der_burgers_mpisplit
    use OPR_Partial, only: der_mode_i, der_mode_j, TYPE_TRANSPOSE, TYPE_SPLIT
    use OPR_Partial, only: fdm_der1_X_split, fdm_der2_X_split, fdm_der1_Y_split, fdm_der2_Y_split
    use OPR_Partial, only: halo_m, halo_p
#endif
    use OPR_Burgers
    implicit none
    private

    public :: NSE_Burgers_Initialize
    public :: NSE_AddBurgers_PerVolume_X
    public :: NSE_AddBurgers_PerVolume_X_Serial_Dev
    public :: NSE_AddBurgers_PerVolume_Y
    public :: NSE_AddBurgers_PerVolume_Z

    ! -----------------------------------------------------------------------
    procedure(nse_burgers_ice) :: NSE_AddBurgers_PerVolume_dt
    abstract interface
        subroutine nse_burgers_ice(is, nx, ny, nz, s, rhs, tmp2, tmp1, rhou_in)
            use TLab_Constants, only: wi, wp
            integer, intent(in) :: is                           ! scalar index; if 0, then velocity
            integer(wi), intent(in) :: nx, ny, nz
            real(wp), intent(in) :: s(nx*ny*nz)
            real(wp), intent(inout) :: rhs(nx*ny*nz)
            real(wp), intent(inout) :: tmp2(nx*ny*nz)
            real(wp), intent(inout) :: tmp1(nx*ny*nz)           ! transposed field s times density
            real(wp), intent(in), optional :: rhou_in(nx*ny*nz) ! transposed field u times density
        end subroutine
    end interface
    procedure(NSE_AddBurgers_PerVolume_dt), pointer :: NSE_AddBurgers_PerVolume_X
    procedure(NSE_AddBurgers_PerVolume_dt), pointer :: NSE_AddBurgers_PerVolume_Y

    ! -----------------------------------------------------------------------
! #ifdef USE_MPI
!     type(der_burgers_mpisplit) :: fdm_burgersX_split, fdm_burgersY_split
! #endif
!     type(der_burgers) :: fdm_burgersX, fdm_burgersY

    ! -----------------------------------------------------------------------
    class(burgers_dt), allocatable :: burgers1d_X(:), burgers1d_Y(:), burgers1d_Z(:)

    ! -----------------------------------------------------------------------
    ! Cache blocking information; tuned to levante in a subdomain 48x32x768
    integer, parameter :: groupSizeX = 1024, groupSizeY = 512, groupSizeZ = 32

contains
    !########################################################################
    !########################################################################
    subroutine NSE_Burgers_Initialize(inifile)
        use TLab_Constants, only: efile
        use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
        use TLab_Memory, only: inb_scal
        use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC, DNS_EQNS_BOUSSINESQ
        use NavierStokes, only: visc, schmidt
        use Thermo_Anelastic, only: rbackground
        use LargeScaleForcing, only: subsidenceProps, TYPE_SUB_CONSTANT, wbackground

        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        character(len=32) bakfile

        integer(wi) is
        real(wp) :: diffusivity

        ! ###################################################################
        ! Read input data
        bakfile = trim(adjustl(inifile))//'.bak'

        ! ###################################################################
        select case (nse_eqns)
        case (DNS_EQNS_ANELASTIC)
            allocate (burgers_XY :: burgers1d_X(0:inb_scal))
            allocate (burgers_XY :: burgers1d_Y(0:inb_scal))
            if (subsidenceProps%type == TYPE_SUB_CONSTANT) then
                allocate (burgers_background_Z :: burgers1d_Z(0:inb_scal))
            else
                allocate (burgers_Z :: burgers1d_Z(0:inb_scal))
            end if

            do is = 0, inb_scal     ! is = 0 corresponds to velocity fields
                if (is == 0) then
                    diffusivity = visc
                else
                    diffusivity = visc/schmidt(is)
                end if
                call burgers1d_X(is)%initialize_setrho(diffusivity, 'x', rbackground)
                call burgers1d_Y(is)%initialize_setrho(diffusivity, 'y', rbackground)
                if (subsidenceProps%type == TYPE_SUB_CONSTANT) then
                    call burgers1d_Z(is)%initialize_setrho(diffusivity, 'z', rbackground, wbackground=wbackground)
                else
                    call burgers1d_Z(is)%initialize_setrho(diffusivity, 'z', rbackground)
                end if

            end do

        case (DNS_EQNS_BOUSSINESQ)
            allocate (burgers :: burgers1d_X(0:inb_scal))
            allocate (burgers :: burgers1d_Y(0:inb_scal))
            if (subsidenceProps%type == TYPE_SUB_CONSTANT) then
                call TLab_Write_ASCII(efile, __FILE__//'. Subsidence in boussinesq not yet implemented.')
                call TLab_Stop(DNS_ERROR_UNDEVELOP)
            else
                allocate (burgers :: burgers1d_Z(0:inb_scal))
            end if

            do is = 0, inb_scal     ! is = 0 corresponds to velocity fields
                if (is == 0) then
                    diffusivity = visc
                else
                    diffusivity = visc/schmidt(is)
                end if
                call burgers1d_X(is)%initialize(diffusivity)
                call burgers1d_Y(is)%initialize(diffusivity)
                if (subsidenceProps%type == TYPE_SUB_CONSTANT) then
                    ! call burgers1d_Z(is)%initialize(diffusivity, wbackground=wbackground)
                else
                    call burgers1d_Z(is)%initialize(diffusivity)
                end if

            end do

        end select

        ! ###################################################################
        ! Setting procedure pointers
#ifdef USE_MPI
        if (xMpi%num_processors > 1) then
            select case (der_mode_i)
            case (TYPE_TRANSPOSE)
                NSE_AddBurgers_PerVolume_X => NSE_AddBurgers_PerVolume_X_MPITranspose
            case (TYPE_SPLIT)
                ! NSE_AddBurgers_PerVolume_X => NSE_AddBurgers_PerVolume_X_MPISplit
                NSE_AddBurgers_PerVolume_X => NSE_AddBurgers_PerVolume_X_Serial
                ! call fdm_burgersX_split%initialize(fdm_der1_X_split, fdm_der2_X_split)
            end select
        else
#endif
            NSE_AddBurgers_PerVolume_X => NSE_AddBurgers_PerVolume_X_Serial
            ! select type (fdm_der1_X)
            ! type is (der1_periodic)
            !     select type (fdm_der2_X)
            !     type is (der2_extended_periodic)
            !         call fdm_burgersX%initialize(fdm_der1_X, fdm_der2_X%der2)
            !     end select
            ! end select

#ifdef USE_MPI
        end if
#endif

#ifdef USE_MPI
        if (yMpi%num_processors > 1) then
            select case (der_mode_j)
            case (TYPE_TRANSPOSE)
                NSE_AddBurgers_PerVolume_Y => NSE_AddBurgers_PerVolume_Y_MPITranspose
            case (TYPE_SPLIT)
                ! NSE_AddBurgers_PerVolume_Y => NSE_AddBurgers_PerVolume_Y_MPISplit
                NSE_AddBurgers_PerVolume_Y => NSE_AddBurgers_PerVolume_Y_Serial
                ! call fdm_burgersY_split%initialize(fdm_der1_Y_split, fdm_der2_Y_split)
            end select
        else
#endif
            NSE_AddBurgers_PerVolume_Y => NSE_AddBurgers_PerVolume_Y_Serial
            ! select type (fdm_der1_Y)
            ! type is (der1_periodic)
            !     select type (fdm_der2_Y)
            !     type is (der2_extended_periodic)
            !         call fdm_burgersY%initialize(fdm_der1_Y, fdm_der2_Y%der2)
            !     end select
            ! end select

#ifdef USE_MPI
        end if
#endif
        return
    end subroutine NSE_Burgers_Initialize

    !########################################################################
    !########################################################################
    subroutine NSE_AddBurgers_PerVolume_X_Serial(is, nx, ny, nz, s, rhs, tmp2, tmp1, rhou_in)
#ifdef USE_MPI
        use TLabMPI_PROCS, only: TLabMPI_Halos_X
#endif
        integer, intent(in) :: is                           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(inout) :: rhs(nx*ny*nz)
        real(wp), intent(inout) :: tmp2(nx*ny*nz)
        real(wp), intent(inout) :: tmp1(nx*ny*nz)           ! transposed field s times density
        real(wp), intent(in), optional :: rhou_in(nx*ny*nz) ! transposed field u times density

        ! -------------------------------------------------------------------
        integer(wi) nlines
#ifdef USE_MPI
        integer np, np1, np2
#endif

        ! ###################################################################
        if (x%size == 1) then ! Set to zero in 2D case
            return
        end if

        ! Transposition: make x-direction the last one
#ifdef USE_ESSL
        call DGETMO(s, nx, nx, ny*nz, tmp1, ny*nz)
#else
        call TLab_Transpose_Real(s, nx, ny*nz, nx, tmp1, ny*nz, locBlock=trans_x_forward)
#endif

        nlines = ny*nz

#ifdef USE_MPI
        np1 = size(fdm_der1_X_split%rhs, 2)/2
        np2 = size(fdm_der2_X_split%rhs, 2)/2
        np = max(np1, np2)
        call TLabMPI_Halos_X(tmp1, nlines, np, halo_m, halo_p)

        call fdm_der1_X_split%compute(nlines, tmp1, halo_m(nlines*(np - np1) + 1:), halo_p, tmp2)
        call fdm_der2_X_split%compute(nlines, tmp1, halo_m(nlines*(np - np2) + 1:), halo_p, wrk3d)
        ! call fdm_burgersX_split%compute(nlines, tmp1, halo_m(1:np*nlines), halo_p, tmp2, wrk3d)

#else
        call fdm_der1_X%compute(nlines, tmp1, tmp2)
        call fdm_der2_X%compute(nlines, tmp1, wrk3d, tmp2)
        ! call fdm_burgersX%compute(nlines, tmp1, tmp2, wrk3d)

#endif
        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            call burgers1d_X(is)%compute(nlines, nx, der1=tmp2, der2=wrk3d, rhou=rhou_in)
        else
            call burgers1d_X(is)%compute_setrhou(nlines, nx, der1=tmp2, der2=wrk3d, rhou=tmp1)
        end if

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(wrk3d, ny*nz, ny*nz, nx, tmp2, nx)
        rhs = rhs + tmp2
#else
        call TLab_AddTranspose(wrk3d, ny*nz, nx, ny*nz, rhs, nx, locBlock=trans_x_backward)
#endif

        return
    end subroutine NSE_AddBurgers_PerVolume_X_Serial

    !########################################################################
    !########################################################################
    subroutine NSE_AddBurgers_PerVolume_X_Serial_Dev(is, nx, ny, nz, s, rhs, tmp2, tmp1, rhou_in)
#ifdef USE_MPI
        use TLabMPI_PROCS, only: TLabMPI_Halos_X
#endif
        integer, intent(in) :: is                           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx, ny*nz)
        real(wp), intent(inout) :: rhs(nx, ny*nz)
        real(wp), intent(inout) :: tmp2(groupSizeX*nx, *)
        real(wp), intent(inout) :: tmp1(groupSizeX*nx, *)           ! transposed field s times density
        real(wp), intent(in), optional :: rhou_in(groupSizeX*nx, *) ! transposed field u times density

        ! -------------------------------------------------------------------
        integer(wi) nlines
        integer(wi) ib, ip
#ifdef USE_MPI
        integer np, np1, np2
#endif

        ! ###################################################################
        if (x%size == 1) then ! Set to zero in 2D case
            return
        end if

        do ib = 1, (nz*ny - 1)/groupSizeX + 1
            ip = (ib - 1)*groupSizeX + 1
            nlines = min(groupSizeX, nz*ny - (ib - 1)*groupSizeX)

            ! Transposition: make x-direction the last one
#ifdef USE_ESSL
            call DGETMO(s(1, ip), nx, nx, nlines, tmp1(1, ib), nlines)
#else
            call TLab_Transpose_Real(s(1, ip), nx, nlines, nx, tmp1(1, ib), nlines, locBlock=trans_x_forward)
#endif

#ifdef USE_MPI
            np1 = size(fdm_der1_X_split%rhs, 2)/2
            np2 = size(fdm_der2_X_split%rhs, 2)/2
            np = max(np1, np2)
            call TLabMPI_Halos_X(tmp1(:, ib), nlines, np, halo_m, halo_p)

            call fdm_der1_X_split%compute(nlines, tmp1(1, ib), halo_m(nlines*(np - np1) + 1:), halo_p, tmp2)
            call fdm_der2_X_split%compute(nlines, tmp1(1, ib), halo_m(nlines*(np - np2) + 1:), halo_p, wrk3d)
            ! call fdm_burgersX_split%compute(nlines, tmp1(1, ib), halo_m, halo_p, tmp2, wrk3d)

#else
            call fdm_der1_X%compute(nlines, tmp1(1, ib), tmp2)
            call fdm_der2_X%compute(nlines, tmp1(1, ib), wrk3d, tmp2)
            ! call fdm_burgersX%compute(nlines, tmp1(1, ib), tmp2, wrk3d)

#endif
            if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
                call burgers1d_X(is)%compute(nlines, nx, der1=tmp2, der2=wrk3d, rhou=rhou_in(1, ib))
            else
                burgers1d_X(is)%offset = ip - 1
                call burgers1d_X(is)%compute_setrhou(nlines, nx, der1=tmp2, der2=wrk3d, rhou=tmp1(1, ib))
                burgers1d_X(is)%offset = 0
            end if

            ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
            call DGETMO(wrk3d, ny*nz, ny*nz, nx, tmp2, nx)
            rhs(ip:ip + nlines*nx - 1, 1) = rhs(ip:ip + nlines*nx - 1) + tmp2(1:nlines*nx)
#else
            call TLab_AddTranspose(wrk3d, nlines, nx, nlines, rhs(1, ip), nx, locBlock=trans_x_backward)
#endif
        end do

        return
    end subroutine NSE_AddBurgers_PerVolume_X_Serial_Dev

    !########################################################################
    !########################################################################
#ifdef USE_MPI
    subroutine NSE_AddBurgers_PerVolume_X_MPITranspose(is, nx, ny, nz, s, rhs, tmp2, tmp1, rhou_in)
        use TLabMPI_Transpose_DerivedTypes, only: tmpi_trp_X
        integer, intent(in) :: is                           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(inout) :: rhs(nx*ny*nz)
        real(wp), intent(inout) :: tmp2(nx*ny*nz)
        real(wp), intent(inout) :: tmp1(nx*ny*nz)           ! transposed field s times density
        real(wp), intent(in), optional :: rhou_in(nx*ny*nz) ! transposed field u times density

        ! -------------------------------------------------------------------
        integer(wi) nlines

        ! ###################################################################
        if (x%size == 1) then ! Set to zero in 2D case
            return
        end if

        nlines = tmpi_trp_X%nlines

        ! Transposition: make x-direction the last one
        call tmpi_trp_X%forward(s, tmp2)
#ifdef USE_ESSL
        call DGETMO(tmp2, x%size, x%size, nlines, tmp1, nlines)
#else
        call TLab_Transpose_Real(tmp2, x%size, nlines, x%size, tmp1, nlines)
#endif

        call fdm_der1_X%compute(nlines, tmp1, wrk3d)
        call fdm_der2_X%compute(nlines, tmp1, tmp2, wrk3d)

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            call burgers1d_X(is)%compute(nlines, nx*xMpi%num_processors, der1=wrk3d, der2=tmp2, rhou=rhou_in)
        else
            call burgers1d_X(is)%compute_setrhou(nlines, nx*xMpi%num_processors, der1=wrk3d, der2=tmp2, rhou=tmp1)
        end if

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(tmp2, nlines, nlines, x%size, wrk3d, x%size)
#else
        call TLab_Transpose_Real(tmp2, nlines, x%size, nlines, wrk3d, x%size)
#endif
        call tmpi_trp_X%backward(wrk3d, tmp2)
        rhs = rhs + tmp2

        return
    end subroutine NSE_AddBurgers_PerVolume_X_MPITranspose

#endif

    !########################################################################
    !########################################################################
    subroutine NSE_AddBurgers_PerVolume_Y_Serial(is, nx, ny, nz, s, rhs, tmp2, tmp1, rhou_in)
#ifdef USE_MPI
        use TLabMPI_PROCS, only: TLabMPI_Halos_Y
#endif
        integer, intent(in) :: is                           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(inout) :: rhs(nx*ny*nz)
        real(wp), intent(inout) :: tmp2(nx*ny*nz)
        real(wp), intent(inout) :: tmp1(nx*ny*nz)           ! transposed field s times density
        real(wp), intent(in), optional :: rhou_in(nx*ny*nz) ! transposed field u times density

        ! -------------------------------------------------------------------
        integer(wi) nlines
        integer(wi) ib, i, k
#ifdef USE_MPI
        integer np, np1, np2
#endif

        ! ###################################################################
        if (y%size == 1) then ! Set to zero in 2D case
            return
        end if

        ib = 1
        do i = 1, nx
            do k = 1, nz, groupSizeY
                nlines = min(groupSizeY, nz - k + 1)

                ! memory arrangement; reduce
                call reduce_y(s, nx, ny, nz, i, k, nlines, tmp1(ib))

#ifdef USE_MPI
                np1 = size(fdm_der1_Y_split%rhs, 2)/2
                np2 = size(fdm_der2_Y_split%rhs, 2)/2
                np = max(np1, np2)
                call TLabMPI_Halos_Y(tmp1(ib:ib + nlines*ny - 1), nlines, np, halo_m, halo_p)

                call fdm_der1_Y_split%compute(nlines, tmp1(ib), halo_m(nlines*(np - np1) + 1:), halo_p, tmp2)
                call fdm_der2_Y_split%compute(nlines, tmp1(ib), halo_m(nlines*(np - np2) + 1:), halo_p, wrk3d)
                ! call fdm_burgersY_split%compute(nlines, tmp1(1, ib), halo_m, halo_p, tmp2, wrk3d)
#else
                call fdm_der1_Y%compute(nlines, tmp1(ib), tmp2)
                call fdm_der2_Y%compute(nlines, tmp1(ib), wrk3d, tmp2)
                ! call fdm_burgersY%compute(nlines, tmp1(1, ib), tmp2, wrk3d)
#endif

                if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
                    call burgers1d_Y(is)%compute(nlines, ny, der1=tmp2, der2=wrk3d, rhou=rhou_in(ib))
                else
                    burgers1d_Y(is)%offset = (i - 1)*nz + k - 1
                    call burgers1d_Y(is)%compute_setrhou(nlines, ny, der1=tmp2, der2=wrk3d, rhou=tmp1(ib))
                    burgers1d_Y(is)%offset = 0
                end if

                ! memory arrangement; spread_add
                call spread_add_y(wrk3d, nx, ny, nz, i, k, nlines, rhs)

                ib = ib + nlines*ny

            end do
        end do

        return
    end subroutine NSE_AddBurgers_PerVolume_Y_Serial

    subroutine reduce_y(a, nx, ny, nz, i, k, nlines, b)
        real(wp), intent(in) :: a(nx, ny, nz)
        integer, intent(in) :: nx, ny, nz, i, k, nlines
        real(wp), intent(out) :: b(*)

        integer ip, jj, kk

        ip = 0
        do jj = 1, ny
            do kk = k, k + nlines - 1
                ip = ip + 1
                b(ip) = a(i, jj, kk)
            end do
        end do

        return
    end subroutine

    subroutine spread_add_y(a, nx, ny, nz, i, k, nlines, b)
        real(wp), intent(in) :: a(*)
        integer, intent(in) :: nx, ny, nz, i, k, nlines
        real(wp), intent(inout) :: b(nx, ny, nz)

        integer ip, jj, kk

        ip = 0
        do jj = 1, ny
            do kk = k, k + nlines - 1
                ip = ip + 1
                b(i, jj, kk) = b(i, jj, kk) + a(ip)
            end do
        end do

        return
    end subroutine

    !########################################################################
    !########################################################################
#ifdef USE_MPI
    subroutine NSE_AddBurgers_PerVolume_Y_MPITranspose(is, nx, ny, nz, s, rhs, tmp2, tmp1, rhou_in)
        use TLabMPI_Transpose_DerivedTypes, only: tmpi_trp_Y
        integer, intent(in) :: is                           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(inout) :: rhs(nx*ny*nz)
        real(wp), intent(inout) :: tmp2(nx*ny*nz)
        real(wp), intent(inout) :: tmp1(nx*ny*nz)           ! transposed field s times density
        real(wp), intent(in), optional :: rhou_in(nx*ny*nz) ! transposed field u times density

        ! -------------------------------------------------------------------
        integer(wi) nlines

        ! ###################################################################
        if (y%size == 1) then ! Set to zero in 2D case
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
        call fdm_der2_Y%compute(nlines, tmp1, tmp2, wrk3d)

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            call burgers1d_Y(is)%compute(nlines, ny*yMpi%num_processors, der1=wrk3d, der2=tmp2, rhou=rhou_in)
        else
            call burgers1d_Y(is)%compute_setrhou(nlines, ny*yMpi%num_processors, der1=wrk3d, der2=tmp2, rhou=tmp1)
        end if

        ! Put arrays back in the order in which they came in
        call tmpi_trp_Y%backward(tmp2, wrk3d)
#ifdef USE_ESSL
        call DGETMO(wrk3d, nz, nz, nx*ny, tmp2, nx*ny)
        rhs = rhs + tmp2
#else
        call TLab_AddTranspose(wrk3d, nz, nx*ny, nz, rhs, nx*ny)
#endif

        return
    end subroutine NSE_AddBurgers_PerVolume_Y_MPITranspose

#endif

    !########################################################################
    !########################################################################
    subroutine NSE_AddBurgers_PerVolume_Z(is, nx, ny, nz, s, rhs, tmp2, tmp1, rhou_in)
        integer, intent(in) :: is                       ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny, nz)
        real(wp), intent(inout) :: rhs(nx*ny, nz)
        real(wp), intent(inout) :: tmp2(groupSizeZ*nz, *) !nx*ny/groupSizeZ)
        real(wp), intent(inout) :: tmp1(groupSizeZ*nz, *) !nx*ny/groupSizeZ)
        real(wp), intent(in), optional :: rhou_in(groupSizeZ*nz, *) !nx*ny/groupSizeZ)

        ! -------------------------------------------------------------------
        integer(wi) nlines
        integer(wi) ib, ip

        ! ###################################################################
        if (z%size == 1) then ! Set to zero in 2D case nx*ny
            return
        end if

        do ib = 1, (nx*ny - 1)/groupSizeZ + 1
            ip = (ib - 1)*groupSizeZ + 1
            nlines = min(groupSizeZ, nx*ny - (ib - 1)*groupSizeZ)

            ! memory arrangement
            call reduce_z(s(ip, 1), nlines, nx*ny, nz, tmp1(1, ib))

            call fdm_der1_Z%compute(nlines, tmp1(1, ib), wrk3d)
            call fdm_der2_Z%compute(nlines, tmp1(1, ib), tmp2, wrk3d)

            if (present(rhou_in)) then      ! velocity (times density) is passed as argument
                call burgers1d_Z(is)%compute(nlines, nz, der1=wrk3d, der2=tmp2, rhou=rhou_in(1, ib))
            else
                call burgers1d_Z(is)%compute_setrhou(nlines, nz, der1=wrk3d, der2=tmp2, rhou=tmp1(1, ib))
            end if

            ! memory arrangement
            call spread_add_z(tmp2, nlines, nx*ny, nz, rhs(ip, 1))

        end do

        return
    end subroutine

    subroutine reduce_z(a, nlines, nca, nmax, b)
        integer(wi), intent(in) :: nlines, nca, nmax
        real(wp), intent(in) :: a(nca, *)
        real(wp), intent(out) :: b(nlines, nmax)

        integer n

        do n = 1, nmax
            b(1:nlines, n) = a(1:nlines, n)
        end do

        return
    end subroutine

    subroutine spread_add_z(a, nlines, ncb, nmax, b)
        integer(wi), intent(in) :: nlines, ncb, nmax
        real(wp), intent(in) :: a(nlines, nmax)
        real(wp), intent(out) :: b(ncb, *)

        integer n

        do n = 1, nmax
            b(1:nlines, n) = b(1:nlines, n) + a(1:nlines, n)
        end do

        return
    end subroutine

end module NSE_Burgers
