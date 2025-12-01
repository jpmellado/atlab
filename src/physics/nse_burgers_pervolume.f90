#include "tlab_error.h"

! Calculate the non-linear operator N(u)(s) = visc* d^2/dx^2 s - u d/dx s

module NSE_Burgers_PerVolume
    use TLab_Constants, only: wp, wi, efile, lfile, BCS_NONE, MAX_VARS
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLab_Arrays, only: wrk2d, wrk3d
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_j
    use TLabMPI_Transpose
#endif
    use TLab_Grid, only: x, y, z
    use FDM, only: fdm_dt, g
    use NavierStokes, only: nse_eqns, DNS_EQNS_ANELASTIC
    use Thermo_Anelastic, only: ribackground, rbackground
    use NavierStokes, only: visc, schmidt
    use OPR_Partial
    use LargeScaleForcing, only: subsidenceProps, TYPE_SUB_CONSTANT, wbackground
    implicit none
    private

    public :: NSE_Burgers_PerVolume_Initialize
    public :: NSE_Burgers_PerVolume_X
    public :: NSE_Burgers_PerVolume_Y
    public :: NSE_Burgers_PerVolume_Z

    ! -----------------------------------------------------------------------
    procedure(NSE_Burgers_PerVolume_interface) :: NSE_Burgers_PerVolume_dt
    abstract interface
        subroutine NSE_Burgers_PerVolume_interface(is, nx, ny, nz, s, result, tmp1, rhou_in)
            use TLab_Constants, only: wi, wp
            integer, intent(in) :: is                           ! scalar index; if 0, then velocity
            integer(wi), intent(in) :: nx, ny, nz
            real(wp), intent(in) :: s(nx*ny*nz)
            real(wp), intent(out) :: result(nx*ny*nz)
            real(wp), intent(out) :: tmp1(nx*ny*nz)             ! transposed field s times density
            real(wp), intent(in), optional :: rhou_in(nx*ny*nz) ! transposed field u times density
        end subroutine
    end interface
    procedure(NSE_Burgers_PerVolume_dt), pointer :: NSE_Burgers_PerVolume_X, NSE_Burgers_PerVolume_Y

    type :: fdm_diffusion_dt
        sequence
        real(wp), allocatable :: lu(:, :, :)
        ! type(rho_anelastic_dt) :: rho_anelastic
    end type fdm_diffusion_dt
    type(fdm_diffusion_dt) :: fdmDiffusion(3)

    type :: rho_anelastic_dt                        ! 1/rho in diffusion term in anelastic formulation
        sequence
        logical :: active = .false.
        real(wp), allocatable :: values(:)
    end type rho_anelastic_dt
    type(rho_anelastic_dt) :: rho_anelastic(3)      ! one for each direction

    real(wp), dimension(:), pointer :: p_vel
    real(wp) :: diffusivity(0:MAX_VARS)

contains
    !########################################################################
    !########################################################################
    subroutine NSE_Burgers_PerVolume_Initialize(inifile)
        use TLab_Memory, only: imax, jmax, kmax
        use TLab_Memory, only: TLab_Allocate_Real
        use TLab_Memory, only: inb_scal
#ifdef USE_MPI
        use TLabMPI_VARS, only: ims_pro_i, ims_npro_i, ims_pro_j, ims_npro_j
#endif

        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        character(len=32) bakfile

        integer(wi) ig, is, ip, j, ndl, idl, ic
        integer(wi) nlines, offset
        real(wp) dummy

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
                    dummy = visc
                    diffusivity(is) = visc
                else
                    dummy = visc/schmidt(is)
                    diffusivity(is) = visc/schmidt(is)
                end if

                fdmDiffusion(ig)%lu(:, :, is) = g(ig)%der2%lu(:, :)                 ! Check routines Thomas3C_LU and Thomas3C_Solve
                fdmDiffusion(ig)%lu(:, idl, is) = g(ig)%der2%lu(:, idl)*dummy
                ! if (g(ig)%periodic) then
                !     fdmDiffusion(ig)%lu(:, ndl + 1, is) = g(ig)%der2%lu(:, ndl + 1)/dummy
                ! end if

            end do
        end do

        ! ###################################################################
        ! Initialize anelastic density correction
        if (nse_eqns == DNS_EQNS_ANELASTIC) then
            call TLab_Write_ASCII(lfile, 'Initialize anelastic density correction in burgers operator.')

            ! -----------------------------------------------------------------------
            ! Density correction term in the burgers operator along X
            rho_anelastic(1)%active = .true.
#ifdef USE_MPI
            if (ims_npro_i > 1) then
                nlines = tmpi_plan_dx%nlines
                offset = nlines*ims_pro_i
            else
#endif
                nlines = jmax*kmax
                offset = 0
#ifdef USE_MPI
            end if
#endif
            allocate (rho_anelastic(1)%values(nlines))
            do j = 1, nlines
                ip = (offset + j - 1)/jmax + 1
                rho_anelastic(1)%values(j) = rbackground(ip)
            end do

            ! -----------------------------------------------------------------------
            ! Density correction term in the burgers operator along Y
            rho_anelastic(2)%active = .true.
#ifdef USE_MPI
            if (ims_npro_j > 1) then
                nlines = tmpi_plan_dy%nlines
                offset = nlines*ims_pro_j
            else
#endif
                nlines = imax*kmax
                offset = 0
#ifdef USE_MPI
            end if
#endif
            allocate (rho_anelastic(2)%values(nlines))
            do j = 1, nlines
                ip = mod(offset + j - 1, z%size) + 1
                rho_anelastic(2)%values(j) = rbackground(ip)
            end do

        end if

        ! ###################################################################
        ! Setting procedure pointers
! #ifdef USE_MPI
!         if (ims_npro_i > 1) then
!             select case (der_mode_i)
!             case (TYPE_TRANSPOSE)
!                 NSE_Burgers_PerVolume_X => NSE_Burgers_PerVolume_X_MPITranspose
!             case (TYPE_SPLIT)
!                 NSE_Burgers_PerVolume_X => NSE_Burgers_PerVolume_X_MPISplit
!             end select
!         else
! #endif
            NSE_Burgers_PerVolume_X => NSE_Burgers_PerVolume_X_Serial
! #ifdef USE_MPI
!         end if
! #endif

! #ifdef USE_MPI
!         if (ims_npro_j > 1) then
!             select case (der_mode_j)
!             case (TYPE_TRANSPOSE)
!                 NSE_Burgers_PerVolume_Y => NSE_Burgers_PerVolume_Y_MPITranspose
!             case (TYPE_SPLIT)
!                 NSE_Burgers_PerVolume_Y => NSE_Burgers_PerVolume_Y_MPISplit
!             end select
!         else
! #endif
            NSE_Burgers_PerVolume_Y => NSE_Burgers_PerVolume_Y_Serial
! #ifdef USE_MPI
!         end if
! #endif
        return
    end subroutine NSE_Burgers_PerVolume_Initialize

    !########################################################################
    !########################################################################
    subroutine NSE_Burgers_PerVolume_X_Serial(is, nx, ny, nz, s, result, tmp1, rhou_in)
        integer, intent(in) :: is                           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
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
        call TLab_Transpose(s, nx, ny*nz, nx, tmp1, ny*nz)
#endif

        nlines = ny*nz

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            call NSE_Burgers_1D(nlines, g(1), fdmDiffusion(1)%lu(:, :, is), &
                                s=tmp1, &
                                result=wrk3d, &
                                dsdx=result, &
                                rhou=rhou_in)
        else
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call NSE_Burgers_1D(nlines, g(1), fdmDiffusion(1)%lu(:, :, is), &
                                    s=tmp1, &
                                    result=wrk3d, &
                                    dsdx=result, &
                                    rhou=tmp1, &
                                    rho=rho_anelastic(1)%values(:))
            else
                call NSE_Burgers_1D(nlines, g(1), fdmDiffusion(1)%lu(:, :, is), &
                                    s=tmp1, &
                                    result=wrk3d, &
                                    dsdx=result, &
                                    rhou=tmp1)
            end if
        end if

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(wrk3d, ny*nz, ny*nz, nx, result, nx)
#else
        call TLab_Transpose(wrk3d, ny*nz, nx, ny*nz, result, nx)
#endif

        return
    end subroutine NSE_Burgers_PerVolume_X_Serial

    !########################################################################
    !########################################################################
! #ifdef USE_MPI
!     subroutine NSE_Burgers_PerVolume_X_MPITranspose(is, nx, ny, nz, s, result, tmp1, u_t)
!         integer, intent(in) :: is                       ! scalar index; if 0, then velocity
!         integer(wi), intent(in) :: nx, ny, nz
!         real(wp), intent(in) :: s(nx*ny*nz)
!         real(wp), intent(out) :: result(nx*ny*nz)
!         real(wp), intent(out), target :: tmp1(nx*ny*nz)         ! transposed field s
!         real(wp), intent(in), optional, target :: u_t(nx*ny*nz)

!         ! -------------------------------------------------------------------
!         integer(wi) nlines

!         ! ###################################################################
!         if (x%size == 1) then ! Set to zero in 2D case
!             result = 0.0_wp
!             return
!         end if

!         nlines = tmpi_plan_dx%nlines

!         ! Transposition: make x-direction the last one
!         call TLabMPI_Trp_ExecI_Forward(s, result, tmpi_plan_dx)
! #ifdef USE_ESSL
!         call DGETMO(result, g(1)%size, g(1)%size, nlines, tmp1, nlines)
! #else
!         call TLab_Transpose(result, g(1)%size, nlines, g(1)%size, tmp1, nlines)
! #endif

!         if (present(u_t)) then  ! transposed velocity is passed as argument
!             p_vel => u_t
!         else
!             p_vel => tmp1
!         end if

!         call NSE_Burgers_1D(nlines, g(1), fdmDiffusion(1)%lu(:, :, is), tmp1, p_vel, result, wrk3d)

!         ! Put arrays back in the order in which they came in
! #ifdef USE_ESSL
!         call DGETMO(result, nlines, nlines, g(1)%size, wrk3d, g(1)%size)
! #else
!         call TLab_Transpose(result, nlines, g(1)%size, nlines, wrk3d, g(1)%size)
! #endif
!         call TLabMPI_Trp_ExecI_Backward(wrk3d, result, tmpi_plan_dx)

!         return
!     end subroutine NSE_Burgers_PerVolume_X_MPITranspose

!     !########################################################################
!     !########################################################################
!     subroutine NSE_Burgers_PerVolume_X_MPISplit(is, nx, ny, nz, s, result, tmp1, u_t)
!         use TLabMPI_PROCS, only: TLabMPI_Halos_X
!         use FDM_Derivative_MPISplit
!         integer, intent(in) :: is                       ! scalar index; if 0, then velocity
!         integer(wi), intent(in) :: nx, ny, nz
!         real(wp), intent(in) :: s(nx*ny*nz)
!         real(wp), intent(out) :: result(nx*ny*nz)
!         real(wp), intent(out), target :: tmp1(nx*ny*nz)         ! transposed field s
!         real(wp), intent(in), optional, target :: u_t(nx*ny*nz)

!         ! -------------------------------------------------------------------
!         integer np, np1, np2

!         ! ###################################################################
!         if (x%size == 1) then ! Set to zero in 2D case
!             result = 0.0_wp
!             return
!         end if

!         ! Transposition: make x-direction the last one
! #ifdef USE_ESSL
!         call DGETMO(s, nx, nx, ny*nz, tmp1, ny*nz)
! #else
!         call TLab_Transpose(s, nx, ny*nz, nx, tmp1, ny*nz)
! #endif

! ! -------------------------------------------------------------------
!         np1 = size(der1_split_x%rhs)/2
!         np2 = size(der2_split_x%rhs)/2
!         np = max(np1, np2)
!         call TLabMPI_Halos_X(tmp1, ny*nz, np, pyz_halo_m(:, 1), pyz_halo_p(:, 1))

!         if (present(u_t)) then  ! transposed velocity is passed as argument
!             p_vel => u_t
!         else
!             p_vel => tmp1
!         end if

!         call FDM_MPISplit_Solve(ny*nz, nx, der2_split_x, tmp1, &
!                                 pyz_halo_m(:, np - np2 + 1:np), pyz_halo_p, wrk3d, wrk2d)
!         call FDM_MPISplit_Solve(ny*nz, nx, der1_split_x, tmp1, &
!                                 pyz_halo_m(:, np - np1 + 1:np), pyz_halo_p, result, wrk2d)

!         if (nse_eqns == DNS_EQNS_ANELASTIC) then
!             call Anelastic(ny, nz, nx, result, wrk3d, p_vel, diffusivity(is))
!         else
!             wrk3d(1:nx*ny*nz) = wrk3d(1:nx*ny*nz)*diffusivity(is) - p_vel*result
!         end if

!         ! Put arrays back in the order in which they came in
! #ifdef USE_ESSL
!         call DGETMO(wrk3d, ny*nz, ny*nz, nx, result, nx)
! #else
!         call TLab_Transpose(wrk3d, ny*nz, nx, ny*nz, result, nx)
! #endif

!         return
!     contains
!         subroutine Anelastic(ny, nz, nx, ds1, ds2, u, diff)
!             integer(wi), intent(in) :: ny, nz, nx
!             real(wp), intent(in) :: u(ny, nz, nx), ds1(ny, nz, nx)
!             real(wp), intent(inout) :: ds2(ny, nz, nx)
!             real(wp), intent(in) :: diff

!             integer(wi) i, k

!             do i = 1, nx
!                 do k = 1, nz
!                     ds2(:, k, i) = ds2(:, k, i)*diff*ribackground(k) - u(:, k, i)*ds1(:, k, i)
!                 end do
!             end do

!             return
!         end subroutine Anelastic

!     end subroutine NSE_Burgers_PerVolume_X_MPISplit

! #endif

    !########################################################################
    !########################################################################
    subroutine NSE_Burgers_PerVolume_Y_Serial(is, nx, ny, nz, s, result, tmp1, rhou_in)
        integer, intent(in) :: is                           ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
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
        call TLab_Transpose(s, nx*ny, nz, nx*ny, tmp1, nz)
#endif

        nlines = nx*nz

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            call NSE_Burgers_1D(nlines, g(2), fdmDiffusion(2)%lu(:, :, is), &
                                s=tmp1, &
                                result=wrk3d, &
                                dsdx=result, &
                                rhou=rhou_in)
        else
            if (nse_eqns == DNS_EQNS_ANELASTIC) then
                call NSE_Burgers_1D(nlines, g(2), fdmDiffusion(2)%lu(:, :, is), &
                                    s=tmp1, &
                                    result=wrk3d, &
                                    dsdx=result, &
                                    rhou=tmp1, &
                                    rho=rho_anelastic(2)%values(:))
            else
                call NSE_Burgers_1D(nlines, g(2), fdmDiffusion(2)%lu(:, :, is), &
                                    s=tmp1, &
                                    result=wrk3d, &
                                    dsdx=result, &
                                    rhou=tmp1)
            end if
        end if

#undef rhou_out

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        call DGETMO(wrk3d, nz, nz, nx*ny, result, nx*ny)
#else
        call TLab_Transpose(wrk3d, nz, nx*ny, nz, result, nx*ny)
#endif

        return
    end subroutine NSE_Burgers_PerVolume_Y_Serial

    !########################################################################
    !########################################################################
! #ifdef USE_MPI
!     subroutine NSE_Burgers_PerVolume_Y_MPITranspose(is, nx, ny, nz, s, result, tmp1, u_t)
!         integer, intent(in) :: is                       ! scalar index; if 0, then velocity
!         integer(wi), intent(in) :: nx, ny, nz
!         real(wp), intent(in) :: s(nx*ny*nz)
!         real(wp), intent(out) :: result(nx*ny*nz)
!         real(wp), intent(out), target :: tmp1(nx*ny*nz)         ! transposed field s
!         real(wp), intent(in), target, optional :: u_t(nx*ny*nz)

!         ! -------------------------------------------------------------------
!         integer(wi) nlines

!         ! ###################################################################
!         if (y%size == 1) then ! Set to zero in 2D case
!             result = 0.0_wp
!             return
!         end if

!         ! Transposition: make y-direction the last one
! #ifdef USE_ESSL
!         call DGETMO(s, nx*ny, nx*ny, nz, wrk3d, nz)
! #else
!         call TLab_Transpose(s, nx*ny, nz, nx*ny, wrk3d, nz)
! #endif
!         call TLabMPI_Trp_ExecJ_Forward(wrk3d, tmp1, tmpi_plan_dy)
!         nlines = tmpi_plan_dy%nlines

!         if (present(u_t)) then  ! transposed velocity is passed as argument
!             p_vel => u_t
!         else
!             p_vel => tmp1
!         end if

!         call NSE_Burgers_1D(nlines, g(2), fdmDiffusion(2)%lu(:, :, is), tmp1, p_vel, result, wrk3d)

!         ! Put arrays back in the order in which they came in
!         call TLabMPI_Trp_ExecJ_Backward(result, wrk3d, tmpi_plan_dy)
! #ifdef USE_ESSL
!         call DGETMO(wrk3d, nz, nz, nx*ny, result, nx*ny)
! #else
!         call TLab_Transpose(wrk3d, nz, nx*ny, nz, result, nx*ny)
! #endif

!         return
!     end subroutine NSE_Burgers_PerVolume_Y_MPITranspose

!     !########################################################################
!     !########################################################################
!     subroutine NSE_Burgers_PerVolume_Y_MPISplit(is, nx, ny, nz, s, result, tmp1, u_t)
!         use TLabMPI_PROCS, only: TLabMPI_Halos_Y
!         use FDM_Derivative_MPISplit
!         integer, intent(in) :: is                       ! scalar index; if 0, then velocity
!         integer(wi), intent(in) :: nx, ny, nz
!         real(wp), intent(in) :: s(nx*ny*nz)
!         real(wp), intent(out) :: result(nx*ny*nz)
!         real(wp), intent(out), target :: tmp1(nx*ny*nz)         ! transposed field s
!         real(wp), intent(in), target, optional :: u_t(nx*ny*nz)

!         ! -------------------------------------------------------------------
!         integer np, np1, np2

!         ! ###################################################################
!         if (y%size == 1) then ! Set to zero in 2D case
!             result = 0.0_wp
!             return
!         end if

!         ! Transposition: make y-direction the last one
! #ifdef USE_ESSL
!         call DGETMO(s, nx*ny, nx*ny, nz, tmp1, nz)
! #else
!         call TLab_Transpose(s, nx*ny, nz, nx*ny, tmp1, nz)
! #endif

!         ! -------------------------------------------------------------------
!         np1 = size(der1_split_y%rhs)/2
!         np2 = size(der2_split_y%rhs)/2
!         np = max(np1, np2)
!         call TLabMPI_Halos_Y(tmp1, nx*nz, np, pxz_halo_m(:, 1), pxz_halo_p(:, 1))

!         if (present(u_t)) then  ! transposed velocity is passed as argument
!             p_vel => u_t
!         else
!             p_vel => tmp1
!         end if

!         call FDM_MPISplit_Solve(nx*nz, ny, der2_split_y, tmp1, &
!                                 pxz_halo_m(:, np - np2 + 1:np), pxz_halo_p, wrk3d, wrk2d)
!         call FDM_MPISplit_Solve(nx*nz, ny, der1_split_y, tmp1, &
!                                 pxz_halo_m(:, np - np1 + 1:np), pxz_halo_p, result, wrk2d)

!         if (nse_eqns == DNS_EQNS_ANELASTIC) then
!             call Anelastic(nz, nx*ny, result, wrk3d, p_vel, diffusivity(is))
!         else
!             wrk3d(1:nx*ny*nz) = wrk3d(1:nx*ny*nz)*diffusivity(is) - p_vel*result
!         end if

!         ! Put arrays back in the order in which they came in
! #ifdef USE_ESSL
!         call DGETMO(wrk3d, nz, nz, nx*ny, result, nx*ny)
! #else
!         call TLab_Transpose(wrk3d, nz, nx*ny, nz, result, nx*ny)
! #endif

!         return
!     contains
!         subroutine Anelastic(nz, nxy, ds1, ds2, u, diff)
!             integer(wi), intent(in) :: nz, nxy
!             real(wp), intent(in) :: u(nz, nxy), ds1(nz, nxy)
!             real(wp), intent(inout) :: ds2(nz, nxy)
!             real(wp), intent(in) :: diff

!             integer(wi) ij

!             do ij = 1, nxy
!                 ds2(:, ij) = ds2(:, ij)*diff*ribackground(:) - u(:, ij)*ds1(:, ij)
!             end do

!             return
!         end subroutine Anelastic

!     end subroutine NSE_Burgers_PerVolume_Y_MPISplit

! #endif

    !########################################################################
    !########################################################################
    subroutine NSE_Burgers_PerVolume_Z(is, nx, ny, nz, s, result, rhou_in, rhou_out)
        use TLab_Pointers_2D, only: pxy_wrk3d
        integer, intent(in) :: is                       ! scalar index; if 0, then velocity
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: s(nx*ny, nz)
        real(wp), intent(out) :: result(nx*ny, nz)
        real(wp), intent(in), optional :: rhou_in(nx*ny, nz)
        real(wp), intent(out), optional :: rhou_out(nx*ny, nz)

        ! -------------------------------------------------------------------
        integer(wi) k, nlines

        ! ###################################################################
        if (z%size == 1) then ! Set to zero in 2D case nx*ny
            result = 0.0_wp
            return
        end if

        nlines = nx*ny

        if (present(rhou_in)) then      ! transposed velocity (times density) is passed as argument
            call NSE_Burgers_1D(nlines, g(3), fdmDiffusion(3)%lu(:, :, is), &
                                s=s, &
                                result=result, &
                                dsdx=wrk3d, &
                                rhou=rhou_in)

        else                            ! Only used in anelastic formulation
            call NSE_Burgers_1D(nlines, g(3), fdmDiffusion(3)%lu(:, :, is), &
                                s=s, &
                                result=result, &
                                dsdx=wrk3d, &
                                rhou=s, &
                                rhou_out=rhou_out, &
                                rho=rho_anelastic(3)%values)

        end if

        if (subsidenceProps%type == TYPE_SUB_CONSTANT) then
            do k = 1, nz
                result(:, k) = result(:, k) + pxy_wrk3d(:, k)*wbackground(k)*rbackground(k)
            end do
        end if

        return
    end subroutine NSE_Burgers_PerVolume_Z

    !########################################################################
    !# Apply the non-linear operator N(u)(s) = visc* d^2/dx^2 s - rho u d/dx s
    !# along generic direction x to nlines lines of data
    !#
    !# Second derivative uses LU decomposition including diffusivity coefficient
    !########################################################################
    subroutine NSE_Burgers_1D(nlines, g, lu2d, s, result, dsdx, rhou, rhou_out, rho)
        use FDM_Derivative, only: FDM_Der1_Solve, FDM_Der2_Solve
        integer(wi), intent(in) :: nlines                           ! # of lines to be solved
        type(fdm_dt), intent(in) :: g
        real(wp), intent(in) :: lu2d(:, :)                          ! LU decomposition
        real(wp), intent(in) :: s(nlines, g%size)                   ! argument field
        real(wp), intent(out) :: result(nlines, g%size)             ! N(u) applied to s
        real(wp), intent(inout) :: dsdx(nlines, g%size)             ! dsdx
        real(wp), intent(in) :: rhou(nlines, g%size)                ! velocity times the density
        real(wp), intent(out), optional :: rhou_out(nlines, g%size) ! velocity times the density
        real(wp), intent(in), optional :: rho(g%size)

        ! -------------------------------------------------------------------
        integer(wi) k

        ! ###################################################################
        ! dsdx: 1st derivative; result: 2nd derivative including diffusivity
        call FDM_Der1_Solve(nlines, g%der1, g%der1%lu, s, dsdx, wrk2d)
        call FDM_Der2_Solve(nlines, g%der2, lu2d, s, result, dsdx, wrk2d)

        if (present(rho)) then
            do k = 1, g%size
                rhou_out(:, k) = s(:, k)*rho(k)
                result(:, k) = result(:, k) - rhou_out(:, k)*dsdx(:, k)
            end do

        else
            result(:, :) = result(:, :) - rhou(:, :)*dsdx(:, :)

        end if

        return
    end subroutine NSE_Burgers_1D

end module NSE_Burgers_PerVolume
