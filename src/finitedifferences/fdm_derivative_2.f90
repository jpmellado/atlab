module FDM_Derivative_2order
    use TLab_Constants, only: wp, wi
    use FDM_Derivative_Base!, only: matmul_halo_thomas_ice, matmul_thomas_ice, thomas_ice
    use Thomas
    use Thomas_Circulant
    ! use MatMul
    ! use MatMul_Halo
    use MatMul_Thomas
    use MatMul_Halo_Thomas
    use FDM_Base
    implicit none
    private

    public :: der_extended_dt
    ! public :: der_dt            ! Made public to make it accessible by loading FDM_Derivative_X and not necessarily FDM_Derivative_Base
    public :: der2_extended_periodic
    public :: der2_extended_biased
    public :: FDM_Der2_ModifyWavenumbers

    ! -----------------------------------------------------------------------
    ! Types for periodic boundary conditions
    type, extends(der_periodic) :: der2_periodic
    contains
        procedure :: initialize => der2_periodic_initialize
        procedure :: compute => der2_periodic_compute
    end type

    ! ! Types for biased boundary conditions
    ! type, extends(der_biased) :: der2_biased
    !     private
    !     real(wp), allocatable :: lu(:, :)               ! LU decomposition
    ! contains
    !     procedure :: initialize => der2_biased_initialize
    !     procedure :: compute => der2_biased_compute
    ! end type

    ! -----------------------------------------------------------------------
    ! The idea is to ude FDM_Base for both 1. and 2. order derivatives, but the current implementation
    ! of some schemes requires u' in the calculation of u'', and we need to redefine most of it
    ! We define wrappers so that we can use most of the general fdm structure

    type, abstract :: der_extended_dt
        integer type                                ! finite-difference method
        real(wp), allocatable :: lhs(:, :)          ! A diagonals of system A u' = B u
        real(wp), allocatable :: rhs(:, :)          ! B diagonals of system A u' = B u
    contains
        procedure(initialize_extended_ice), deferred :: initialize
        procedure(compute_extended_ice), deferred :: compute
    end type
    abstract interface
        subroutine initialize_extended_ice(self, x, fdm_type, fdm_der1, uniform)
            use FDM_Derivative_1order, only: der_dt
            import der_extended_dt, wp
            class(der_extended_dt), intent(out) :: self
            real(wp), intent(in) :: x(:)
            integer, intent(in) :: fdm_type
            class(der_dt), intent(in), optional :: fdm_der1
            logical, intent(in), optional :: uniform
        end subroutine
        subroutine compute_extended_ice(self, nlines, u, result, du)
            import der_extended_dt, wp, wi
            class(der_extended_dt), intent(in) :: self
            integer(wi), intent(in) :: nlines
            real(wp), intent(in) :: u(nlines, size(self%lhs, 1))
            real(wp), intent(out) :: result(nlines, size(self%lhs, 1))
            real(wp), intent(in), optional :: du(nlines, size(self%lhs, 1))
        end subroutine
    end interface

    type, extends(der_extended_dt), abstract :: der_extended_periodic
        ! procedure(matmul_halo_ice), pointer, nopass :: matmul => null()
        procedure(matmul_halo_thomas_ice), pointer, nopass :: matmul => null()
        procedure(thomas_ice), pointer, nopass :: thomasU => null()
        real(wp), allocatable :: lu(:, :)               ! LU decomposition
        real(wp), allocatable :: z(:, :)                ! boundary corrections
    contains
    end type

    type, extends(der_extended_dt), abstract :: der_extended_biased
        ! procedure(matmul_ice), pointer, nopass :: matmul => null()
        procedure(matmul_thomas_ice), pointer, nopass :: matmul => null()
        procedure(thomas_ice), pointer, nopass :: thomasU => null()
    contains
    end type

    ! -----------------------------------------------------------------------
    ! Types for periodic boundary conditions
    type, extends(der_extended_periodic) :: der2_extended_periodic
        ! private
        type(der2_periodic) :: der2
    contains
        procedure :: initialize => der2_extended_periodic_initialize
        procedure :: compute => der2_extended_periodic_compute
    end type

    ! -----------------------------------------------------------------------
    ! Types for biased boundary conditions

    ! I need to consider a case with Jacobian, nonuniform
    type :: bcs
        private
        procedure(thomas_ice), pointer, nopass :: thomasU => null()
        real(wp), allocatable :: lu(:, :)
        real(wp), pointer :: rhs(:, :) => null()
    end type

    type, extends(bcs) :: bcsDD
        ! procedure(matmul_ice), pointer, nopass :: matmul => null()
        procedure(matmul_thomas_ice), pointer, nopass :: matmul => null()
    contains
        private
        procedure :: initialize => bcsDD_initialize
        procedure, public :: compute => bcsDD_compute
    end type

    type, extends(bcs) :: bcsDD_jacobian
        ! procedure(matmul_add_ice), pointer, nopass :: matmul => null()
        procedure(matmul_add_thomas_ice), pointer, nopass :: matmul => null()
        real(wp), allocatable :: rhs_d1(:, :)
    contains
        private
        procedure :: initialize => bcsDD_jacobian_initialize
        procedure, public :: compute => bcsDD_jacobian_compute
    end type

    type, extends(der_extended_biased) :: der2_extended_biased
        private
        logical :: need_1der = .false.  ! In nonuniform, Jacobian formulation, we need 1. order derivative for the 2. order one
        type(bcsDD), public :: bcsDD
        type(bcsDD_jacobian), public :: bcsDD_jacobian
    contains
        procedure :: initialize => der2_extended_biased_initialize
        procedure :: compute => der2_extended_biased_compute
    end type

    ! -----------------------------------------------------------------------
    abstract interface
        ! subroutine matmul_add_ice(rhs, rhs_b, rhs_t, u, f, rhs_add, u_add, bcs_b, bcs_t)
        !     use TLab_Constants, only: wp
        !     real(wp), intent(in) :: rhs(:, :)
        !     real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        !     real(wp), intent(in) :: u(:, :)
        !     real(wp), intent(out) :: f(:, :)
        !     real(wp), intent(in) :: rhs_add(:, :)
        !     real(wp), intent(in) :: u_add(:, :)
        !     real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        ! end subroutine

        subroutine matmul_add_thomas_ice(rhs, rhs_b, rhs_t, u, f, rhs_add, u_add, L, bcs_b, bcs_t)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: rhs(:, :)
            real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
            real(wp), intent(in) :: u(:, :)
            real(wp), intent(out) :: f(:, :)
            real(wp), intent(in) :: rhs_add(:, :)
            real(wp), intent(in) :: u_add(:, :)
            real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
            real(wp), intent(in) :: L(:, :)
        end subroutine

    end interface

contains
    ! ###################################################################
    ! ###################################################################
    ! Wrappers
    subroutine der2_extended_periodic_initialize(self, x, fdm_type, fdm_der1, uniform)
        use FDM_Derivative_1order, only: der_dt
        class(der2_extended_periodic), intent(out) :: self
        real(wp), intent(in) :: x(:)
        integer, intent(in) :: fdm_type
        class(der_dt), optional, intent(in) :: fdm_der1
        logical, intent(in), optional :: uniform

        call self%der2%initialize(x, fdm_type)

        ! I need it in elliptic operators
        allocate (self%lhs, source=self%der2%lhs)
        allocate (self%rhs, source=self%der2%rhs)

        return
    end subroutine der2_extended_periodic_initialize

    subroutine der2_extended_periodic_compute(self, nlines, u, result, du)
        class(der2_extended_periodic), intent(in) :: self
        integer(wi), intent(in) :: nlines
        real(wp), intent(in) :: u(nlines, size(self%lhs, 1))
        real(wp), intent(out) :: result(nlines, size(self%lhs, 1))
        real(wp), intent(in), optional :: du(nlines, size(self%lhs, 1))

        call self%der2%compute(nlines, u, result)

        return
    end subroutine der2_extended_periodic_compute

    ! ###################################################################
    ! ###################################################################
    subroutine der2_periodic_initialize(self, x, fdm_type)
        use FDM_ComX_Direct
        use FDM_Com2_Jacobian
        use Preconditioning
        class(der2_periodic), intent(out) :: self
        real(wp), intent(in) :: x(:)
        integer, intent(in) :: fdm_type

        integer nx, ndl

        ! ###################################################################
        self%type = fdm_type

        select case (fdm_type)              ! periodic implies uniform grid, direct schemes coincide with Jacobian ones
        case (FDM_COM4_DIRECT)
            self%type = FDM_COM4_JACOBIAN
        case (FDM_COM6_DIRECT)
            self%type = FDM_COM6_JACOBIAN
        case (FDM_COM6_DIRECT_HYPER)
            self%type = FDM_COM6_JACOBIAN_HYPER
        end select

        select case (self%type)
        case (FDM_COM4_JACOBIAN)
            call FDM_C2N4_Jacobian(size(x), self%lhs, self%rhs, periodic=.true.)
            self%matmul => MatMul_Halo_5_sym_ThomasL_3      ! MatMul_Halo_3_sym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM6_JACOBIAN)
            call FDM_C2N6_Jacobian(size(x), self%lhs, self%rhs, periodic=.true.)
            self%matmul => MatMul_Halo_5_sym_ThomasL_3      ! MatMul_Halo_3_sym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM6_JACOBIAN_HYPER)
            call FDM_C2N6_Hyper_Jacobian(size(x), self%lhs, self%rhs, periodic=.true.)
            self%matmul => MatMul_Halo_7_sym_ThomasL_3      ! MatMul_Halo_7_sym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        end select

        ! Jacobian
        self%lhs = self%lhs*(x(2) - x(1))*(x(2) - x(1))

        ! Preconditioning; Jacobian linear procedures assume 1 in the first upper diagonal.
        call NormalizeByDiagonal(self%rhs, &
                                 1, &                       ! use 1. upper diagonal in rhs to normalize system
                                 self%lhs, &
                                 switchAtBoundary=.false.)

        ! Construct LU decomposition
        nx = size(self%lhs, 1)
        ndl = size(self%lhs, 2)

        allocate (self%lu, source=self%lhs)

        allocate (self%z(ndl/2, nx))

        select case (ndl)
        case (3)
            call ThomasCirculant_3_Initialize(self%lu(:, 1:ndl/2), &
                                              self%lu(:, ndl/2 + 1:ndl), &
                                              self%z)
        end select

        return
    end subroutine der2_periodic_initialize

    ! ###################################################################
    ! ###################################################################
    subroutine der2_periodic_compute(self, nlines, u, result)
        use TLab_Arrays, only: wrk2d
        class(der2_periodic), intent(in) :: self
        integer(wi), intent(in) :: nlines
        real(wp), intent(in) :: u(nlines, size(self%lhs, 1))
        real(wp), intent(out) :: result(nlines, size(self%lhs, 1))

        integer nx, ndl, ndr

        ! ###################################################################
        nx = size(self%lhs, 1)
        ndl = size(self%lhs, 2)
        ndr = size(self%rhs, 2)

        ! Calculate RHS in system of equations A u' = B u
        ! call self%matmul(rhs=self%rhs(1, 1:ndr), &
        !                    u=u, &
        !                    u_halo_m=u(:, nx - ndr/2 + 1:nx), &
        !                    u_halo_p=u(:, 1:ndr/2), &
        !                    f=result)
        call self%matmul(rhs=self%rhs(1, 1:ndr), &
                         u=u, &
                         u_halo_m=u(:, nx - ndr/2 + 1:nx), &
                         u_halo_p=u(:, 1:ndr/2), &
                         f=result, &
                         L=self%lu(:, 1:ndl/2))

        ! Solve for u' in system of equations A u' = B u
        ! call self%thomasL(self%lu(:, 1:ndl/2), result)
        call self%thomasU(self%lu(:, ndl/2 + 1:ndl), result)
        select case (ndl)
        case (3)
            call ThomasCirculant_3_Reduce(self%lu(:, 1:ndl/2), &
                                          self%lu(:, ndl/2 + 1:ndl), &
                                          self%z(1, :), &
                                          result, wrk2d(:, 1))
        end select

        return
    end subroutine der2_periodic_compute

    ! ###################################################################
    ! ###################################################################
    subroutine der2_extended_biased_initialize(self, x, fdm_type, fdm_der1, uniform)
        use FDM_Base, only: MultiplyByDiagonal
        use FDM_Derivative_1order, only: der_dt
        use FDM_ComX_Direct
        use FDM_Com2_Jacobian
        use Preconditioning
        class(der2_extended_biased), intent(out) :: self
        real(wp), intent(in) :: x(:)
        integer, intent(in) :: fdm_type
        class(der_dt), intent(in), optional :: fdm_der1
        logical, intent(in), optional :: uniform

        real(wp), allocatable :: x_aux(:, :), dx(:, :)
        integer i

        ! ###################################################################
        self%type = fdm_type

        select case (fdm_type)
        case (FDM_COM6_DIRECT_HYPER)                    ! Not yet implemented; fall back to Jacobian version
            self%type = FDM_COM6_JACOBIAN_HYPER
        end select

        select case (self%type)
        case (FDM_COM4_JACOBIAN)
            call FDM_C2N4_Jacobian(size(x), self%lhs, self%rhs, periodic=.false.)
            self%matmul => MatMul_5_sym_ThomasL_3       ! MatMul_5_sym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM6_JACOBIAN)
            call FDM_C2N6_Jacobian(size(x), self%lhs, self%rhs, periodic=.false.)
            self%matmul => MatMul_5_sym_ThomasL_3       ! MatMul_5_sym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM6_JACOBIAN_HYPER)
            call FDM_C2N6_Hyper_Jacobian(size(x), self%lhs, self%rhs, periodic=.false.)
            self%matmul => MatMul_7_sym_ThomasL_3       ! MatMul_7_sym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM4_DIRECT)
            call FDM_C2N4_Direct(x, self%lhs, self%rhs)
            self%matmul => MatMul_5_ThomasL_3           ! MatMul_5 together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM6_DIRECT)
            call FDM_C2N6_Direct(x, self%lhs, self%rhs)
            self%matmul => MatMul_5_ThomasL_3           ! MatMul_5 together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        end select

        ! Preconditioning; Jacobian linear procedures assume 1 in the first upper diagonal.
        call NormalizeByDiagonal(self%rhs, &
                                 1, &                   ! use 1. upper diagonal in rhs
                                 self%lhs, &
                                 switchAtBoundary=.true.)

        ! Construct LU decomposition
        call self%bcsDD%initialize(self)

        ! Jacobian, if needed
        select case (self%type)
        case (FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_HYPER)
            if (allocated(x_aux)) deallocate (x_aux)
            allocate (x_aux(1, size(x)))
            if (allocated(dx)) deallocate (dx)
            allocate (dx(1, size(x)))

            ! Calculating derivative dxds; calculate dsdx and invert it
            x_aux(1, :) = [(real(i - 1, wp), i=1, size(x))]
            call fdm_der1%compute(1, x_aux, dx)
            dx(1, :) = 1.0_wp/dx(1, :)

            call MultiplyByDiagonal(self%lhs, dx(1, :))     ! multiply by the Jacobians
            call MultiplyByDiagonal(self%lhs, dx(1, :))

            call self%bcsDD%initialize(self)                ! Reconstruct LU decomposition

            ! Contribution from 1. order derivative in nonuniform grids
            if (.not. uniform) then
                self%need_1der = .true.
                call self%bcsDD_jacobian%initialize(self, x)
            end if

        end select

        return
    end subroutine der2_extended_biased_initialize

    subroutine der2_extended_biased_compute(self, nlines, u, result, du)
        class(der2_extended_biased), intent(in) :: self
        integer(wi), intent(in) :: nlines
        real(wp), intent(in) :: u(nlines, size(self%lhs, 1))
        real(wp), intent(out) :: result(nlines, size(self%lhs, 1))
        real(wp), intent(in), optional :: du(nlines, size(self%lhs, 1))

        ! ###################################################################
        if (self%need_1der) then           ! add Jacobian correction A_2 dx2 du
            call self%bcsDD_jacobian%compute(nlines, u, result, du)
        else
            call self%bcsDD%compute(nlines, u, result)
        end if

        return
    end subroutine der2_extended_biased_compute

    ! ###################################################################
    ! ###################################################################
    subroutine bcsDD_initialize(self, ref)
        class(bcsDD), intent(out) :: self
        class(der2_extended_biased), intent(in), target :: ref

        integer ndl

        ! ###################################################################
        self%matmul => ref%matmul
        self%thomasU => ref%thomasU
        self%rhs => ref%rhs

        allocate (self%lu, source=ref%lhs)

        ndl = size(ref%lhs, 2)

        call Thomas_FactorLU_InPlace(self%lu(:, 1:ndl/2), &
                                     self%lu(:, ndl/2 + 1:ndl))

        return
    end subroutine bcsDD_initialize

    subroutine bcsDD_compute(self, nlines, u, result)
        class(bcsDD), intent(in) :: self
        integer(wi), intent(in) :: nlines
        real(wp), intent(in) :: u(nlines, size(self%lu, 1))
        real(wp), intent(out) :: result(nlines, size(self%lu, 1))

        integer nx, ndl, ndr

        ! ###################################################################
        nx = size(self%lu, 1)
        ndl = size(self%lu, 2)
        ndr = size(self%rhs, 2)

        ! call self%matmul(rhs=self%rhs, &
        !                    rhs_b=self%rhs(1:ndr/2, 1:ndr), &
        !                    rhs_t=self%rhs(nx - ndr/2 + 1:nx, 1:ndr), &
        !                    u=u, &
        !                    f=result)
        call self%matmul(rhs=self%rhs, &
                         rhs_b=self%rhs(1:ndr/2, 1:ndr), &
                         rhs_t=self%rhs(nx - ndr/2 + 1:nx, 1:ndr), &
                         u=u, &
                         f=result, &
                         L=self%lu(:, 1:ndl/2))

        ! Solve for u' in system of equations A u' = B u
        ! call self%thomasL(lu(:, 1:ndl/2), result)
        call self%thomasU(self%lu(:, ndl/2 + 1:ndl), result)

        return
    end subroutine bcsDD_compute

    ! ###################################################################
    ! ###################################################################
    subroutine bcsDD_jacobian_initialize(self, ref, x)
        use Preconditioning
        use FDM_Base, only: MultiplyByDiagonal
        class(bcsDD_jacobian), intent(out) :: self
        class(der2_extended_biased), intent(in), target :: ref
        real(wp), intent(in) :: x(:)

        real(wp), allocatable :: dx(:, :)

        ! ###################################################################
        select case (ref%type)
        case (FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN)
            self%matmul => MatMul_5_sym_add_3_ThomasL_3     ! MatMul_5_sym together with self%thomasL => Thomas3_SolveL
        case (FDM_COM6_JACOBIAN_HYPER)
            self%matmul => MatMul_7_sym_add_3_ThomasL_3     ! MatMul_7_sym together with self%thomasL => Thomas3_SolveL
        end select
        self%thomasU => ref%thomasU
        self%rhs => ref%rhs

        allocate (self%lu, source=ref%bcsDD%lu)

        ! Contribution from 1. order derivative in nonuniform grids
        allocate (self%rhs_d1, mold=ref%lhs)
        self%rhs_d1 = -ref%lhs

        if (allocated(dx)) deallocate (dx)                  ! Calculate dx2ds2
        allocate (dx(1, size(x)))
        call ref%bcsDD%compute(1, x, dx)

        call MultiplyByDiagonal(self%rhs_d1, dx(1, :))      ! multiply by the Jacobian

        return
    end subroutine bcsDD_jacobian_initialize

    subroutine bcsDD_jacobian_compute(self, nlines, u, result, du)
        class(bcsDD_jacobian), intent(in) :: self
        integer(wi), intent(in) :: nlines
        real(wp), intent(in) :: u(nlines, size(self%lu, 1))
        real(wp), intent(out) :: result(nlines, size(self%lu, 1))
        real(wp), intent(in) :: du(nlines, size(self%lu, 1))

        integer nx, ndl, ndr

        ! ###################################################################
        nx = size(self%lu, 1)
        ndl = size(self%lu, 2)
        ndr = size(self%rhs, 2)

        ! call self%matmul(rhs=self%rhs, &
        !                        rhs_b=self%rhs(1:ndr/2, 1:ndr), &
        !                        rhs_t=self%rhs(nx - ndr/2 + 1:nx, 1:ndr), &
        !                        u=u, &
        !                        rhs_add=self%rhs(:, ip + 1:ip + 3), &
        !                        u_add=du, &
        !                        f=result)
        call self%matmul(rhs=self%rhs, &
                         rhs_b=self%rhs(1:ndr/2, 1:ndr), &
                         rhs_t=self%rhs(nx - ndr/2 + 1:nx, 1:ndr), &
                         u=u, &
                         rhs_add=self%rhs_d1, &
                         u_add=du, &
                         f=result, &
                         L=self%lu(:, 1:ndl/2))

        ! Solve for u' in system of equations A u' = B u
        ! call self%thomasL(lu(:, 1:ndl/2), result)
        call self%thomasU(self%lu(:, ndl/2 + 1:ndl), result)

        return
    end subroutine bcsDD_jacobian_compute

    ! #######################################################################
    ! #######################################################################
    subroutine FDM_Der2_ModifyWavenumbers(nx, c_lhs, c_rhs, modified_wn)
        use TLab_Constants, only: pi_wp
        integer, intent(in) :: nx
        real(wp), intent(in) :: c_lhs(:), c_rhs(:)              ! coefficients of A u' = B u
        real(wp), allocatable, intent(out) :: modified_wn(:)

        integer i, ic, ndl, idl, ndr, idr
        real(wp) num, den

        ! #######################################################################
        allocate (modified_wn(nx), source=0.0_wp)

        ndl = size(c_lhs, 1)
        idl = ndl/2 + 1
        ndr = size(c_rhs, 1)
        idr = ndr/2 + 1

#define wn(i) modified_wn(i)

        do i = 1, nx
            ! wavenumbers, the independent variable to construct the modified ones
            if (i <= nx/2 + 1) then
                wn(i) = 2.0_wp*pi_wp*real(i - 1, wp)/real(nx, wp)
            else
                wn(i) = 2.0_wp*pi_wp*real(i - 1 - nx, wp)/real(nx, wp)
            end if

            ! compute modified wavenumbers
            num = 0.0_wp
            do ic = 1, ndr/2
                num = num + 2.0_wp*c_rhs(idr + ic)*(1.0_wp - cos(real(ic, wp)*wn(i)))
            end do
            den = c_lhs(idl)
            do ic = 1, ndl/2
                den = den + 2.0_wp*c_lhs(idl + ic)*cos(real(ic, wp)*wn(i))
            end do
            modified_wn(i) = num/den

        end do

#undef wn

        return
    end subroutine FDM_Der2_ModifyWavenumbers

end module FDM_Derivative_2order
