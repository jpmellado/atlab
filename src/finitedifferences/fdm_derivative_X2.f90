module FDM_Derivative_2order_X
    use TLab_Constants, only: wp, wi
    use FDM_Base_X, only: matmul_halo_thomas_ice, matmul_thomas_ice, thomas_ice
    use Thomas
    use Thomas_Circulant
    ! use MatMul
    ! use MatMul_Halo
    use MatMul_Thomas
    use MatMul_Halo_Thomas
    use FDM_Derivative          ! to be removed
    implicit none
    private

    public :: der2_dt
    ! public :: der_dt            ! Made public to make it accessible by loading FDM_Derivative_X and not necessarily FDM_Base_X
    public :: der2_periodic
    public :: der2_biased

    !

    ! -----------------------------------------------------------------------
    ! ! Types for periodic boundary conditions
    ! type, extends(der_periodic) :: der2_periodic
    ! contains
    !     procedure :: initialize => der2_periodic_initialize
    !     procedure :: compute => der2_periodic_compute
    ! end type

    ! ! Types for biased boundary conditions
    ! type, extends(der_biased) :: der2_biased
    !     private
    !     real(wp), allocatable :: lu(:, :)               ! LU decomposition
    ! contains
    !     procedure :: initialize => der2_biased_initialize
    !     procedure :: compute => der2_biased_compute
    ! end type

    ! The idea is to ude FDM_Base for both 1. and 2. order derivatives, but the current implementation
    ! of some schemes requires u' in the calculation of u'', and we need to redefine most of it

    type, abstract :: der2_dt
        integer type                                ! finite-difference method
        real(wp), allocatable :: lhs(:, :)          ! A diagonals of system A u' = B u
        real(wp), allocatable :: rhs(:, :)          ! B diagonals of system A u' = B u
    contains
        procedure(initialize_ice), deferred :: initialize
        procedure(compute_ice), deferred :: compute
    end type
    abstract interface
        subroutine initialize_ice(self, x, dx, fdm_type, uniform)
            import der2_dt, wp
            class(der2_dt), intent(out) :: self
            real(wp), intent(in) :: x(:), dx(:, :)
            integer, intent(in) :: fdm_type
            logical :: uniform
        end subroutine
        subroutine compute_ice(self, u, result, du)
            import der2_dt, wp
            class(der2_dt), intent(in) :: self
            real(wp), intent(in) :: u(:, :)
            real(wp), intent(out) :: result(:, :)
            real(wp), intent(in), optional :: du(:, :)
        end subroutine
    end interface

    type, extends(der2_dt), abstract :: der_periodic
        ! procedure(matmul_halo_ice), pointer, nopass :: matmul => null()
        procedure(matmul_halo_thomas_ice), pointer, nopass :: matmul => null()
        procedure(thomas_ice), pointer, nopass :: thomasU => null()
        real(wp), allocatable :: mwn(:)                 ! modified wavenumbers
        real(wp), allocatable :: lu(:, :)               ! LU decomposition
        real(wp), allocatable :: z(:, :)                ! boundary corrections
    contains
    end type

    type, extends(der2_dt), abstract :: der_biased
        ! procedure(matmul_ice), pointer, nopass :: matmul => null()
        procedure(matmul_thomas_ice), pointer, nopass :: matmul => null()
        procedure(thomas_ice), pointer, nopass :: thomasU => null()
    contains
    end type

    ! Types for periodic boundary conditions
    type, extends(der_periodic) :: der2_periodic
    contains
        procedure :: initialize => der2_periodic_initialize
        procedure :: compute => der2_periodic_compute
    end type

    ! Types for biased boundary conditions
    type, extends(der_biased) :: der2_biased
        private
        ! procedure(matmul_add_ice), pointer, nopass :: matmul_add => null()
        procedure(matmul_add_thomas_ice), pointer, nopass :: matmul_add => null()
        logical :: need_1der = .false.                  ! In nonuniform, Jacobian formulation, we need 1. order derivative for the 2. order one
        real(wp), allocatable :: rhs_d1(:, :)           ! 1. order derivative correction in nonuniform case
        real(wp), allocatable :: lu(:, :)               ! LU decomposition
    contains
        procedure :: initialize => der2_biased_initialize
        procedure :: compute => der2_biased_compute
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
    subroutine der2_periodic_initialize(self, x, dx, fdm_type, uniform)
        use Preconditioning
        class(der2_periodic), intent(out) :: self
        real(wp), intent(in) :: x(:), dx(:, :)
        integer, intent(in) :: fdm_type
        logical :: uniform

        integer nx, ndl

        ! ###################################################################
        print *, 'iniPer'

        self%type = fdm_type

        select case (fdm_type)              ! periodic implies uniform grid, direct schemes coincide with Jacobian ones
        case (FDM_COM4_DIRECT)
            self%type = FDM_COM4_JACOBIAN
        case (FDM_COM6_DIRECT)
            self%type = FDM_COM6_JACOBIAN
        case (FDM_COM6_DIRECT_HYPER)
            self%type = FDM_COM6_JACOBIAN_HYPER
        end select

        call FDM_Der2_CreateSystem(x, self, periodic=.true.)

        call NormalizeByDiagonal(self%rhs, &
                                 1, &                           ! use 1. upper diagonal in rhs to normalize system
                                 self%lhs, &
                                 switchAtBoundary=.false.)

        ! -------------------------------------------------------------------
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

        ! -------------------------------------------------------------------
        ! Procedure pointers to matrix multiplication to calculate the right-hand side
        select case (self%type)
        case (FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN)
            self%matmul => MatMul_Halo_5_sym_ThomasL_3      ! MatMul_Halo_3_sym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM6_JACOBIAN_HYPER)
            self%matmul => MatMul_Halo_7_sym_ThomasL_3      ! MatMul_Halo_7_sym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        end select

        return
    end subroutine der2_periodic_initialize

    subroutine der2_periodic_compute(self, u, result, du)
        use TLab_Arrays, only: wrk2d
        class(der2_periodic), intent(in) :: self
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: result(:, :)
        real(wp), intent(in), optional :: du(:, :)

        integer nx, ndl, ndr

        ! ###################################################################
        print *, 'comPer'
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
    subroutine der2_biased_initialize(self, x, dx, fdm_type, uniform)
        use Preconditioning
        use FDM_Base, only: MultiplyByDiagonal
        class(der2_biased), intent(out) :: self
        real(wp), intent(in) :: x(:), dx(:, :)
        integer, intent(in) :: fdm_type
        logical :: uniform

        integer ndl

        ! ###################################################################
        print *, 'iniDD'

        self%type = fdm_type

        select case (fdm_type)
        case (FDM_COM6_DIRECT_HYPER)                ! Not yet implemented; fall back to Jacobian version
            self%type = FDM_COM6_JACOBIAN_HYPER
        end select

        call FDM_Der2_CreateSystem(x, self, periodic=.false.)

        ! Jacobian, if needed
        select case (self%type)
        case (FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_HYPER)
            if (.not. uniform) then
                self%need_1der = .true.
                if (allocated(self%rhs_d1)) deallocate (self%rhs_d1)  ! Contribution from 1. order derivative in nonuniform grids
                allocate (self%rhs_d1, mold=self%lhs)
                self%rhs_d1 = -self%lhs
                call MultiplyByDiagonal(self%rhs_d1, dx(:, 2))
            end if
            call MultiplyByDiagonal(self%lhs, dx(:, 1))            ! multiply by the Jacobians
            call MultiplyByDiagonal(self%lhs, dx(:, 1))
        end select

        ! Preconditioning
        if (self%need_1der) then
            call NormalizeByDiagonal(self%rhs, &
                                     1, &                           ! use 1. upper diagonal in rhs
                                     self%lhs, &
                                     self%rhs_d1, &
                                     switchAtBoundary=.true.)
        else
            call NormalizeByDiagonal(self%rhs, &
                                     1, &                           ! use 1. upper diagonal in rhs
                                     self%lhs, &
                                     switchAtBoundary=.true.)
        end if

        ! -------------------------------------------------------------------
        ! Construct LU decomposition
        allocate (self%lu, source=self%lhs)

        ndl = size(self%lhs, 2)
        call Thomas_FactorLU_InPlace(self%lu(:, 1:ndl/2), &
                                     self%lu(:, ndl/2 + 1:ndl))

        ! -------------------------------------------------------------------
        ! Procedure pointers to linear solvers
        select case (self%type)
        case (FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN)
            if (self%need_1der) then                                ! add Jacobian correction A_2 dx2 du
                self%matmul_add => MatMul_5_sym_add_3_ThomasL_3     ! MatMul_5_sym together with self%thomasL => Thomas3_SolveL
            else
                self%matmul => MatMul_5_sym_ThomasL_3               ! MatMul_5_sym together with self%thomasL => Thomas3_SolveL
            end if
            self%thomasU => Thomas3_SolveU

        case (FDM_COM6_JACOBIAN_HYPER)
            if (self%need_1der) then                                ! add Jacobian correction A_2 dx2 du
                self%matmul_add => MatMul_7_sym_add_3_ThomasL_3     ! MatMul_7_sym together with self%thomasL => Thomas3_SolveL
            else
                self%matmul => MatMul_7_sym_ThomasL_3               ! MatMul_7_sym together with self%thomasL => Thomas3_SolveL
            end if
            self%thomasU => Thomas3_SolveU

        case (FDM_COM4_DIRECT, FDM_COM6_DIRECT)
            self%matmul => MatMul_5_ThomasL_3           ! MatMul_5 together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        end select

        return
    end subroutine der2_biased_initialize

    subroutine der2_biased_compute(self, u, result, du)
        class(der2_biased), intent(in) :: self
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: result(:, :)
        real(wp), intent(in), optional :: du(:, :)

        integer nx, ndl, ndr

        ! ###################################################################
        print *, 'comDD'
        nx = size(self%lu, 1)
        ndl = size(self%lu, 2)
        ndr = size(self%rhs, 2)

        if (self%need_1der) then           ! add Jacobian correction A_2 dx2 du
            ! call self%matmul(rhs=self%rhs, &
            !                        rhs_b=self%rhs(1:ndr/2, 1:ndr), &
            !                        rhs_t=self%rhs(nx - ndr/2 + 1:nx, 1:ndr), &
            !                        u=u, &
            !                        rhs_add=self%rhs(:, ip + 1:ip + 3), &
            !                        u_add=du, &
            !                        f=result)
            call self%matmul_add(rhs=self%rhs, &
                                 rhs_b=self%rhs(1:ndr/2, 1:ndr), &
                                 rhs_t=self%rhs(nx - ndr/2 + 1:nx, 1:ndr), &
                                 u=u, &
                                 rhs_add=self%rhs_d1, &
                                 u_add=du, &
                                 f=result, &
                                 L=self%lu(:, 1:ndl/2))
        else
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
        end if

        ! Solve for u' in system of equations A u' = B u
        ! call self%thomasL(lu(:, 1:ndl/2), result)
        call self%thomasU(self%lu(:, ndl/2 + 1:ndl), result)

        return
    end subroutine der2_biased_compute

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Der2_CreateSystem(x, g, periodic)
        use TLab_Constants, only: pi_wp
        use FDM_ComX_Direct
        use FDM_Com2_Jacobian
        real(wp), intent(in) :: x(:)            ! node positions
        class(der2_dt), intent(inout) :: g      ! fdm plan for 2. order derivative
        logical, intent(in) :: periodic

        ! -------------------------------------------------------------------
        real(wp) :: coef(5)
        integer(wi) nx

        ! ###################################################################
        nx = size(x)                            ! # grid points

        ! -------------------------------------------------------------------
        select case (g%type)
        case (FDM_COM4_JACOBIAN)
            call FDM_C2N4_Jacobian(nx, g%lhs, g%rhs, coef, periodic)

        case (FDM_COM6_JACOBIAN)
            call FDM_C2N6_Jacobian(nx, g%lhs, g%rhs, coef, periodic)

        case (FDM_COM6_JACOBIAN_HYPER)
            call FDM_C2N6_Hyper_Jacobian(nx, g%lhs, g%rhs, coef, periodic)

        case (FDM_COM4_DIRECT)
            call FDM_C2N4_Direct(x, g%lhs, g%rhs)

        case (FDM_COM6_DIRECT)
            call FDM_C2N6_Direct(x, g%lhs, g%rhs)

        end select

        select type (g)
        type is (der2_periodic)
            call FDM_Der2_ModifyWavenumbers(nx, coef, g%mwn)
        end select

        return
    end subroutine FDM_Der2_CreateSystem

    subroutine FDM_Der2_ModifyWavenumbers(nx, coef, modified_wn)
        use TLab_Constants, only: pi_wp
        integer, intent(in) :: nx
        real(wp), intent(in) :: coef(:)
        real(wp), allocatable, intent(out) :: modified_wn(:)

        integer i

        allocate (modified_wn(nx), source=0.0_wp)

#define wn(i) modified_wn(i)

        do i = 1, nx        ! wavenumbers, the independent variable to construct the modified ones
            if (i <= nx/2 + 1) then
                wn(i) = 2.0_wp*pi_wp*real(i - 1, wp)/real(nx, wp)
            else
                wn(i) = 2.0_wp*pi_wp*real(i - 1 - nx, wp)/real(nx, wp)
            end if
        end do

        modified_wn(:) = 2.0_wp*(coef(3)*(1.0_wp - cos(wn(:))) + coef(4)*(1.0_wp - cos(2.0_wp*wn(:))) + coef(5)*(1.0_wp - cos(3.0_wp*wn(:)))) &
                         /(1.0_wp + 2.0_wp*coef(1)*cos(wn(:)) + 2.0_wp*coef(2)*cos(2.0_wp*wn(:)))

#undef wn

        return
    end subroutine FDM_Der2_ModifyWavenumbers

end module FDM_Derivative_2order_X

! ! ###################################################################
! ! ###################################################################
! program test2
!     use TLab_Constants, only: wp, wi, pi_wp
!     use TLab_Arrays, only: wrk2d
!     use FDM_Derivative_2order_X
!     use FDM_Derivative

!     integer, parameter :: nx = 32
!     real(wp) x(nx), dx(nx, 2), u(1, nx), du(1, nx), du_a(1, nx)

!     class(der2_dt), allocatable :: derX

!     integer :: cases1(5) = [FDM_COM4_JACOBIAN, &
!                             FDM_COM6_JACOBIAN, &
!                             FDM_COM6_JACOBIAN_HYPER, &
!                             FDM_COM4_DIRECT, &
!                             FDM_COM6_DIRECT]

!     ! ###################################################################
!     x = [(real(i, wp), i=1, nx)]
!     dx(:, 1) = [(1.0_wp, i=1, nx)]
!     dx(:, 2) = [(0.0_wp, i=1, nx)]
!     allocate (wrk2d(nx, 2))

!     ! u(1, :) = x(:)**2
!     ! du_a(1, :) = 2.0_wp*x(:)
!     u(1, :) = [(cos(2.0*pi_wp/(x(nx) - x(1))*(x(i) - x(1))), i=1, nx)]
!     ! du_a(1, :) = -[(sin(2.0*pi_wp/(x(nx) - x(1))*(x(i) - x(1))), i=1, nx)]*2.0*pi_wp/(x(nx) - x(1))
!     du_a(1, :) = -u(1, :)*(2.0*pi_wp/(x(nx) - x(1)))**2

!     allocate (der2_biased :: derX)
!     do ic = 1, size(cases1)
!         call derX%initialize(x, dx, cases1(ic), uniform=.true.)
!         call derX%compute(u, du)
!         print *, maxval(abs(du - du_a))
!     end do

!     if (allocated(derX)) deallocate (derX)
!     allocate (der2_periodic :: derX)
!     do ic = 1, size(cases1)
!         call derX%initialize(x(:nx - 1), dx(:nx - 1, :), cases1(ic), uniform=.true.)
!         call derX%compute(u(:, :nx - 1), du(:, :nx - 1))
!         print *, maxval(abs(du(:, :nx - 1) - du_a(:, :nx - 1)))
!     end do

!     stop
! end program
