module FDM_Derivative_2order_X
    use TLab_Constants, only: wp, wi
    use FDM_Base_X
    use Thomas
    use Thomas_Circulant
    ! use MatMul
    ! use MatMul_Halo
    use MatMul_Thomas
    use MatMul_Halo_Thomas
    use FDM_Derivative          ! to be removed
    implicit none
    private

    public :: der_dt            ! Made public to make it accessible by loading FDM_Derivative_X and not necessarily FDM_Base_X
    public :: der2_periodic
    public :: der2_biased

    ! -----------------------------------------------------------------------
    ! Types for periodic boundary conditions
    type, extends(der_periodic) :: der2_periodic
        private
        logical :: need_1der = .false.                  ! In nonuniform, Jacobian formulation, we need 1. order derivative for the 2. order one
        real(wp), allocatable :: rhs_d1(:, :)           ! 1. order derivative correction in nonuniform case
    contains
        procedure :: initialize => der2_periodic_initialize
        procedure :: compute => der2_periodic_compute
    end type

    ! -----------------------------------------------------------------------
    ! Types for biased boundary conditions
    type, extends(der_biased) :: der2_biased
        private
        logical :: need_1der = .false.                  ! In nonuniform, Jacobian formulation, we need 1. order derivative for the 2. order one
        real(wp), allocatable :: lu(:, :)               ! LU decomposition
        real(wp), allocatable :: rhs_d1(:, :)           ! 1. order derivative correction in nonuniform case
    contains
        procedure :: initialize => der2_biased_initialize
        procedure :: compute => der2_biased_compute
    end type

contains
    ! ###################################################################
    ! ###################################################################
    subroutine der2_periodic_initialize(self, x, dx, fdm_type)
        use Preconditioning
        class(der2_periodic), intent(out) :: self
        real(wp), intent(in) :: x(:), dx(:)
        integer, intent(in) :: fdm_type

        integer nx, ndl

        ! ###################################################################
        print *, 'iniPer'

        ! -------------------------------------------------------------------
        self%type = fdm_type
        select case (fdm_type)
        case (FDM_COM4_DIRECT)
            self%type = FDM_COM4_JACOBIAN
        case (FDM_COM6_DIRECT)
            self%type = FDM_COM6_JACOBIAN
        case (FDM_COM6_DIRECT_HYPER)
            self%type = FDM_COM6_JACOBIAN_HYPER
        end select

        call FDM_Der2_CreateSystem(x, dx, self, periodic=.true., uniform=.true.)

        call NormalizeByDiagonal(self%rhs, &
                                 1, &                           ! use 1. upper diagonal in rhs
                                 self%lhs, &
                                 switchAtBoundary=.false.)

        ! ! For code readability later in the code
        ! g%nb_diag = [size(g%lhs, 2), size(g%rhs, 2)]

        nx = size(self%lhs, 1)
        ndl = size(self%lhs, 2)

        ! -------------------------------------------------------------------
        ! Construct LU decomposition
        allocate (self%lu(nx, ndl))
        self%lu(:, :) = self%lhs(:, :)

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
        case (FDM_COM4_JACOBIAN)
            self%matmul => MatMul_Halo_3_sym_ThomasL_3      ! MatMul_Halo_3_sym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM6_JACOBIAN)
            self%matmul => MatMul_Halo_5_sym_ThomasL_3      ! MatMul_Halo_5_sym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        end select

        return
    end subroutine der2_periodic_initialize

    subroutine der2_periodic_compute(self, u, result)
        use TLab_Arrays, only: wrk2d
        class(der2_periodic), intent(in) :: self
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: result(:, :)

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
    subroutine der2_biased_initialize(self, x, dx, fdm_type)
        use Preconditioning
        class(der2_biased), intent(out) :: self
        real(wp), intent(in) :: x(:), dx(:)
        integer, intent(in) :: fdm_type

        ! ###################################################################
        self%type = fdm_type
        select case (fdm_type)
        case (FDM_COM6_DIRECT_HYPER)
            self%type = FDM_COM6_JACOBIAN_HYPER
        end select

        call FDM_Der2_CreateSystem(x, dx, self, periodic=.false., uniform=.true.)   ! uniform to be fixed

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

        ! ! For code readability later in the code
        ! g%nb_diag = [size(g%lhs, 2), size(g%rhs, 2)]

        ! -------------------------------------------------------------------
        ! Construct LU decomposition
        allocate (self%lu, mold=self%lhs)
        self%lu(:, :) = self%lhs(:, :)

        ! -------------------------------------------------------------------
        ! Procedure pointers to linear solvers
        select case (self%type)
        case (FDM_COM4_JACOBIAN)
            self%matmul => MatMul_5_sym_ThomasL_3       ! MatMul_5_sym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM6_JACOBIAN)
            self%matmul => MatMul_5_sym_ThomasL_3       ! MatMul_5_sym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM6_JACOBIAN_HYPER)
            self%matmul => MatMul_7_sym_ThomasL_3       ! MatMul_7_sym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM4_DIRECT)
            self%matmul => MatMul_3_ThomasL_3           ! MatMul_3 together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM6_DIRECT)
            self%matmul => MatMul_5_ThomasL_3           ! MatMul_5 together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        end select

        return
    end subroutine der2_biased_initialize

    subroutine der2_biased_compute(self, u, result)
        class(der2_biased), intent(in) :: self
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: result(:, :)

        result = u

        return
    end subroutine der2_biased_compute

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Der2_CreateSystem(x, dx, g, periodic, uniform)
        use TLab_Constants, only: pi_wp
        use FDM_Base, only: MultiplyByDiagonal
        use FDM_ComX_Direct
        use FDM_Com2_Jacobian
        real(wp), intent(in) :: x(:)                    ! node positions
        real(wp), intent(in) :: dx(size(x, 1), 2)       ! Jacobians
        class(der_dt), intent(inout) :: g               ! fdm plan for 2. order derivative
        logical, intent(in) :: periodic, uniform

        ! -------------------------------------------------------------------
        real(wp) :: coef(5)
        integer(wi) i, nx

        ! ###################################################################
        nx = size(x)                    ! # grid points

        ! -------------------------------------------------------------------
        select type (g)
        type is (der2_periodic)
            select case (g%type)
            case (FDM_COM4_JACOBIAN)
                call FDM_C2N4_Jacobian(nx, g%lhs, g%rhs, coef, periodic)

            case (FDM_COM6_JACOBIAN)
                call FDM_C2N6_Jacobian(nx, g%lhs, g%rhs, coef, periodic)

            case (FDM_COM6_JACOBIAN_HYPER)
                call FDM_C2N6_Hyper_Jacobian(nx, g%lhs, g%rhs, coef, periodic)

            case (FDM_COM4_DIRECT)
                call FDM_C2N4_Direct(nx, x, g%lhs, g%rhs)

            case (FDM_COM6_DIRECT)
                call FDM_C2N6_Direct(nx, x, g%lhs, g%rhs)

            end select

            select case (g%type)
            case (FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_HYPER)
                if (.not. uniform) then
                    g%need_1der = .true.
                    if (allocated(g%rhs_d1)) deallocate (g%rhs_d1)  ! Contribution from 1. order derivative in nonuniform grids
                    allocate (g%rhs_d1, mold=g%lhs)
                    g%rhs_d1 = -g%lhs
                    call MultiplyByDiagonal(g%rhs_d1, dx(:, 2))
                end if
                call MultiplyByDiagonal(g%lhs, dx(:, 1))            ! multiply by the Jacobians
                call MultiplyByDiagonal(g%lhs, dx(:, 1))
            end select

        end select

        ! -------------------------------------------------------------------
        ! modified wavenumbers
        select type (g)
        type is (der2_periodic)
            if (allocated(g%mwn)) deallocate (g%mwn)
            allocate (g%mwn(nx))
#define wn(i) g%mwn(i)

            do i = 1, nx        ! wavenumbers, the independent variable to construct the modified ones
                if (i <= nx/2 + 1) then
                    wn(i) = 2.0_wp*pi_wp*real(i - 1, wp)/real(nx, wp)
                else
                    wn(i) = 2.0_wp*pi_wp*real(i - 1 - nx, wp)/real(nx, wp)
                end if
            end do

            g%mwn(:) = 2.0_wp*(coef(3)*(1.0_wp - cos(wn(:))) + coef(4)*(1.0_wp - cos(2.0_wp*wn(:))) + coef(5)*(1.0_wp - cos(3.0_wp*wn(:)))) &
                       /(1.0_wp + 2.0_wp*coef(1)*cos(wn(:)) + 2.0_wp*coef(2)*cos(2.0_wp*wn(:)))

#undef wn

        end select

        return
    end subroutine FDM_Der2_CreateSystem

end module FDM_Derivative_2order_X

! ###################################################################
! ###################################################################
program test2
    use TLab_Constants, only: wp, wi, pi_wp
    use TLab_Arrays, only: wrk2d
    use FDM_Derivative_2order_X
    use FDM_Derivative

    integer, parameter :: nx = 32
    real(wp) x(nx), dx(nx), u(1, nx), du(1, nx), du_a(1, nx)

    class(der_dt), allocatable :: derX

    integer :: cases1(5) = [FDM_COM4_JACOBIAN, &
                            FDM_COM6_JACOBIAN, &
                            FDM_COM6_JACOBIAN_HYPER, &
                            FDM_COM4_DIRECT, &
                            FDM_COM6_DIRECT]

    ! ###################################################################
    x = [(real(i, wp), i=1, nx)]
    dx = [(1.0_wp, i=1, nx)]
    allocate (wrk2d(nx, 2))

    ! u(1, :) = x(:)**2
    ! du_a(1, :) = 2.0_wp*x(:)
    u(1, :) = [(cos(2.0*pi_wp/(x(nx) - x(1))*(x(i) - x(1))), i=1, nx)]
    du_a(1, :) = -[(sin(2.0*pi_wp/(x(nx) - x(1))*(x(i) - x(1))), i=1, nx)]*2.0*pi_wp/(x(nx) - x(1))

    allocate (der2_biased :: derX)
    do ic = 1, size(cases1)
        call derX%initialize(x, dx, cases1(ic))
        call derX%compute(u, du)
        print *, maxval(abs(du - du_a))
    end do

    if (allocated(derX)) deallocate (derX)
    allocate (der2_periodic :: derX)
    do ic = 1, size(cases1)
        call derX%initialize(x(:nx - 1), dx(:nx - 1), cases1(ic))
        call derX%compute(u(:, :nx - 1), du(:, :nx - 1))
        print *, maxval(abs(du(:, :nx - 1) - du_a(:, :nx - 1)))

    end do

    stop
end program
