module FDM_Base_X
    use TLab_Constants, only: wp, wi
    implicit none
    ! everything is public, so no private statement

    ! -----------------------------------------------------------------------
    type, abstract :: der_dt
        integer type                                ! finite-difference method
        real(wp), allocatable :: lhs(:, :)          ! A diagonals of system A u' = B u
        real(wp), allocatable :: rhs(:, :)          ! B diagonals of system A u' = B u
    contains
        procedure(initialize_ice), deferred :: initialize
        procedure(compute_ice), deferred :: compute
    end type
    abstract interface
        subroutine initialize_ice(self, x, dx, fdm_type)
            import der_dt, wp
            class(der_dt), intent(out) :: self
            real(wp), intent(in) :: x(:), dx(:)
            integer, intent(in) :: fdm_type
        end subroutine
        subroutine compute_ice(self, u, result)
            import der_dt, wp
            class(der_dt), intent(in) :: self
            real(wp), intent(in) :: u(:, :)
            real(wp), intent(out) :: result(:, :)
        end subroutine
    end interface

    type, extends(der_dt), abstract :: der_periodic
        ! procedure(matmul_halo_ice), pointer, nopass :: matmul => null()
        procedure(matmul_halo_thomas_ice), pointer, nopass :: matmul => null()
        procedure(thomas_ice), pointer, nopass :: thomasU => null()
        real(wp), allocatable :: mwn(:)                 ! modified wavenumbers
        real(wp), allocatable :: lu(:, :)               ! LU decomposition
        real(wp), allocatable :: z(:, :)                ! boundary corrections
    contains
    end type

    type, extends(der_dt), abstract :: der_biased
        ! procedure(matmul_ice), pointer, nopass :: matmul => null()
        procedure(matmul_thomas_ice), pointer, nopass :: matmul => null()
        procedure(thomas_ice), pointer, nopass :: thomasU => null()
    contains
    end type

    ! -----------------------------------------------------------------------
    abstract interface
        ! subroutine matmul_halo_ice(rhs, u, u_halo_m, u_halo_p, f)
        !     use TLab_Constants, only: wp
        !     real(wp), intent(in) :: rhs(:)              ! diagonals of B
        !     real(wp), intent(in) :: u(:, :)             ! vector u
        !     real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
        !     real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
        !     real(wp), intent(out) :: f(:, :)            ! vector f = B u
        ! end subroutine

        subroutine matmul_halo_thomas_ice(rhs, u, u_halo_m, u_halo_p, f, L)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: rhs(:)              ! diagonals of B
            real(wp), intent(in) :: u(:, :)             ! vector u
            real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
            real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
            real(wp), intent(out) :: f(:, :)            ! vector f = B u
            real(wp), intent(in) :: L(:, :)
        end subroutine

        ! subroutine matmul_ice(rhs, rhs_b, rhs_t, u, f, bcs_b, bcs_t)
        !     use TLab_Constants, only: wp
        !     real(wp), intent(in) :: rhs(:, :)
        !     real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
        !     real(wp), intent(in) :: u(:, :)
        !     real(wp), intent(out) :: f(:, :)
        !     real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        ! end subroutine

        subroutine matmul_thomas_ice(rhs, rhs_b, rhs_t, u, f, L, bcs_b, bcs_t)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: rhs(:, :)
            real(wp), intent(in) :: rhs_b(:, :), rhs_t(:, :)
            real(wp), intent(in) :: u(:, :)
            real(wp), intent(out) :: f(:, :)
            real(wp), intent(in) :: L(:, :)
            real(wp), intent(inout), optional :: bcs_b(:), bcs_t(:)
        end subroutine

        subroutine thomas_ice(A, f)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: A(:, :)
            real(wp), intent(inout) :: f(:, :)          ! RHS and solution
        end subroutine

    end interface

end module FDM_Base_X

! ###################################################################
! ###################################################################
module FDM_Derivative_1order_X
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_DD, BCS_ND, BCS_DN, BCS_NN
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
    public :: der1_periodic
    public :: der1_biased

    ! -----------------------------------------------------------------------
    ! Types for periodic boundary conditions
    type, extends(der_periodic) :: der1_periodic
    contains
        procedure :: initialize => der1_periodic_initialize
        procedure :: compute => der1_periodic_compute
    end type

    ! -----------------------------------------------------------------------
    ! Types for biased boundary conditions

    ! for the first-order derivative, we consider different types of boundary conditions
    type :: bcs
        private
        ! procedure(matmul_ice), pointer, nopass :: matmul => null()
        procedure(matmul_thomas_ice), pointer, nopass :: matmul => null()
        procedure(thomas_ice), pointer, nopass :: thomasU => null()
        real(wp), allocatable :: lu(:, :)
        real(wp), pointer :: rhs(:, :) => null()
    end type

    type, extends(bcs) :: bcsDD
    contains
        private
        procedure :: initialize => bcsDD_initialize
        procedure, public :: compute => bcsDD_compute
    end type

    type, extends(bcs) :: bcsND
        real(wp), allocatable :: rhs_b(:, :)
    contains
        private
        procedure :: initialize => bcsND_initialize
        procedure, public :: compute => bcsND_compute
    end type

    type, extends(bcs) :: bcsDN
        real(wp), allocatable :: rhs_t(:, :)
    contains
        private
        procedure :: initialize => bcsDN_initialize
        procedure, public :: compute => bcsDN_compute
    end type

    type, extends(bcs) :: bcsNN
        real(wp), allocatable :: rhs_b(:, :)
        real(wp), allocatable :: rhs_t(:, :)
    contains
        private
        procedure :: initialize => bcsNN_initialize
        procedure, public :: compute => bcsNN_compute
    end type

    type, extends(der_biased) :: der1_biased
        private
        type(bcsDD), public :: bcsDD
        type(bcsDN), public :: bcsDN
        type(bcsND), public :: bcsND
        type(bcsNN), public :: bcsNN
    contains
        procedure :: initialize => der1_biased_initialize
        procedure :: compute => der1_biased_compute
    end type

    integer(wi) nx, nlines
    integer ndl, ndr, idl, idr

contains
    ! ###################################################################
    ! ###################################################################
    subroutine der1_periodic_initialize(self, x, dx, fdm_type)
        use Preconditioning
        class(der1_periodic), intent(out) :: self
        real(wp), intent(in) :: x(:), dx(:)
        integer, intent(in) :: fdm_type

        ! ###################################################################
        print *, 'iniPer'

        ! -------------------------------------------------------------------
        self%type = fdm_type
        select case (fdm_type)
        case (FDM_COM4_DIRECT)
            self%type = FDM_COM4_JACOBIAN
        case (FDM_COM6_DIRECT)
            self%type = FDM_COM6_JACOBIAN
        end select

        call FDM_Der1_CreateSystem(x, self, periodic=.true.)

        call Precon_Rhs(self%lhs, self%rhs, periodic=.true.)

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
        case (5)
            call ThomasCirculant_5_Initialize(self%lu(:, 1:ndl/2), &
                                              self%lu(:, ndl/2 + 1:ndl), &
                                              self%z)
        end select

        ! -------------------------------------------------------------------
        ! Procedure pointers to matrix multiplication to calculate the right-hand side
        select case (self%type)
        case (FDM_COM4_JACOBIAN)
            self%matmul => MatMul_Halo_3_antisym_ThomasL_3   ! MatMul_Halo_3_antisym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM6_JACOBIAN)
            self%matmul => MatMul_Halo_5_antisym_ThomasL_3   ! MatMul_Halo_5_antisym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM6_JACOBIAN_PENTA)
            self%matmul => MatMul_Halo_7_antisym_ThomasL_5   ! MatMul_Halo_7_antisym together with self%thomasL => Thomas5_SolveL
            self%thomasU => Thomas5_SolveU

        end select

        return
    end subroutine der1_periodic_initialize

    ! ###################################################################
    ! ###################################################################
    subroutine der1_periodic_compute(self, u, result)
        use TLab_Arrays, only: wrk2d
        class(der1_periodic), intent(in) :: self
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: result(:, :)

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
        case (5)
            call ThomasCirculant_5_Reduce(self%lu(:, 1:ndl/2), &
                                          self%lu(:, ndl/2 + 1:ndl), &
                                          self%z, &
                                          result)!, wrk2d)
        end select

        return
    end subroutine der1_periodic_compute

    ! ###################################################################
    ! ###################################################################
    subroutine der1_biased_initialize(self, x, dx, fdm_type)
        use Preconditioning
        use FDM_Base, only: MultiplyByDiagonal
        class(der1_biased), intent(out) :: self
        real(wp), intent(in) :: x(:), dx(:)
        integer, intent(in) :: fdm_type

        ! ###################################################################
        self%type = fdm_type
        call FDM_Der1_CreateSystem(x, self, periodic=.false.)

        select case (self%type)
        case (FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA)
            call MultiplyByDiagonal(self%lhs, dx)    ! multiply by the Jacobian
        end select

        call Precon_Rhs(self%lhs, self%rhs, periodic=.false.)

        ! -------------------------------------------------------------------
        ! Procedure pointers to linear solvers
        select case (self%type)
        case (FDM_COM4_JACOBIAN)
            self%matmul => MatMul_3_antisym_ThomasL_3   ! MatMul_3_antisym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM6_JACOBIAN)
            self%matmul => MatMul_5_antisym_ThomasL_3   ! MatMul_5_antisym together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM6_JACOBIAN_PENTA)
            self%matmul => MatMul_7_antisym_ThomasL_5   ! MatMul_7_antisym together with self%thomasL => Thomas5_SolveL
            self%thomasU => Thomas5_SolveU

        case (FDM_COM4_DIRECT)
            self%matmul => MatMul_3_ThomasL_3           ! MatMul_3 together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        case (FDM_COM6_DIRECT)
            self%matmul => MatMul_5_ThomasL_3           ! MatMul_5 together with self%thomasL => Thomas3_SolveL
            self%thomasU => Thomas3_SolveU

        end select

        ! -------------------------------------------------------------------
        ! Construct LU decomposition
        call self%bcsDD%initialize(self)
        call self%bcsND%initialize(self)
        call self%bcsDN%initialize(self)
        call self%bcsNN%initialize(self)

        nullify (self%matmul)
        nullify (self%thomasU)

        return
    end subroutine der1_biased_initialize

    subroutine der1_biased_compute(self, u, result)
        class(der1_biased), intent(in) :: self
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: result(:, :)

        call self%bcsDD%compute(u, result)

        return
    end subroutine der1_biased_compute

    ! ###################################################################
    ! ###################################################################
    subroutine bcsDD_initialize(self, ref)
        class(bcsDD), intent(out) :: self
        class(der1_biased), intent(in), target :: ref

        ! ###################################################################
        print *, 'iniDD'
        self%matmul => ref%matmul
        self%thomasU => ref%thomasU
        self%rhs => ref%rhs

        nx = size(ref%lhs, 1)
        ndl = size(ref%lhs, 2)

        allocate (self%lu(nx, ndl))
        self%lu(:, :) = ref%lhs(:, :)

        call Thomas_FactorLU_InPlace(self%lu(:, 1:ndl/2), &
                                     self%lu(:, ndl/2 + 1:ndl))

        return
    end subroutine bcsDD_initialize

    subroutine bcsDD_compute(self, u, result)
        class(bcsDD), intent(in) :: self
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: result(:, :)

        ! ###################################################################
        print *, 'comDD'
        nx = size(self%lu, 1)
        ndl = size(self%lu, 2)
        ndr = size(self%rhs, 2)

        ! Calculate RHS in A u' = B u
        call self%matmul(rhs=self%rhs, &
                         rhs_b=self%rhs(1:ndr/2, 1:ndr), &
                         rhs_t=self%rhs(nx - ndr/2 + 1:nx, 1:ndr), &
                         u=u, &
                         f=result, &
                         L=self%lu(:, 1:ndl/2))

        ! Solve for u' in system of equations A u' = B u
        ! call self%thomasL(self%lu(:,1:ndl/2), result)
        call self%thomasU(self%lu(:, ndl/2 + 1:ndl), result)

        return
    end subroutine bcsDD_compute

    ! ###################################################################
    ! ###################################################################
    subroutine bcsND_initialize(self, ref)
        class(bcsND), intent(out) :: self
        class(der1_biased), intent(in), target :: ref

        ! ###################################################################
        print *, 'iniND'
        self%matmul => ref%matmul
        self%thomasU => ref%thomasU
        self%rhs => ref%rhs

        nx = size(ref%lhs, 1)
        ndl = size(ref%lhs, 2)
        idl = ndl/2 + 1
        ndr = size(ref%rhs, 2)
        idr = ndr/2 + 1

        allocate (self%lu(nx, ndl))
        self%lu(:, :) = ref%lhs(:, :)

        allocate (self%rhs_b(max(idl, idr + 1), 1:ndr + 2), source=0.0_wp)

        call FDM_Der1_Neumann_Reduce(ref%lhs, ref%rhs, &
                                     self%lu, r_rhs_b=self%rhs_b)

        call Thomas_FactorLU_InPlace(self%lu(2:nx, 1:ndl/2), &
                                     self%lu(2:nx, ndl/2 + 1:ndl))

        return
    end subroutine bcsND_initialize

    subroutine bcsND_compute(self, u, result)
        use TLab_Arrays, only: wrk2d
        class(bcsND), intent(in) :: self
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: result(:, :)

        ! ###################################################################
        print *, 'comND'
        nx = size(self%lu, 1)
        ndl = size(self%lu, 2)
        idl = ndl/2 + 1
        ndr = size(self%rhs, 2)
        idr = ndr/2 + 1

#define bcs_hb(i) wrk2d(i,1)

        nlines = size(result, 1)

        ! homogeneous Neumann bcs
        result(:, 1) = 0.0_wp
        bcs_hb(1:nlines) = 0.0_wp

        ! Calculate RHS in A u' = B u
        call self%matmul(rhs=self%rhs, &
                         rhs_b=self%rhs_b, &
                         rhs_t=self%rhs(nx - ndr/2 + 1:nx, 1:ndr), &
                         u=u, &
                         f=result, &
                         L=self%lu(:, 1:ndl/2), &
                         bcs_b=bcs_hb(1:nlines))

        ! Solve for u' in system of equations A u' = B u
        ! call self%thomasL(self%lu(:,1:ndl/2), result)
        call self%thomasU(self%lu(2:nx, ndl/2 + 1:ndl), result(:, 2:nx))

#undef bcs_hb

        return
    end subroutine bcsND_compute

    ! ###################################################################
    ! ###################################################################
    subroutine bcsDN_initialize(self, ref)
        class(bcsDN), intent(out) :: self
        class(der1_biased), intent(in), target :: ref

        ! ###################################################################
        print *, 'iniDN'
        self%matmul => ref%matmul
        self%thomasU => ref%thomasU
        self%rhs => ref%rhs

        nx = size(ref%lhs, 1)
        ndl = size(ref%lhs, 2)
        idl = ndl/2 + 1
        ndr = size(ref%rhs, 2)
        idr = ndr/2 + 1

        allocate (self%lu(nx, ndl))
        self%lu(:, :) = ref%lhs(:, :)

        allocate (self%rhs_t(max(idl, idr + 1), 1:ndr + 2), source=0.0_wp)

        call FDM_Der1_Neumann_Reduce(ref%lhs, ref%rhs, &
                                     self%lu, r_rhs_t=self%rhs_t)

        call Thomas_FactorLU_InPlace(self%lu(1:nx - 1, 1:ndl/2), &
                                     self%lu(1:nx - 1, ndl/2 + 1:ndl))

        return
    end subroutine bcsDN_initialize

    subroutine bcsDN_compute(self, u, result)
        use TLab_Arrays, only: wrk2d
        class(bcsDN), intent(in) :: self
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: result(:, :)

        ! ###################################################################
        print *, 'comDN'
        nx = size(self%lu, 1)
        ndl = size(self%lu, 2)
        idl = ndl/2 + 1
        ndr = size(self%rhs, 2)
        idr = ndr/2 + 1

#define bcs_ht(i) wrk2d(i,2)

        nlines = size(result, 1)

        ! homogeneous Neumann bcs
        result(:, nx) = 0.0_wp
        bcs_ht(1:nlines) = 0.0_wp

        ! Calculate RHS in A u' = B u
        call self%matmul(rhs=self%rhs, &
                         rhs_b=self%rhs(1:ndr/2, 1:ndr), &
                         rhs_t=self%rhs_t, &
                         u=u, &
                         f=result, &
                         L=self%lu(:, 1:ndl/2), &
                         bcs_t=bcs_ht(1:nlines))

        ! Solve for u' in system of equations A u' = B u
        ! call self%thomasL(self%lu(:,1:ndl/2), result)
        call self%thomasU(self%lu(1:nx - 1, ndl/2 + 1:ndl), result(:, 1:nx - 1))

#undef bcs_ht

        return
    end subroutine bcsDN_compute

    ! ###################################################################
    ! ###################################################################
    subroutine bcsNN_initialize(self, ref)
        class(bcsNN), intent(out) :: self
        class(der1_biased), intent(in), target :: ref

        ! ###################################################################
        print *, 'iniNN'
        self%matmul => ref%matmul
        self%thomasU => ref%thomasU
        self%rhs => ref%rhs

        nx = size(ref%lhs, 1)
        ndl = size(ref%lhs, 2)
        idl = ndl/2 + 1
        ndr = size(ref%rhs, 2)
        idr = ndr/2 + 1

        allocate (self%lu(nx, ndl))
        self%lu(:, :) = ref%lhs(:, :)

        allocate (self%rhs_b(max(idl, idr + 1), 1:ndr + 2), source=0.0_wp)
        allocate (self%rhs_t(max(idl, idr + 1), 1:ndr + 2), source=0.0_wp)

        call FDM_Der1_Neumann_Reduce(ref%lhs, ref%rhs, &
                                     self%lu, r_rhs_b=self%rhs_b, r_rhs_t=self%rhs_t)

        call Thomas_FactorLU_InPlace(self%lu(2:nx - 1, 1:ndl/2), &
                                     self%lu(2:nx - 1, ndl/2 + 1:ndl))

        return
    end subroutine bcsNN_initialize

    subroutine bcsNN_compute(self, u, result)
        use TLab_Arrays, only: wrk2d
        class(bcsNN), intent(in) :: self
        real(wp), intent(in) :: u(:, :)
        real(wp), intent(out) :: result(:, :)

        ! ###################################################################
        print *, 'comNN'
        nx = size(self%lu, 1)
        ndl = size(self%lu, 2)
        idl = ndl/2 + 1
        ndr = size(self%rhs, 2)
        idr = ndr/2 + 1

#define bcs_hb(i) wrk2d(i,1)
#define bcs_ht(i) wrk2d(i,2)

        nlines = size(result, 1)

        ! homogeneous Neumann bcs
        result(:, 1) = 0.0_wp
        bcs_hb(1:nlines) = 0.0_wp
        result(:, nx) = 0.0_wp
        bcs_ht(1:nlines) = 0.0_wp

        ! Calculate RHS in A u' = B u
        call self%matmul(rhs=self%rhs, &
                         rhs_b=self%rhs_b, &
                         rhs_t=self%rhs_t, &
                         u=u, &
                         f=result, &
                         L=self%lu(:, 1:ndl/2), &
                         bcs_b=bcs_hb(1:nlines), &
                         bcs_t=bcs_ht(1:nlines))

        ! Solve for u' in system of equations A u' = B u
        ! call self%thomasL(self%lu(2:nx - 1,1:ndl/2), result)
        call self%thomasU(self%lu(2:nx - 1, ndl/2 + 1:ndl), result(:, 2:nx - 1))

#undef bcs_hb
#undef bcs_ht

        return
    end subroutine bcsNN_compute

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_Der1_CreateSystem(x, g, periodic)
        use TLab_Constants, only: pi_wp
        use FDM_ComX_Direct
        use FDM_Com1_Jacobian
        real(wp), intent(in) :: x(:)                    ! node positions
        class(der_dt), intent(inout) :: g
        logical, intent(in) :: periodic

        ! -------------------------------------------------------------------
        real(wp) :: coef(5)
        integer(wi) i, nx

        ! ###################################################################
        nx = size(x)                    ! # grid points

        ! -------------------------------------------------------------------
        select case (g%type)
        case (FDM_COM4_JACOBIAN)
            call FDM_C1N4_Jacobian(nx, g%lhs, g%rhs, coef, periodic)

        case (FDM_COM6_JACOBIAN)
            call FDM_C1N6_Jacobian(nx, g%lhs, g%rhs, coef, periodic)

        case (FDM_COM6_JACOBIAN_PENTA)
            call FDM_C1N6_Jacobian_Penta(nx, g%lhs, g%rhs, coef, periodic)

        case (FDM_COM4_DIRECT)
            call FDM_C1N4_Direct(nx, x, g%lhs, g%rhs)

        case (FDM_COM6_DIRECT)
            call FDM_C1N6_Direct(nx, x, g%lhs, g%rhs)

        end select

        ! -------------------------------------------------------------------
        ! modified wavenumbers
        select type (g)
        type is (der1_periodic)
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

            g%mwn(:) = 2.0_wp*(coef(3)*sin(wn(:)) + coef(4)*sin(2.0_wp*wn(:)) + coef(5)*sin(3.0_wp*wn(:))) &
                       /(1.0_wp + 2.0_wp*coef(1)*cos(wn(:)) + 2.0_wp*coef(2)*cos(wn(:)))

#undef wn

        end select

        return
    end subroutine FDM_Der1_CreateSystem

! #######################################################################
! #######################################################################
    subroutine FDM_Der1_Neumann_Reduce(lhs, rhs, r_lhs, r_rhs_b, r_rhs_t)
        use TLab_Constants, only: BCS_MIN, BCS_MAX
        use FDM_Base, only: FDM_Bcs_Reduce
        real(wp), intent(in) :: lhs(:, :)
        real(wp), intent(in) :: rhs(:, :)
        real(wp), intent(out) :: r_lhs(:, :)                        ! new, reduced lhs
        real(wp), intent(inout), optional :: r_rhs_b(:, :), r_rhs_t(:, :)     ! new, reduced rhs, extended diagonals

        ! -------------------------------------------------------------------
        integer(wi) idl, ndl, idr, ndr, ir, nx
        integer(wi) idr_t, ndr_t, idr_b, ndr_b, nx_t
        real(wp), allocatable :: aux(:, :)
        real(wp) locRhs_b(7, 8), locRhs_t(7, 8)

        ! ###################################################################
        nx = size(lhs, 1)           ! # grid points

        ndl = size(lhs, 2)
        idl = ndl/2 + 1             ! center diagonal in lhs
        ndr = size(rhs, 2)
        idr = ndr/2 + 1             ! center diagonal in rhs

        ! ! For A_22, we need idl >= idr -1
        ! if (idl < idr - 1) then
        !     call TLab_Write_ASCII(efile, __FILE__//'. LHS array is too small.')
        !     call TLab_Stop(DNS_ERROR_UNDEVELOP)
        ! end if
        ! ! For b_21, we need idr >= idl
        ! if (idr < idl) then
        !     call TLab_Write_ASCII(efile, __FILE__//'. RHS array is too small.')
        !     call TLab_Stop(DNS_ERROR_UNDEVELOP)
        ! end if

        if (allocated(aux)) deallocate (aux)
        allocate (aux(1:nx, 1:ndr))

        ! -------------------------------------------------------------------
        r_lhs(:, 1:ndl) = lhs(:, 1:ndl)

        ! reorganize data
        if (present(r_rhs_b)) then
            ndr_b = size(r_rhs_b, 2)                ! can have a different # of diagonals than rhs
            idr_b = ndr_b/2 + 1

            aux(1:nx, 1:ndr) = rhs(1:nx, 1:ndr)     ! array changed in FDM_Bcs_Reduce

            locRhs_b = 0.0_wp
            call FDM_Bcs_Reduce(BCS_MIN, aux, lhs, &
                                rhs_b=locRhs_b(1:max(idl, idr + 1), 1:ndr_b))

            r_rhs_b(:, :) = 0.0_wp
            r_rhs_b(1:idr + 1, idr_b - ndr/2:idr_b + ndr/2) = aux(1:idr + 1, 1:ndr)
            r_rhs_b(1, idr_b) = lhs(1, idl)         ! save a_11 for nonzero bc
            do ir = 1, idr - 1                      ! save -a^R_{21} for nonzero bc
                r_rhs_b(1 + ir, idr_b - ir) = -locRhs_b(1 + ir, idr_b - ir)
            end do

            r_lhs(2:idl + 1, 1:ndl) = locRhs_b(2:idl + 1, idr_b - ndl/2:idr_b + ndl/2)
            r_lhs(1, idl) = rhs(1, idr)

            ! moving extended stencil in first element of old array to natural position
            r_rhs_b(1, idr_b + ndr/2 + 1) = r_rhs_b(1, idr_b - ndr/2)
            r_rhs_b(1, idr_b - ndr/2) = 0.0_wp

        end if

        if (present(r_rhs_t)) then
            ndr_t = size(r_rhs_t, 2)                ! can have a different # of diagonals than rhs
            idr_t = ndr_t/2 + 1
            nx_t = max(idl, idr + 1)

            aux(1:nx, 1:ndr) = rhs(1:nx, 1:ndr)     ! array changed in FDM_Bcs_Reduce

            locRhs_t = 0.0_wp
            call FDM_Bcs_Reduce(BCS_MAX, aux, lhs, &
                                rhs_t=locRhs_t(1:max(idl, idr + 1), 1:ndr_t))

            r_rhs_t(:, :) = 0.0_wp
            r_rhs_t(nx_t - idr:nx_t, idr_t - ndr/2:idr_t + ndr/2) = aux(nx - idr:nx, 1:ndr)
            r_rhs_t(nx_t, idr_t) = lhs(nx, idl)
            do ir = 1, idr - 1                      ! change sign in a^R_{21} for nonzero bc
                r_rhs_t(nx_t - ir, idr_t + ir) = -locRhs_t(nx_t - ir, idr_t + ir)
            end do

            r_lhs(nx - idl:nx - 1, 1:ndl) = locRhs_t(nx_t - idl:nx_t - 1, idr_t - ndl/2:idr_t + ndl/2)
            r_lhs(nx, idl) = rhs(nx, idr)

            ! moving extended stencil in last element of old array to natural position
            r_rhs_t(nx_t, idr_t - ndr/2 - 1) = r_rhs_t(nx_t, idr_t + ndr/2)
            r_rhs_t(nx_t, idr_t + ndr/2) = 0.0_wp

        end if

        if (allocated(aux)) deallocate (aux)

        return
    end subroutine FDM_Der1_Neumann_Reduce

end module FDM_Derivative_1order_X

! ! ###################################################################
! ! ###################################################################
! program test1
!     use TLab_Constants, only: wp, wi, pi_wp
!     use TLab_Arrays, only: wrk2d
!     use FDM_Derivative_1order_X
!     use FDM_Derivative

!     integer, parameter :: nx = 32
!     real(wp) x(nx), dx(nx), u(1, nx), du(1, nx), du_a(1, nx)

!     class(der_dt), allocatable :: derX

!     integer :: cases1(5) = [FDM_COM4_JACOBIAN, &
!                             FDM_COM6_JACOBIAN, &
!                             FDM_COM6_JACOBIAN_PENTA, &
!                             FDM_COM4_DIRECT, &
!                             FDM_COM6_DIRECT]

!     ! ###################################################################
!     x = [(real(i, wp), i=1, nx)]
!     dx = [(1.0_wp, i=1, nx)]
!     allocate (wrk2d(nx, 2))

!     ! u(1, :) = x(:)**2
!     ! du_a(1, :) = 2.0_wp*x(:)
!     u(1, :) = [(cos(2.0*pi_wp/(x(nx) - x(1))*(x(i) - x(1))), i=1, nx)]
!     du_a(1, :) = -[(sin(2.0*pi_wp/(x(nx) - x(1))*(x(i) - x(1))), i=1, nx)]*2.0*pi_wp/(x(nx) - x(1))

!     allocate (der1_biased :: derX)
!     do ic = 1, size(cases1)
!         call derX%initialize(x, dx, cases1(ic))
!         call derX%compute(u, du)
!         print *, maxval(abs(du - du_a))
!         select type (derX)
!         type is (der1_biased)
!             call derX%bcsDD%compute(u, du)
!             print *, maxval(abs(du - du_a))
!             call derX%bcsND%compute(u, du)
!             print *, maxval(abs(du - du_a))
!             call derX%bcsDN%compute(u, du)
!             print *, maxval(abs(du - du_a))
!             call derX%bcsNN%compute(u, du)
!             print *, maxval(abs(du - du_a))
!         end select
!     end do

!     if (allocated(derX)) deallocate (derX)
!     allocate (der1_periodic :: derX)
!     do ic = 1, size(cases1)
!         call derX%initialize(x(:nx - 1), dx(:nx - 1), cases1(ic))
!         call derX%compute(u(:, :nx - 1), du(:, :nx - 1))
!         print *, maxval(abs(du(:, :nx - 1) - du_a(:, :nx - 1)))

!     end do

!     stop
! end program
