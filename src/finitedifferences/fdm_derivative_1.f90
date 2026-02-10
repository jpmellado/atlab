module FDM_Derivative_1order
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_DD, BCS_ND, BCS_DN, BCS_NN
    use FDM_Derivative_Base, only: matmul_halo_thomas_ice, matmul_thomas_ice
    use FDM_Derivative_Base, only: der_periodic, der_biased, der_dt
    use Thomas, only: thomas_dt, Thomas_FactorLU_InPlace
    ! use MatMul
    ! use MatMul_Halo
    use MatMul_Thomas
    use MatMul_Halo_Thomas
    use FDM_Base
    implicit none
    private

    public :: der_dt            ! Accessible by loading FDM_Derivative_* without FDM_Derivative_Base
    public :: der1_periodic
    public :: der1_biased
    public :: der1_biased_extended
    public :: FDM_Der1_ModifyWavenumbers

    ! -----------------------------------------------------------------------
    type, extends(der_periodic) :: der1_periodic
    contains
        procedure :: initialize => der1_periodic_initialize
    end type

    type, extends(der_biased) :: der1_biased
    contains
        procedure :: initialize => der1_biased_initialize
    end type

    ! -----------------------------------------------------------------------
    ! Types for biased boundary conditions

    ! for the first-order derivative, we consider different types of boundary conditions
    type :: bcs
        ! private; I need public for vpartial
        ! procedure(matmul_ice), pointer, nopass :: matmul => null()
        procedure(matmul_thomas_ice), pointer, nopass :: matmul => null()
        real(wp), pointer :: rhs(:, :) => null()
        type(thomas_dt) :: thomas
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

    type, extends(der1_biased) :: der1_biased_extended
        private
        type(bcsDN), public :: bcsDN
        type(bcsND), public :: bcsND
        type(bcsNN), public :: bcsNN
    end type

    integer(wi) nx
    integer ndl, ndr, idl, idr

contains
    ! ###################################################################
    ! ###################################################################
    subroutine der1_periodic_initialize(self, x, fdm_type)
        use FDM_ComX_Direct
        use FDM_Com1_Jacobian
        use Preconditioning
        class(der1_periodic), intent(out) :: self
        real(wp), intent(in) :: x(:)
        integer, intent(in) :: fdm_type

        ! ###################################################################
        self%type = fdm_type

        select case (fdm_type)              ! periodic implies uniform grid, direct schemes coincide with Jacobian ones
        case (FDM_COM4_DIRECT)
            self%type = FDM_COM4_JACOBIAN
        case (FDM_COM6_DIRECT)
            self%type = FDM_COM6_JACOBIAN
        end select

        select case (self%type)
        case (FDM_COM4_JACOBIAN)
            call FDM_C1N4_Jacobian(size(x), self%lhs, self%rhs, periodic=.true.)
            self%matmul => MatMul_Halo_3_antisym_ThomasL_3   ! MatMul_Halo_3_antisym + Thomas3_SolveL

        case (FDM_COM6_JACOBIAN)
            call FDM_C1N6_Jacobian(size(x), self%lhs, self%rhs, periodic=.true.)
            self%matmul => MatMul_Halo_5_antisym_ThomasL_3   ! MatMul_Halo_5_antisym + Thomas3_SolveL

        case (FDM_COM6_JACOBIAN_PENTA)
            call FDM_C1N6_Jacobian_Penta(size(x), self%lhs, self%rhs, periodic=.true.)
            self%matmul => MatMul_Halo_7_antisym_ThomasL_5   ! MatMul_Halo_7_antisym + Thomas5_SolveL

        end select

        ! Jacobian
        self%lhs = self%lhs*(x(2) - x(1))

        ! Preconditioning; the linear procedures assume 1 in the first upper diagonal.
        call Precon_Rhs(self%lhs, self%rhs, periodic=.true.)

        ! Construct LU decomposition
        call self%thomas%initialize(self%lhs)

        return
    end subroutine der1_periodic_initialize

    ! ###################################################################
    ! ###################################################################
    subroutine der1_biased_initialize(self, x, fdm_type)
        use FDM_Base, only: MultiplyByDiagonal
        use FDM_ComX_Direct
        use FDM_Com1_Jacobian
        use Preconditioning
        class(der1_biased), intent(out) :: self
        real(wp), intent(in) :: x(:)
        integer, intent(in) :: fdm_type

        real(wp), allocatable :: dx(:, :)

        ! ###################################################################
        self%type = fdm_type

        select case (self%type)
        case (FDM_COM4_JACOBIAN)
            call FDM_C1N4_Jacobian(size(x), self%lhs, self%rhs, periodic=.false.)
            self%matmul => MatMul_3_antisym_ThomasL_3   ! MatMul_3_antisym + Thomas3_SolveL

        case (FDM_COM6_JACOBIAN)
            call FDM_C1N6_Jacobian(size(x), self%lhs, self%rhs, periodic=.false.)
            self%matmul => MatMul_5_antisym_ThomasL_3   ! MatMul_5_antisym + Thomas3_SolveL

        case (FDM_COM6_JACOBIAN_PENTA)
            call FDM_C1N6_Jacobian_Penta(size(x), self%lhs, self%rhs, periodic=.false.)
            self%matmul => MatMul_7_antisym_ThomasL_5   ! MatMul_7_antisym + Thomas5_SolveL

        case (FDM_COM4_DIRECT)
            call FDM_C1N4_Direct(x, self%lhs, self%rhs)
            self%matmul => MatMul_3_ThomasL_3           ! MatMul_3 together + Thomas3_SolveL

        case (FDM_COM6_DIRECT)
            call FDM_C1N6_Direct(x, self%lhs, self%rhs)
            self%matmul => MatMul_5_ThomasL_3           ! MatMul_5 together + Thomas3_SolveL

        end select

        ! Preconditioning; Jacobian linear procedures assume 1 in the first upper diagonal.
        call Precon_Rhs(self%lhs, self%rhs, periodic=.false.)

        ! Construct LU decomposition
        call self%thomas%initialize(self%lhs)

        ! Jacobian, if needed
        select case (self%type)
        case (FDM_COM4_JACOBIAN, FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_PENTA)
            if (allocated(dx)) deallocate (dx)
            allocate (dx(1, size(x)))
            call self%compute(1, x, dx)

            call MultiplyByDiagonal(self%lhs, dx(1, :))     ! multiply by the Jacobian

            call self%thomas%initialize(self%lhs)           ! Reconstruct LU decomposition

        end select

        ! Construct LU decomposition for different types of bcs
        select type (self)
        type is (der1_biased_extended)
            call self%bcsND%initialize(self)
            call self%bcsDN%initialize(self)
            call self%bcsNN%initialize(self)
        end select

        return
    end subroutine der1_biased_initialize

    ! ###################################################################
    ! ###################################################################
    subroutine bcsND_initialize(self, ref)
        class(bcsND), intent(out) :: self
        class(der1_biased), intent(in), target :: ref

        real(wp), allocatable :: r_lhs_aux(:, :)

        ! ###################################################################
        self%matmul => ref%matmul
        self%rhs => ref%rhs

        call self%thomas%initialize(ref%lhs)

        ! adapting for bcs
        nx = size(ref%lhs, 1)
        ndl = size(ref%lhs, 2)
        idl = ndl/2 + 1
        ndr = size(ref%rhs, 2)
        idr = ndr/2 + 1

        allocate (r_lhs_aux, source=ref%lhs)
        allocate (self%rhs_b(max(idl, idr + 1), 1:ndr + 2), source=0.0_wp)
        call FDM_Der1_Neumann_Reduce(ref%lhs, ref%rhs, &
                                     r_lhs_aux, r_rhs_b=self%rhs_b)

        self%thomas%L(:, :) = r_lhs_aux(:, 1:ndl/2)
        self%thomas%U(:, :) = r_lhs_aux(:, ndl/2 + 1:ndl)
        call Thomas_FactorLU_InPlace(self%thomas%L(2:, :), self%thomas%U(2:, :))

        deallocate (r_lhs_aux)

        return
    end subroutine bcsND_initialize

    subroutine bcsND_compute(self, nlines, u, result, bcs_b)
        class(bcsND), intent(in) :: self
        integer(wi), intent(in) :: nlines
        real(wp), intent(in) :: u(nlines, size(self%rhs, 1))
        real(wp), intent(out) :: result(nlines, size(self%rhs, 1))
        real(wp), intent(inout) :: bcs_b(nlines)        ! Normal derivative as input, function value as output

        integer ic

        ! ###################################################################
        nx = size(self%rhs, 1)
        ndr = size(self%rhs, 2)
        idr = ndr/2 + 1

        result(:, 1) = bcs_b(:)

        ! Calculate RHS in A u' = B u
        call self%matmul(rhs=self%rhs, &
                         rhs_b=self%rhs_b, &
                         rhs_t=self%rhs(nx - ndr/2 + 1:nx, 1:ndr), &
                         u=u, &
                         f=result, &
                         L=self%thomas%L, &
                         bcs_b=bcs_b)

        ! Solve for u' in system of equations A u' = B u
        ! call self%thomas%ptr_solveL(self%thomas%L(2:nx, :), result(:, 2:nx))
        call self%thomas%ptr_solveU(self%thomas%U(2:nx, :), result(:, 2:nx))

        ! Calculate boundary value of u; u is not overwritten, this should be done outside if needed
        idl = size(self%thomas%U, 2)
        do ic = 1, idl - 1
            bcs_b(:) = bcs_b(:) + self%thomas%U(1, 1 + ic)*result(:, 1 + ic)
        end do
        bcs_b(:) = bcs_b(:)/self%rhs(1, idr)

        return
    end subroutine bcsND_compute

    ! ###################################################################
    ! ###################################################################
    subroutine bcsDN_initialize(self, ref)
        class(bcsDN), intent(out) :: self
        class(der1_biased), intent(in), target :: ref

        real(wp), allocatable :: r_lhs_aux(:, :)

        ! ###################################################################
        self%matmul => ref%matmul
        self%rhs => ref%rhs

        call self%thomas%initialize(ref%lhs)

        ! adapting for bcs
        nx = size(ref%lhs, 1)
        ndl = size(ref%lhs, 2)
        idl = ndl/2 + 1
        ndr = size(ref%rhs, 2)
        idr = ndr/2 + 1

        allocate (r_lhs_aux, source=ref%lhs)
        allocate (self%rhs_t(max(idl, idr + 1), 1:ndr + 2), source=0.0_wp)
        call FDM_Der1_Neumann_Reduce(ref%lhs, ref%rhs, &
                                     r_lhs_aux, r_rhs_t=self%rhs_t)

        self%thomas%L(:, :) = r_lhs_aux(:, 1:ndl/2)
        self%thomas%U(:, :) = r_lhs_aux(:, ndl/2 + 1:ndl)
        call Thomas_FactorLU_InPlace(self%thomas%L(:nx - 1, :), self%thomas%U(:nx - 1, :))

        deallocate (r_lhs_aux)

        return
    end subroutine bcsDN_initialize

    subroutine bcsDN_compute(self, nlines, u, result, bcs_t)
        class(bcsDN), intent(in) :: self
        integer(wi), intent(in) :: nlines
        real(wp), intent(in) :: u(nlines, size(self%rhs, 1))
        real(wp), intent(out) :: result(nlines, size(self%rhs, 1))
        real(wp), intent(inout) :: bcs_t(nlines)

        integer ic

        ! ###################################################################
        nx = size(self%rhs, 1)
        ndr = size(self%rhs, 2)
        idr = ndr/2 + 1

        result(:, nx) = bcs_t(:)

        ! Calculate RHS in A u' = B u
        call self%matmul(rhs=self%rhs, &
                         rhs_b=self%rhs(1:ndr/2, 1:ndr), &
                         rhs_t=self%rhs_t, &
                         u=u, &
                         f=result, &
                         L=self%thomas%L, &
                         bcs_t=bcs_t)

        ! Solve for u' in system of equations A u' = B u
        ! call self%thomas%ptr_solveL(self%thomas%L(:nx - 1, :), result(:, :nx - 1))
        call self%thomas%ptr_solveU(self%thomas%U(:nx - 1, :), result(:, :nx - 1))

        ! Calculate boundary value of u; u is not overwritten, this should be done outside if needed
        idl = size(self%thomas%U, 2)
        do ic = 1, idl - 1
            bcs_t(:) = bcs_t(:) + self%thomas%L(nx, idl - ic)*result(:, nx - ic)
        end do
        bcs_t(:) = bcs_t(:)/self%rhs(nx, idr)

        return
    end subroutine bcsDN_compute

    ! ###################################################################
    ! ###################################################################
    subroutine bcsNN_initialize(self, ref)
        class(bcsNN), intent(out) :: self
        class(der1_biased), intent(in), target :: ref

        real(wp), allocatable :: r_lhs_aux(:, :)

        ! ###################################################################
        self%matmul => ref%matmul
        self%rhs => ref%rhs

        call self%thomas%initialize(ref%lhs)

        ! adapting for bcs
        nx = size(ref%lhs, 1)
        ndl = size(ref%lhs, 2)
        idl = ndl/2 + 1
        ndr = size(ref%rhs, 2)
        idr = ndr/2 + 1

        allocate (r_lhs_aux, source=ref%lhs)
        allocate (self%rhs_b(max(idl, idr + 1), 1:ndr + 2), source=0.0_wp)
        allocate (self%rhs_t(max(idl, idr + 1), 1:ndr + 2), source=0.0_wp)
        call FDM_Der1_Neumann_Reduce(ref%lhs, ref%rhs, &
                                     r_lhs_aux, r_rhs_b=self%rhs_b, r_rhs_t=self%rhs_t)

        self%thomas%L(:, :) = r_lhs_aux(:, 1:ndl/2)
        self%thomas%U(:, :) = r_lhs_aux(:, ndl/2 + 1:ndl)
        call Thomas_FactorLU_InPlace(self%thomas%L(2:nx - 1, :), self%thomas%U(2:nx - 1, :))

        deallocate (r_lhs_aux)

        return
    end subroutine bcsNN_initialize

    subroutine bcsNN_compute(self, nlines, u, result, bcs_b, bcs_t)
        class(bcsNN), intent(in) :: self
        integer(wi), intent(in) :: nlines
        real(wp), intent(in) :: u(nlines, size(self%rhs, 1))
        real(wp), intent(out) :: result(nlines, size(self%rhs, 1))
        real(wp), intent(inout) :: bcs_b(nlines), bcs_t(nlines)

        integer ic

        ! ###################################################################
        nx = size(self%rhs, 1)
        ndr = size(self%rhs, 2)
        idr = ndr/2 + 1

        result(:, 1) = bcs_b(:)
        result(:, nx) = bcs_t(:)

        ! Calculate RHS in A u' = B u
        call self%matmul(rhs=self%rhs, &
                         rhs_b=self%rhs_b, &
                         rhs_t=self%rhs_t, &
                         u=u, &
                         f=result, &
                         L=self%thomas%L, &
                         bcs_b=bcs_b, &
                         bcs_t=bcs_t)

        ! Solve for u' in system of equations A u' = B u
        ! call self%thomas%ptr_solveL(self%thomas%L(2:nx - 1, :), result(:, 2:nx - 1))
        call self%thomas%ptr_solveU(self%thomas%U(2:nx - 1, :), result(:, 2:nx - 1))

        ! Calculate boundary value of u; u is not overwritten, this should be done outside if needed
        idl = size(self%thomas%U, 2)
        do ic = 1, idl - 1
            bcs_b(:) = bcs_b(:) + self%thomas%U(1, 1 + ic)*result(:, 1 + ic)
        end do
        bcs_b(:) = bcs_b(:)/self%rhs(1, idr)

        do ic = 1, idl - 1
            bcs_t(:) = bcs_t(:) + self%thomas%L(nx, idl - ic)*result(:, nx - ic)
        end do
        bcs_t(:) = bcs_t(:)/self%rhs(nx, idr)

        return
    end subroutine bcsNN_compute

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
        real(wp), allocatable :: locRhs_b(:, :), locRhs_t(:, :)

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

            if (allocated(locRhs_b)) deallocate (locRhs_b)
            allocate (locRhs_b(max(idl, idr + 1), ndr_b), source=0.0_wp)
            call FDM_Bcs_Reduce(BCS_MIN, aux, lhs, rhs_b=locRhs_b)

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

            if (allocated(locRhs_t)) deallocate (locRhs_t)
            allocate (locRhs_t(max(idl, idr + 1), ndr_t), source=0.0_wp)
            call FDM_Bcs_Reduce(BCS_MAX, aux, lhs, rhs_t=locRhs_t)

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

    ! #######################################################################
    ! #######################################################################
    subroutine FDM_Der1_ModifyWavenumbers(nx, c_lhs, c_rhs, modified_wn)
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
            num = 2.0_wp*c_rhs(idr)
            do ic = 1, ndr/2
                num = num + 2.0_wp*c_rhs(idr + ic)*sin(real(ic, wp)*wn(i))
            end do
            den = c_lhs(idl)
            do ic = 1, ndl/2
                den = den + 2.0_wp*c_lhs(idl + ic)*cos(real(ic, wp)*wn(i))
            end do
            modified_wn(i) = num/den

        end do

#undef wn

        return
    end subroutine FDM_Der1_ModifyWavenumbers

end module FDM_Derivative_1order
