module FDM_Derivative_Burgers
    use TLab_Constants, only: wp, wi
    use FDM_Base, only: FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_HYPER
    use FDM_Derivative_Base
    use FDM_Derivative_1order, only: der1_periodic
    use FDM_Derivative_2order, only: der2_periodic
    implicit none
    private

    public :: der1_periodic, der2_periodic

    type, public :: der_burgers
        procedure(matmul_halo_thomas_combined_ice), pointer, nopass :: matmul => null()
        procedure(matmul_halo_thomas_ice), pointer, nopass :: matmul1 => null()
        procedure(matmul_halo_thomas_ice), pointer, nopass :: matmul2 => null()
        procedure(thomas_ice), pointer, nopass :: thomasU1 => null()
        procedure(thomas_ice), pointer, nopass :: thomasU2 => null()
        real(wp), pointer :: lu1(:, :) => null()
        real(wp), pointer :: z1(:, :) => null()
        real(wp), pointer :: rhs1(:, :) => null()
        real(wp), pointer :: lu2(:, :) => null()
        real(wp), pointer :: z2(:, :) => null()
        real(wp), pointer :: rhs2(:, :) => null()

    contains
        procedure :: initialize => der_burgers_periodic_initialize
        procedure :: compute => der_burgers_periodic_compute

    end type

    ! -----------------------------------------------------------------------
    abstract interface
        subroutine matmul_halo_thomas_combined_ice(rhs1, rhs2, u, u_halo_m, u_halo_p, f, L1, g, L2)
            import wp
            real(wp), intent(in) :: rhs1(:)             ! diagonals of B1
            real(wp), intent(in) :: rhs2(:)             ! diagonals of B2
            real(wp), intent(in) :: u(:, :)             ! vector u
            real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
            real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
            real(wp), intent(out) :: f(:, :)            ! vector f = B1 u
            real(wp), intent(in) :: L1(:, :)
            real(wp), intent(out) :: g(:, :)            ! vector g = B2 u
            real(wp), intent(in) :: L2(:, :)
        end subroutine
    end interface

contains
    ! ###################################################################
    ! ###################################################################
    subroutine der_burgers_periodic_initialize(self, fdm_der1, fdm_der2)
        use Matmul_Halo_Thomas
        class(der_burgers), intent(out) :: self
        type(der1_periodic), intent(in), target :: fdm_der1
        type(der2_periodic), intent(in), target :: fdm_der2

        ! ###################################################################
        self%lu1 => fdm_der1%lu
        self%z1 => fdm_der1%z
        self%rhs1 => fdm_der1%rhs
        self%matmul1 => fdm_der1%matmul
        self%thomasU1 => fdm_der1%thomasU

        self%lu2 => fdm_der2%lu
        self%z2 => fdm_der2%z
        self%rhs2 => fdm_der2%rhs
        self%matmul2 => fdm_der2%matmul
        self%thomasU2 => fdm_der2%thomasU

        if (fdm_der1%type == FDM_COM6_JACOBIAN .and. fdm_der2%type == FDM_COM6_JACOBIAN_HYPER) then
            self%matmul => MatMul_Halo_5_antisym_7_sym_ThomasL_3
        else
            print *, 'undeveloped'
            stop
        end if

        return
    end subroutine der_burgers_periodic_initialize

    ! ###################################################################
    ! ###################################################################
    subroutine der_burgers_periodic_compute(self, nlines, u, du1, du2)
        use TLab_Arrays, only: wrk2d
        use Thomas_Circulant
        class(der_burgers), intent(in) :: self
        integer(wi), intent(in) :: nlines
        real(wp), intent(in) :: u(nlines, size(self%lu1, 1))
        real(wp), intent(out) :: du1(nlines, size(self%lu1, 1))
        real(wp), intent(out) :: du2(nlines, size(self%lu1, 1))

        integer nx, ndl1, ndr1, ndl2, ndr2

        ! ###################################################################
        nx = size(self%lu1, 1)
        ndl1 = size(self%lu1, 2)
        ndr1 = size(self%rhs1, 2)
        ndl2 = size(self%lu2, 2)
        ndr2 = size(self%rhs2, 2)

        ! Calculate RHS in system of equations A u' = B u
        call self%matmul(rhs1=self%rhs1(1, :), &
                         rhs2=self%rhs2(1, :), &
                         u=u, &
                         u_halo_m=u(:, nx - ndr2/2 + 1:nx), &
                         u_halo_p=u(:, 1:ndr2/2), &
                         f=du1, &
                         L1=self%lu1(:, 1:ndl1/2), &
                         g=du2, &
                         L2=self%lu2(:, 1:ndl2/2))

        ! call self%matmul1(rhs=self%rhs1(1, :), &
        !                   u=u, &
        !                   u_halo_m=u(:, nx - ndr1/2 + 1:nx), &
        !                   u_halo_p=u(:, 1:ndr1/2), &
        !                   f=du1, &
        !                   L=self%lu1(:, 1:ndl1/2))

        ! call self%matmul2(rhs=self%rhs2(1, :), &
        !                   u=u, &
        !                   u_halo_m=u(:, nx - ndr2/2 + 1:nx), &
        !                   u_halo_p=u(:, 1:ndr2/2), &
        !                   f=du2, &
        !                   L=self%lu2(:, 1:ndl2/2))

        ! Solve for u' in system of equations A u' = B u1
        call self%thomasU1(self%lu1(:, ndl1/2 + 1:ndl1), du1)
        call self%thomasU2(self%lu2(:, ndl2/2 + 1:ndl2), du2)
        call ThomasCirculant_3_Reduce(self%lu1(:, 1:ndl1/2), &
                                      self%lu1(:, ndl1/2 + 1:ndl1), &
                                      self%z1(1, :), &
                                      du1, wrk2d(:, 1))
        call ThomasCirculant_3_Reduce(self%lu2(:, 1:ndl2/2), &
                                      self%lu2(:, ndl2/2 + 1:ndl2), &
                                      self%z2(1, :), &
                                      du2, wrk2d(:, 1))

        return
    end subroutine der_burgers_periodic_compute

end module FDM_Derivative_Burgers
