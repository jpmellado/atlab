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
        type(der1_periodic), pointer :: der1 => null()
        type(der2_periodic), pointer :: der2 => null()
    contains
        procedure :: initialize => der_burgers_periodic_initialize
        procedure :: compute => der_burgers_periodic_compute
    end type

contains
    ! ###################################################################
    ! ###################################################################
    subroutine der_burgers_periodic_initialize(self, fdm_der1, fdm_der2)
        use Matmul_Halo_Thomas
        class(der_burgers), intent(out) :: self
        class(der1_periodic), intent(in), target :: fdm_der1
        class(der2_periodic), intent(in), target :: fdm_der2

        ! ###################################################################
        self%der1 => fdm_der1
        self%der2 => fdm_der2

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
        real(wp), intent(in) :: u(nlines, size(self%der1%lu, 1))
        real(wp), intent(out) :: du1(nlines, size(self%der1%lu, 1))
        real(wp), intent(out) :: du2(nlines, size(self%der1%lu, 1))

        integer nx, ndl1, ndr1, ndl2, ndr2

        ! ###################################################################
        nx = size(self%der1%lu, 1)
        ndl1 = size(self%der1%lu, 2)
        ndr1 = size(self%der1%rhs, 2)
        ndl2 = size(self%der2%lu, 2)
        ndr2 = size(self%der2%rhs, 2)

        ! Calculate RHS in system of equations A u' = B u
        call self%matmul(rhs1=self%der1%rhs(1, :), &
                         rhs2=self%der2%rhs(1, :), &
                         u=u, &
                         u_halo_m=u(:, nx - ndr2/2 + 1:nx), &
                         u_halo_p=u(:, 1:ndr2/2), &
                         f=du1, &
                         L1=self%der1%lu(:, 1:ndl1/2), &
                         g=du2, &
                         L2=self%der2%lu(:, 1:ndl2/2))

        ! call self%der1%matmul(rhs=self%der1%rhs(1, :), &
        !                       u=u, &
        !                       u_halo_m=u(:, nx - ndr1/2 + 1:nx), &
        !                       u_halo_p=u(:, 1:ndr1/2), &
        !                       f=du1, &
        !                       L=self%der1%lu(:, 1:ndl1/2))

        ! call self%der2%matmul(rhs=self%der2%rhs(1, :), &
        !                       u=u, &
        !                       u_halo_m=u(:, nx - ndr2/2 + 1:nx), &
        !                       u_halo_p=u(:, 1:ndr2/2), &
        !                       f=du2, &
        !                       L=self%der2%lu(:, 1:ndl2/2))

        ! Solve for u' in system of equations A u' = B u1
        call self%der1%thomasU(self%der1%lu(:, ndl1/2 + 1:ndl1), du1)
        call self%der2%thomasU(self%der2%lu(:, ndl2/2 + 1:ndl2), du2)
        call ThomasCirculant_3_Reduce(self%der1%lu(:, 1:ndl1/2), &
                                      self%der1%lu(:, ndl1/2 + 1:ndl1), &
                                      self%der1%z(1, :), &
                                      du1, wrk2d(:, 1))
        call ThomasCirculant_3_Reduce(self%der2%lu(:, 1:ndl2/2), &
                                      self%der2%lu(:, ndl2/2 + 1:ndl2), &
                                      self%der2%z(1, :), &
                                      du2, wrk2d(:, 1))

        return
    end subroutine der_burgers_periodic_compute

end module FDM_Derivative_Burgers
