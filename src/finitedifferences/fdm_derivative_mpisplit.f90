#include "tlab_error.h"

! Split algorithm assumes periodic case over uniform grid
module FDM_Derivative_MPISplit
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use FDM_Derivative_Base, only: der_periodic
    use FDM_Derivative_Base, only: matmul_halo_thomas_ice, matmul_halo_thomas_combined_ice
    use FDM_Derivative_1order, only: der1_periodic
    use FDM_Derivative_2order, only: der2_periodic
    use Thomas_Parallel, only: thomas_parallel_dt
    implicit none
    private

    public :: der_periodic_mpisplit
    public :: der_burgers_mpisplit

    ! -----------------------------------------------------------------------
    type :: der_mpisplit
        integer type                                ! finite-difference method
        real(wp), pointer :: rhs(:, :) => null()
    contains
    end type

    type, extends(der_mpisplit) :: der_periodic_mpisplit
        ! procedure(matmul_halo_ice), pointer, nopass :: matmul => null()
        procedure(matmul_halo_thomas_ice), pointer, nopass :: matmul => null()
        type(thomas_parallel_dt) thomas3
    contains
        procedure :: initialize => der_periodic_initialize
        procedure :: compute => der_periodic_compute
    end type

    type :: der_burgers_mpisplit
        procedure(matmul_halo_thomas_combined_ice), pointer, nopass :: matmul => null()
        type(der_periodic_mpisplit), pointer :: der1 => null()
        type(der_periodic_mpisplit), pointer :: der2 => null()
    contains
        procedure :: initialize => der_burgers_periodic_initialize
        procedure :: compute => der_burgers_periodic_compute
    end type

contains
    ! ###################################################################
    ! ###################################################################
    subroutine der_periodic_initialize(self, ref, mpiAxis)
        use TLabMPI_VARS, only: mpi_axis_dt
        class(der_periodic_mpisplit), intent(out) :: self
        class(der_periodic), intent(in), target :: ref
        type(mpi_axis_dt) mpiAxis

        integer nx, k, np

        ! ###################################################################
        if (size(ref%lhs, 2) /= 3) then
            call TLab_Write_ASCII(efile, __FILE__//'. Matrix splitting only for tridiagonal.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        self%type = ref%type

        self%matmul => ref%matmul
        self%rhs => ref%rhs

        nx = size(ref%lhs, 1)

        np = mpiAxis%num_processors     ! for clarity below
        call self%thomas3%initialize(ref%lhs, &
                                     [(k, k=nx/np, nx, nx/np)], &
                                     block_id=mpiAxis%rank + 1, &
                                     circulant=.true.)
        self%thomas3%mpi = mpiAxis

        return
    end subroutine

    ! ###################################################################
    ! ###################################################################
    subroutine der_periodic_compute(self, nlines, u, u_halo_m, u_halo_p, result)
        use TLab_Arrays, only: wrk2d
        class(der_periodic_mpisplit), intent(in) :: self
        integer(wi), intent(in) :: nlines
        real(wp), intent(in) :: u(nlines, size(self%thomas3%L, 1))
        real(wp), intent(in) :: u_halo_m(:, :)
        real(wp), intent(in) :: u_halo_p(:, :)
        real(wp), intent(out) :: result(nlines, size(self%thomas3%L, 1))

        ! ###################################################################
        call self%matmul(rhs=self%rhs(1, :), &
                         u=u, &
                         u_halo_m=u_halo_m, &
                         u_halo_p=u_halo_p, &
                         f=result, &
                         L=self%thomas3%L)

        ! call self%thomas3%SolveL(result)
        call self%thomas3%SolveU(result)
        call self%thomas3%reduce(result, wrk2d(:, 1), wrk2d(:, 2))

        return
    end subroutine

    ! ###################################################################
    ! ###################################################################
    subroutine der_burgers_periodic_initialize(self, fdm_der1, fdm_der2)
        use FDM_Base, only: FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_HYPER
        use Matmul_Halo_Thomas, only: MatMul_Halo_5_antisym_7_sym_ThomasL_3
        class(der_burgers_mpisplit), intent(out) :: self
        class(der_periodic_mpisplit), intent(in), target :: fdm_der1
        class(der_periodic_mpisplit), intent(in), target :: fdm_der2

        ! ###################################################################
        self%der1 => fdm_der1
        self%der2 => fdm_der2

        if (fdm_der1%type == FDM_COM6_JACOBIAN .and. &
            fdm_der2%type == FDM_COM6_JACOBIAN_HYPER) then
            self%matmul => MatMul_Halo_5_antisym_7_sym_ThomasL_3
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Burgers splitting only for tridiagonal.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        return
    end subroutine

    ! ###################################################################
    ! ###################################################################
    subroutine der_burgers_periodic_compute(self, nlines, u, u_halo_m, u_halo_p, du1, du2)
        use TLab_Arrays, only: wrk2d
        class(der_burgers_mpisplit), intent(in) :: self
        integer(wi), intent(in) :: nlines
        real(wp), intent(in) :: u(nlines, size(self%der1%thomas3%L, 1))
        real(wp), intent(in) :: u_halo_m(:, :)
        real(wp), intent(in) :: u_halo_p(:, :)
        real(wp), intent(out) :: du1(nlines, size(self%der1%thomas3%L, 1))
        real(wp), intent(out) :: du2(nlines, size(self%der1%thomas3%L, 1))

        integer ndr1, ndr2, np

        ! ###################################################################
        ndr1 = size(self%der1%rhs, 2)
        ndr2 = size(self%der2%rhs, 2)
        np = size(u_halo_m, 2)
        call self%matmul(rhs1=self%der1%rhs(1, :), &
                         rhs2=self%der2%rhs(1, :), &
                         u=u, &
                         u_halo_m=u_halo_m(:, np - ndr2/2 + 1:np), &
                         u_halo_p=u_halo_p, &
                         f=du1, &
                         L1=self%der1%thomas3%L, &
                         g=du2, &
                         L2=self%der2%thomas3%L)

        ! call self%der1%matmul(rhs=self%der1%rhs(1, :), &
        !                       u=u, &
        !                       u_halo_m=u_halo_m(:, np - ndr1/2 + 1:np), &
        !                       u_halo_p=u_halo_p, &
        !                       f=du1, &
        !                       L=self%der1%thomas3%lhs(:, 1:1))

        ! call self%der2%matmul(rhs=self%der2%rhs(1, :), &
        !                       u=u, &
        !                       u_halo_m=u_halo_m(:, np - ndr2/2 + 1:np), &
        !                       u_halo_p=u_halo_p, &
        !                       f=du2, &
        !                       L=self%der2%thomas3%lhs(:, 1:1))

        call self%der1%thomas3%solveU(du1)
        call self%der2%thomas3%solveU(du2)
        call self%der1%thomas3%reduce(du1, wrk2d(:, 1), wrk2d(:, 2))
        call self%der2%thomas3%reduce(du2, wrk2d(:, 1), wrk2d(:, 2))

        return
    end subroutine

end module FDM_Derivative_MPISplit
