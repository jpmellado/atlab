#include "tlab_error.h"

! Split algorithm assumes periodic case over uniform grid
module FDM_Derivative_MPISplit
    use TLab_Constants, only: wp, wi, efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLabMPI_VARS
    use FDM_Derivative_Base, only: der_dt, der_periodic
    use FDM_Derivative_Base, only: matmul_halo_thomas_ice, matmul_halo_thomas_combined_ice
    use FDM_Derivative_1order, only: der1_periodic
    use FDM_Derivative_2order, only: der_extended_dt, der2_extended_periodic
    use Thomas
    use Thomas_Split
    implicit none
    private

    public :: der_periodic_mpisplit
    public :: der_burgers_mpisplit

    ! -----------------------------------------------------------------------
    type :: der_mpisplit
        ! procedure(matmul_halo_ice), pointer, nopass :: matmul
        procedure(matmul_halo_thomas_ice), pointer, nopass :: matmul
        real(wp), pointer :: rhs(:, :) => null()
        type(thomas3_split_dt) thomas3
    contains
    end type

    type, extends(der_mpisplit) :: der_periodic_mpisplit
    contains
        procedure :: initialize => der_periodic_initialize
        procedure :: compute => der_periodic_compute
    end type

    type :: der_burgers_mpisplit
        procedure(matmul_halo_thomas_combined_ice), pointer, nopass :: matmul => null()
        type(der_periodic_mpisplit) der1
        type(der_periodic_mpisplit) der2
    contains
        procedure :: initialize => der_burgers_periodic_initialize
        procedure :: compute => der_burgers_periodic_compute
    end type

contains
    ! ###################################################################
    ! ###################################################################
    subroutine der_periodic_initialize(self, ref, axis)
        class(der_periodic_mpisplit), intent(out) :: self
        class(der_periodic), intent(in), target :: ref
        character(len=*), intent(in) :: axis

        integer nx, k
        real(wp), allocatable :: lhs_loc(:, :)

        ! ###################################################################
        if (size(ref%lhs, 2) /= 3) then
            call TLab_Write_ASCII(efile, __FILE__//'. Matrix splitting only for tridiagonal.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        self%matmul => ref%matmul
        self%rhs => ref%rhs

        select case (trim(adjustl(axis)))
        case ('x')
            self%thomas3%communicator = ims_comm_x
            self%thomas3%rank = ims_pro_i
            self%thomas3%n_ranks = ims_npro_i

        case ('y')
            self%thomas3%communicator = ims_comm_y
            self%thomas3%rank = ims_pro_j
            self%thomas3%n_ranks = ims_npro_j

        end select

        self%thomas3%circulant = .true.
        self%thomas3%block_id = self%thomas3%rank + 1

        nx = size(ref%lhs, 1)
        allocate (lhs_loc, source=ref%lhs)
        call Thomas_Split_3_Initialize(lhs_loc(:, 1:1), lhs_loc(:, 2:3), &
                                       [(k, k=nx/self%thomas3%n_ranks, nx, nx/self%thomas3%n_ranks)], &
                                       self%thomas3)

        deallocate (lhs_loc)

        return
    end subroutine

    ! ###################################################################
    ! ###################################################################
    subroutine der_periodic_compute(self, nlines, u, u_halo_m, u_halo_p, result)
        use TLab_Arrays, only: wrk2d
        class(der_periodic_mpisplit), intent(in) :: self
        integer(wi), intent(in) :: nlines
        real(wp), intent(in) :: u(nlines, size(self%thomas3%lhs, 1))
        real(wp), intent(in) :: u_halo_m(:, :)
        real(wp), intent(in) :: u_halo_p(:, :)
        real(wp), intent(out) :: result(nlines, size(self%thomas3%lhs, 1))

        ! ###################################################################
        ! call self%matmul(self%rhs, u, u_halo_m, u_halo_p, result)
        call self%matmul(rhs=self%rhs(1, :), &
                         u=u, &
                         u_halo_m=u_halo_m, &
                         u_halo_p=u_halo_p, &
                         f=result, &
                         L=self%thomas3%lhs(:, 1:1))

        ! call Thomas3_SolveL(self%thomas3%lhs(:, 1:1), result)
        call Thomas3_SolveU(self%thomas3%lhs(:, 2:3), result)
        call ThomasSplit_3_Reduce_MPI(self%thomas3, result, wrk2d(:, 1), wrk2d(:, 2))

        return
    end subroutine

    ! ###################################################################
    ! ###################################################################
    subroutine der_burgers_periodic_initialize(self, ref1, ref2, axis)
        use FDM_Base, only: FDM_COM6_JACOBIAN, FDM_COM6_JACOBIAN_HYPER
        use Matmul_Halo_Thomas, only: MatMul_Halo_5_antisym_7_sym_ThomasL_3
        class(der_burgers_mpisplit), intent(out) :: self
        class(der_dt), intent(in) :: ref1
        class(der_extended_dt), intent(in) :: ref2
        character(len=*), intent(in) :: axis

        select type (ref1)
        type is (der1_periodic)
            call der_periodic_initialize(self%der1, ref1, axis)
        end select

        select type (ref2)
        type is (der2_extended_periodic)
            call der_periodic_initialize(self%der2, ref2%der2, axis)
        end select

        if (ref1%type == FDM_COM6_JACOBIAN .and. ref2%type == FDM_COM6_JACOBIAN_HYPER) then
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
        real(wp), intent(in) :: u(nlines, size(self%der1%thomas3%lhs, 1))
        real(wp), intent(in) :: u_halo_m(:, :)
        real(wp), intent(in) :: u_halo_p(:, :)
        real(wp), intent(out) :: du1(nlines, size(self%der1%thomas3%lhs, 1))
        real(wp), intent(out) :: du2(nlines, size(self%der1%thomas3%lhs, 1))

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
                         L1=self%der1%thomas3%lhs(:, 1:1), &
                         g=du2, &
                         L2=self%der2%thomas3%lhs(:, 1:1))

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

        call Thomas3_SolveU(self%der1%thomas3%lhs(:, 2:3), du1)
        call Thomas3_SolveU(self%der2%thomas3%lhs(:, 2:3), du2)
        call ThomasSplit_3_Reduce_MPI(self%der1%thomas3, du1, wrk2d(:, 1), wrk2d(:, 2))
        call ThomasSplit_3_Reduce_MPI(self%der2%thomas3, du2, wrk2d(:, 1), wrk2d(:, 2))

        return
    end subroutine

end module FDM_Derivative_MPISplit
