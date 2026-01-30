#include "tlab_error.h"

! Split algorithm assumes periodic case over uniform grid
module FDM_Derivative_MPISplit
    use TLab_Constants, only: wp, wi, efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, stagger_on
    use TLabMPI_VARS
    use FDM_Derivative_Base
    use Thomas
    use Thomas_Split
    implicit none
    private

    ! -----------------------------------------------------------------------
    type, public :: der_periodic_mpisplit
        ! procedure(matmul_halo_ice), pointer, nopass :: matmul
        procedure(matmul_halo_thomas_ice), pointer, nopass :: matmul
        real(wp), pointer :: rhs(:, :) => null()
        type(thomas3_split_dt) thomas3
    contains
        procedure :: initialize => der_periodic_mpisplit_initialize
        procedure :: compute => der_periodic_mpisplit_compute
    end type

contains
    ! ###################################################################
    ! ###################################################################
    subroutine der_periodic_mpisplit_initialize(self, ref, axis)
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
    subroutine der_periodic_mpisplit_compute(self, nlines, u, u_halo_m, u_halo_p, result)
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

end module FDM_Derivative_MPISplit
