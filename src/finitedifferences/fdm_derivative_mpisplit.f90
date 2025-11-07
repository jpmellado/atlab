#include "tlab_error.h"

! Split algorithm assumes periodic case
module FDM_Derivative_MPISplit
    use TLab_Constants, only: wp, wi, efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, stagger_on
    use TLabMPI_VARS
    use Thomas
    use Thomas_Split
    use Matmul_Halo
    use Matmul_Halo_Thomas
    implicit none
    private

    public :: FDM_MPISplit_Initialize
    public :: FDM_MPISplit_Solve

    type, public :: fdm_derivative_split_dt
        real(wp), allocatable :: rhs(:)
        procedure(matmul_halo_ice), pointer, nopass :: matmul_halo
        procedure(matmul_halo_thomas_ice), pointer, nopass :: matmul_halo_thomas
        type(thomas3_split_dt) thomas3
    end type fdm_derivative_split_dt

    ! -----------------------------------------------------------------------
    abstract interface
        subroutine matmul_halo_ice(rhs, u, u_halo_m, u_halo_p, f)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: rhs(:)              ! diagonals of B
            real(wp), intent(in) :: u(:, :)             ! vector u
            real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
            real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
            real(wp), intent(out) :: f(:, :)            ! vector f = B u
        end subroutine
    end interface

    abstract interface
        subroutine matmul_halo_thomas_ice(rhs, u, u_halo_m, u_halo_p, f, L)
            use TLab_Constants, only: wp
            real(wp), intent(in) :: rhs(:)              ! diagonals of B
            real(wp), intent(in) :: u(:, :)             ! vector u
            real(wp), intent(in) :: u_halo_m(:, :)      ! minus, coming from left
            real(wp), intent(in) :: u_halo_p(:, :)      ! plus, coming from right
            real(wp), intent(out) :: f(:, :)            ! vector f = B u
            real(wp), intent(in) :: L(:, :)
        end subroutine
    end interface

contains
    ! ###################################################################
    ! ###################################################################
    subroutine FDM_MPISplit_Initialize(order, g, gSplit, axis)
        use FDM_Derivative, only: fdm_derivative_dt
        integer, intent(in) :: order                            ! order of the derivative
        type(fdm_derivative_dt), intent(in) :: g                ! original information about derivative
        type(fdm_derivative_split_dt), intent(out) :: gSplit    ! split information
        character(len=*), intent(in) :: axis

        integer(wi) nsize, nd
        integer k
        real(wp), allocatable :: lhs_loc(:, :)

        ! ###################################################################
        nsize = size(g%lhs, 1)
        nd = g%nb_diag(1)
        if (nd /= 3) then
            call TLab_Write_ASCII(efile, __FILE__//'. Matrix splitting only for tridiagonal.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)
        end if

        if (allocated(gSplit%rhs)) deallocate (gSplit%rhs)
        allocate (gSplit%rhs(g%nb_diag(2)))
        gSplit%rhs(:) = g%rhs(1, 1:g%nb_diag(2))    ! assuming periodic over uniform grid

        select case (trim(adjustl(axis)))
        case ('x')
            gSplit%thomas3%communicator = ims_comm_x
            gSplit%thomas3%rank = ims_pro_i
            gSplit%thomas3%n_ranks = ims_npro_i

        case ('y')
            gSplit%thomas3%communicator = ims_comm_y
            gSplit%thomas3%rank = ims_pro_j
            gSplit%thomas3%n_ranks = ims_npro_j

        end select

        gSplit%thomas3%circulant = .true.
        gSplit%thomas3%block_id = gSplit%thomas3%rank + 1

        allocate (lhs_loc(nsize, nd))
        lhs_loc(1:nsize, 1:nd) = g%lhs(1:nsize, 1:nd)
        call Thomas_Split_3_Initialize(lhs_loc(:, 1:1), lhs_loc(:, 2:3), &
                                       [(k, k=nsize/gSplit%thomas3%n_ranks, nsize, nsize/gSplit%thomas3%n_ranks)], &
                                       gSplit%thomas3)

        select case (g%nb_diag(2))
        case (3)
            ! if (order == 1) gSplit%matmul_halo => MatMul_Halo_3d_antisym
            ! if (order == 2) gSplit%matmul_halo => MatMul_Halo_3d_sym
            if (order == 1) gSplit%matmul_halo_thomas => MatMul_Halo_3_antisym_ThomasL_3
            if (order == 2) gSplit%matmul_halo_thomas => MatMul_Halo_3_sym_ThomasL_3

        case (5)
            ! if (order == 1) gSplit%matmul_halo => MatMul_Halo_5d_antisym
            ! if (order == 2) gSplit%matmul_halo => MatMul_Halo_5d_sym
            if (order == 1) gSplit%matmul_halo_thomas => MatMul_Halo_5_antisym_ThomasL_3
            if (order == 2) gSplit%matmul_halo_thomas => MatMul_Halo_5_sym_ThomasL_3

        case (7)
            ! if (order == 1) gSplit%matmul_halo => MatMul_Halo_7d_antisym
            ! if (order == 2) gSplit%matmul_halo => MatMul_Halo_7d_sym
            if (order == 1) gSplit%matmul_halo_thomas => MatMul_Halo_7_antisym_ThomasL_3
            if (order == 2) gSplit%matmul_halo_thomas => MatMul_Halo_7_sym_ThomasL_3

        case default
            call TLab_Write_ASCII(efile, __FILE__//'. Undeveloped number of diagonals.')
            call TLab_Stop(DNS_ERROR_UNDEVELOP)

        end select

        deallocate (lhs_loc)

        return
    end subroutine FDM_MPISplit_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine FDM_MPISplit_Solve(nlines, nsize, gSplit, u, u_halo_m, u_halo_p, result, wrk2d)
        integer(wi), intent(in) :: nlines                   ! # of lines to be solved
        integer(wi), intent(in) :: nsize                    ! local size of the system
        type(fdm_derivative_split_dt), intent(in) :: gSplit
        real(wp), intent(in) :: u(nlines, nsize)
        real(wp), intent(in) :: u_halo_m(:, :)
        real(wp), intent(in) :: u_halo_p(:, :)
        real(wp), intent(out) :: result(nlines, nsize)      ! derivative of u
        real(wp), intent(inout) :: wrk2d(nlines, 2)

        ! call gSplit%matmul_halo(gSplit%rhs, u, u_halo_m, u_halo_p, result)
        call gSplit%matmul_halo_thomas(gSplit%rhs, u, u_halo_m, u_halo_p, result, gSplit%thomas3%lhs(:, 1:1))

        ! call Thomas3_SolveL(gSplit%thomas3%lhs(:, 1:1), result)
        call Thomas3_SolveU(gSplit%thomas3%lhs(:, 2:3), result)
        call ThomasSplit_3_Reduce_MPI(gSplit%thomas3, result, wrk2d(:, 1), wrk2d(:, 2))

        return
    end subroutine FDM_MPISplit_Solve

end module FDM_Derivative_MPISplit
