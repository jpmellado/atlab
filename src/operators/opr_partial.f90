#include "tlab_error.h"

module OPR_Partial
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_NONE
    use TLab_Arrays, only: wrk2d, wrk3d
    use TLab_Transpose
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_j
    use TLabMPI_Transpose
    use FDM_Derivative_MPISplit
#endif
    use TLab_Grid, only: x, y, z
    use FDM, only: g
    use FDM, only: fdm_der1_X, fdm_der1_Y, fdm_der1_Z
    use FDM_Derivative, only: FDM_Der1_Solve, FDM_Der2_Solve
    use Thomas_Split
    implicit none
    private

    public :: OPR_Partial_Initialize
    public :: OPR_Partial_X     ! These first 2 could be written in terms of the last one...
    public :: OPR_Partial_Y
    public :: OPR_Partial_Z

    integer, parameter, public :: OPR_P1 = 1                ! 1. order derivative
    integer, parameter, public :: OPR_P2 = 2                ! 2. order derivative
    integer, parameter, public :: OPR_P2_P1 = 3             ! 2. and 1.order derivatives
    integer, parameter, public :: OPR_P1_ADD = 4
    integer, parameter, public :: OPR_P1_SUBTRACT = 5

    ! -----------------------------------------------------------------------
    procedure(OPR_Partial_interface) :: OPR_Partial_dt
    abstract interface
        subroutine OPR_Partial_interface(type, nx, ny, nz, u, result, tmp1, ibc)
            use TLab_Constants, only: wi, wp
            integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
            integer(wi), intent(in) :: nx, ny, nz
            real(wp), intent(in) :: u(nx*ny*nz)
            real(wp), intent(out) :: result(nx*ny*nz)
            real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
            integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative
        end subroutine
    end interface
    procedure(OPR_Partial_dt), pointer :: OPR_Partial_X, OPR_Partial_Y

#ifdef USE_MPI
    type(fdm_derivative_split_dt), public, protected :: der1_split_x, der2_split_x
    type(fdm_derivative_split_dt), public, protected :: der1_split_y, der2_split_y
    real(wp), allocatable, target :: halo_m(:), halo_p(:)
    real(wp), pointer, public :: pyz_halo_m(:, :) => null(), pyz_halo_p(:, :) => null()
    real(wp), pointer, public :: pxz_halo_m(:, :) => null(), pxz_halo_p(:, :) => null()

    integer, public :: der_mode_i, der_mode_j
    integer, parameter, public :: TYPE_TRANSPOSE = 1
    integer, parameter, public :: TYPE_SPLIT = 2

#endif

contains
    ! ###################################################################
    ! ###################################################################
    subroutine OPR_Partial_Initialize(inifile)
        use TLab_Constants, only: efile
#ifdef USE_MPI
        use TLab_Memory, only: imax, jmax, kmax
        use TLab_Memory, only: TLab_Allocate_Real
#endif
        use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop

        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
#ifdef USE_MPI
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        integer np
#endif

#ifdef USE_MPI
        ! #######################################################################
        ! Read data
        bakfile = trim(adjustl(inifile))//'.bak'

        block = 'Parallel'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call ScanFile_Char(bakfile, inifile, block, 'DerivativeModeI', 'split', sRes)
        if (trim(adjustl(sRes)) == 'transpose') then; der_mode_i = TYPE_TRANSPOSE
        elseif (trim(adjustl(sRes)) == 'split') then; der_mode_i = TYPE_SPLIT
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong DerivativeModeI option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'DerivativeModeJ', 'split', sRes)
        if (trim(adjustl(sRes)) == 'transpose') then; der_mode_j = TYPE_TRANSPOSE
        elseif (trim(adjustl(sRes)) == 'split') then; der_mode_j = TYPE_SPLIT
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong DerivativeModeJ option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if
#endif

        ! ###################################################################
        ! Setting procedure pointers
#ifdef USE_MPI
        np = 0

        if (ims_npro_i > 1) then
            select case (der_mode_i)
            case (TYPE_TRANSPOSE)
                OPR_Partial_X => OPR_Partial_X_MPITranspose
            case (TYPE_SPLIT)
                OPR_Partial_X => OPR_Partial_X_MPISplit
                call FDM_MPISplit_Initialize(1, g(1)%der1, der1_split_x, 'x')
                call FDM_MPISplit_Initialize(2, g(1)%der2, der2_split_x, 'x')
                np = max(np, size(der1_split_x%rhs)/2)
                np = max(np, size(der2_split_x%rhs)/2)
            end select

        else
#endif
            OPR_Partial_X => OPR_Partial_X_Serial
#ifdef USE_MPI
        end if
#endif

#ifdef USE_MPI
        if (ims_npro_j > 1) then
            select case (der_mode_j)
            case (TYPE_TRANSPOSE)
                OPR_Partial_Y => OPR_Partial_Y_MPITranspose
            case (TYPE_SPLIT)
                OPR_Partial_Y => OPR_Partial_Y_MPISplit
                call FDM_MPISplit_Initialize(1, g(2)%der1, der1_split_y, 'y')
                call FDM_MPISplit_Initialize(2, g(2)%der2, der2_split_y, 'y')
                np = max(np, size(der1_split_y%rhs)/2)
                np = max(np, size(der2_split_y%rhs)/2)
            end select

        else
#endif
            OPR_Partial_Y => OPR_Partial_Y_Serial
#ifdef USE_MPI
        end if

        if (np > 0) then
            allocate (halo_m(max(imax*kmax, jmax*kmax)*np))
            pyz_halo_m(1:jmax*kmax, 1:np) => halo_m(1:jmax*kmax*np)
            pxz_halo_m(1:imax*kmax, 1:np) => halo_m(1:imax*kmax*np)

            allocate (halo_p(max(imax*kmax, jmax*kmax)*np))
            pyz_halo_p(1:jmax*kmax, 1:np) => halo_p(1:jmax*kmax*np)
            pxz_halo_p(1:imax*kmax, 1:np) => halo_p(1:imax*kmax*np)
        end if

#endif

        return
    end subroutine OPR_Partial_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine OPR_Partial_X_Serial(type, nx, ny, nz, u, result, tmp1, ibc)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
        integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative

        ! -------------------------------------------------------------------
        integer ibc_loc

        ! ###################################################################
        if (x%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (type == OPR_P2_P1) tmp1 = 0.0_wp
            return
        end if

        ! Transposition: make x-direction the last one
#ifdef USE_ESSL
        call DGETMO(u, nx, nx, ny*nz, result, ny*nz)
#else
        call TLab_Transpose_Real(u, nx, ny*nz, nx, result, ny*nz, locBlock=trans_x_forward)
#endif

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        select case (type)
        case (OPR_P2)
            ! if (g(1)%der2%need_1der) call FDM_Der1_Solve(ny*nz, g(1)%der1, g(1)%der1%lu, result, tmp1, wrk2d, ibc_loc)
            if (.not. x%uniform) call fdm_der1_X%compute(ny*nz, result, tmp1)
            call FDM_Der2_Solve(ny*nz, g(1)%der2, g(1)%der2%lu, result, wrk3d, tmp1, wrk2d)

        case (OPR_P2_P1)
            ! call FDM_Der1_Solve(ny*nz, g(1)%der1, g(1)%der1%lu, result, wrk3d, wrk2d, ibc_loc)
            call fdm_der1_X%compute(ny*nz, result, wrk3d)
            call FDM_Der2_Solve(ny*nz, g(1)%der2, g(1)%der2%lu, result, tmp1, wrk3d, wrk2d)

        case (OPR_P1, OPR_P1_ADD, OPR_P1_SUBTRACT)
            ! call FDM_Der1_Solve(ny*nz, g(1)%der1, g(1)%der1%lu, result, wrk3d, wrk2d, ibc_loc)
            call fdm_der1_X%compute(ny*nz, result, wrk3d)

        end select

        ! Put arrays back in the order in which they came in
        select case (type)
        case (OPR_P2_P1)
            call TLab_Transpose_Real(tmp1, ny*nz, nx, ny*nz, result, nx, locBlock=trans_x_backward)
            call TLab_Transpose_Real(wrk3d, ny*nz, nx, ny*nz, tmp1, nx, locBlock=trans_x_backward)

        case (OPR_P1_ADD)
            call TLab_AddTranspose(wrk3d, ny*nz, nx, ny*nz, tmp1, nx, locBlock=trans_x_backward)

        case (OPR_P1_SUBTRACT)
            call TLab_SubtractTranspose(wrk3d, ny*nz, nx, ny*nz, tmp1, nx, locBlock=trans_x_backward)

        case default
            call TLab_Transpose_Real(wrk3d, ny*nz, nx, ny*nz, result, nx, locBlock=trans_x_backward)

        end select

        return
    end subroutine OPR_Partial_X_Serial

    ! ###################################################################
    ! ###################################################################
#ifdef USE_MPI
    subroutine OPR_Partial_X_MPITranspose(type, nx, ny, nz, u, result, tmp1, ibc)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
        integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative

        ! -------------------------------------------------------------------
        integer(wi) nlines
        integer ibc_loc

        ! ###################################################################
        if (x%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (type == OPR_P2_P1) tmp1 = 0.0_wp
            return
        end if

        nlines = tmpi_plan_dx%nlines

        ! Transposition: make x-direction the last one
        call TLabMPI_Trp_ExecI_Forward(u, result, tmpi_plan_dx)
#ifdef USE_ESSL
        call DGETMO(result, g(1)%size, g(1)%size, nlines, wrk3d, nlines)
#else
        call TLab_Transpose_Real(result, g(1)%size, nlines, g(1)%size, wrk3d, nlines)
#endif

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        select case (type)
        case (OPR_P2)
            ! if (g(1)%der2%need_1der) call FDM_Der1_Solve(nlines, g(1)%der1, g(1)%der1%lu, wrk3d, tmp1, wrk2d, ibc_loc)
            if (.not. x%uniform) call fdm_der1_X%compute(nlines, wrk3d, tmp1)
            call FDM_Der2_Solve(nlines, g(1)%der2, g(1)%der2%lu, wrk3d, result, tmp1, wrk2d)

        case (OPR_P2_P1)
            ! call FDM_Der1_Solve(nlines, g(1)%der1, g(1)%der1%lu, wrk3d, tmp1, wrk2d, ibc_loc)
            call fdm_der1_X%compute(nlines, wrk3d, tmp1)
            call FDM_Der2_Solve(nlines, g(1)%der2, g(1)%der2%lu, wrk3d, result, tmp1, wrk2d)

        case (OPR_P1, OPR_P1_ADD, OPR_P1_SUBTRACT)
            ! call FDM_Der1_Solve(nlines, g(1)%der1, g(1)%der1%lu, wrk3d, result, wrk2d, ibc_loc)
            call fdm_der1_X%compute(nlines, wrk3d, result)

        end select

        ! Put arrays back in the order in which they came in
        select case (type)
        case (OPR_P2_P1)
            call TLab_Transpose_Real(tmp1, nlines, g(1)%size, nlines, wrk3d, g(1)%size)
            call TLab_Transpose_Real(result, nlines, g(1)%size, nlines, tmp1, g(1)%size)
            call TLabMPI_Trp_ExecI_Backward(tmp1, result, tmpi_plan_dx)
            call TLabMPI_Trp_ExecI_Backward(wrk3d, tmp1, tmpi_plan_dx)

        case (OPR_P1_ADD)
            call TLab_Transpose_Real(result, nlines, g(1)%size, nlines, wrk3d, g(1)%size)
            call TLabMPI_Trp_ExecI_Backward(wrk3d, result, tmpi_plan_dx)
            tmp1 = tmp1 + result

        case (OPR_P1_SUBTRACT)
            call TLab_Transpose_Real(result, nlines, g(1)%size, nlines, wrk3d, g(1)%size)
            call TLabMPI_Trp_ExecI_Backward(wrk3d, result, tmpi_plan_dx)
            tmp1 = tmp1 - result

        case default
            call TLab_Transpose_Real(result, nlines, g(1)%size, nlines, wrk3d, g(1)%size)
            call TLabMPI_Trp_ExecI_Backward(wrk3d, result, tmpi_plan_dx)

        end select

        return
    end subroutine OPR_Partial_X_MPITranspose

    !########################################################################
    !########################################################################
    subroutine OPR_Partial_X_MPISplit(type, nx, ny, nz, u, result, tmp1, ibc)
        use TLabMPI_PROCS, only: TLabMPI_Halos_X
        ! use TLab_Pointers_2D, only: pyz_wrk3d
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
        integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative

        ! -------------------------------------------------------------------
        integer np, np1, np2

        ! ###################################################################
        if (x%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (type == OPR_P2_P1) tmp1 = 0.0_wp
            return
        end if

        ! Transposition: make x-direction the last one
#ifdef USE_ESSL
        call DGETMO(u, nx, nx, ny*nz, result, ny*nz)
#else
        call TLab_Transpose_Real(u, nx, ny*nz, nx, result, ny*nz, locBlock=trans_x_forward)
#endif

        np1 = size(der1_split_x%rhs)/2
        np2 = size(der2_split_x%rhs)/2
        np = max(np1, np2)
        call TLabMPI_Halos_X(result, ny*nz, np, pyz_halo_m(:, 1), pyz_halo_p(:, 1))

        select case (type)
        case (OPR_P2)
            call FDM_MPISplit_Solve(ny*nz, nx, der2_split_x, result, &
                                    pyz_halo_m(:, np - np2 + 1:np), pyz_halo_p, wrk3d, wrk2d)

        case (OPR_P2_P1)
            call FDM_MPISplit_Solve(ny*nz, nx, der2_split_x, result, &
                                    pyz_halo_m(:, np - np2 + 1:np), pyz_halo_p, tmp1, wrk2d)
            call FDM_MPISplit_Solve(ny*nz, nx, der1_split_x, result, &
                                    pyz_halo_m(:, np - np1 + 1:np), pyz_halo_p, wrk3d, wrk2d)

        case (OPR_P1, OPR_P1_ADD, OPR_P1_SUBTRACT)
            call FDM_MPISplit_Solve(ny*nz, nx, der1_split_x, result, &
                                    pyz_halo_m(:, np - np1 + 1:np), pyz_halo_p, wrk3d, wrk2d)

        end select

        ! Put arrays back in the order in which they came in
        select case (type)
        case (OPR_P2_P1)
            call TLab_Transpose_Real(tmp1, ny*nz, nx, ny*nz, result, nx, locBlock=trans_x_backward)
            call TLab_Transpose_Real(wrk3d, ny*nz, nx, ny*nz, tmp1, nx, locBlock=trans_x_backward)

        case (OPR_P1_ADD)
            call TLab_AddTranspose(wrk3d, ny*nz, nx, ny*nz, tmp1, nx, locBlock=trans_x_backward)

        case (OPR_P1_SUBTRACT)
            call TLab_SubtractTranspose(wrk3d, ny*nz, nx, ny*nz, tmp1, nx, locBlock=trans_x_backward)

        case default
            call TLab_Transpose_Real(wrk3d, ny*nz, nx, ny*nz, result, nx, locBlock=trans_x_backward)

        end select

        return
    end subroutine OPR_Partial_X_MPISplit

#endif

    !########################################################################
    !########################################################################
    subroutine OPR_Partial_Y_Serial(type, nx, ny, nz, u, result, tmp1, ibc)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
        integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative

        ! -------------------------------------------------------------------
        integer(wi) nlines
        integer ibc_loc

        ! ###################################################################
        if (y%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (type == OPR_P2_P1) tmp1 = 0.0_wp
            return
        end if

        ! Transposition: make y-direction the last one
#ifdef USE_ESSL
        call DGETMO(u, nx*ny, nx*ny, nz, result, nz)
#else
        call TLab_Transpose_Real(u, nx*ny, nz, nx*ny, result, nz, locBlock=trans_y_forward)
#endif
        nlines = nx*nz

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        select case (type)
        case (OPR_P2)
            ! if (g(2)%der2%need_1der) call FDM_Der1_Solve(nlines, g(2)%der1, g(2)%der1%lu, result, tmp1, wrk2d, ibc_loc)
            if (.not. y%uniform) call fdm_der1_Y%compute(nlines, result, tmp1)
            call FDM_Der2_Solve(nlines, g(2)%der2, g(2)%der2%lu, result, wrk3d, tmp1, wrk2d)

        case (OPR_P2_P1)
            ! call FDM_Der1_Solve(nlines, g(2)%der1, g(2)%der1%lu, result, wrk3d, wrk2d, ibc_loc)
            call fdm_der1_Y%compute(nlines, result, wrk3d)
            call FDM_Der2_Solve(nlines, g(2)%der2, g(2)%der2%lu, result, tmp1, wrk3d, wrk2d)

        case (OPR_P1, OPR_P1_ADD, OPR_P1_SUBTRACT)
            ! call FDM_Der1_Solve(nlines, g(2)%der1, g(2)%der1%lu, result, wrk3d, wrk2d, ibc_loc)
            call fdm_der1_Y%compute(nlines, result, wrk3d)

        end select

        ! Put arrays back in the order in which they came in
        select case (type)
        case (OPR_P2_P1)
            call TLab_Transpose_Real(tmp1, nz, nx*ny, nz, result, nx*ny, locBlock=trans_y_backward)
            call TLab_Transpose_Real(wrk3d, nz, nx*ny, nz, tmp1, nx*ny, locBlock=trans_y_backward)

        case (OPR_P1_ADD)
            call TLab_AddTranspose(wrk3d, nz, nx*ny, nz, tmp1, nx*ny, locBlock=trans_y_backward)

        case (OPR_P1_SUBTRACT)
            call TLab_SubtractTranspose(wrk3d, nz, nx*ny, nz, tmp1, nx*ny, locBlock=trans_y_backward)

        case default
            call TLab_Transpose_Real(wrk3d, nz, nx*ny, nz, result, nx*ny, locBlock=trans_y_backward)

        end select

        return
    end subroutine OPR_Partial_Y_Serial

    !########################################################################
    !########################################################################
#ifdef USE_MPI
    subroutine OPR_Partial_Y_MPITranspose(type, nx, ny, nz, u, result, tmp1, ibc)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
        integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative

        ! -------------------------------------------------------------------
        integer(wi) nlines
        integer ibc_loc

        ! ###################################################################
        if (y%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (type == OPR_P2_P1) tmp1 = 0.0_wp
            return
        end if

        ! Transposition: make y-direction the last one
#ifdef USE_ESSL
        call DGETMO(u, nx*ny, nx*ny, nz, result, nz)
#else
        call TLab_Transpose_Real(u, nx*ny, nz, nx*ny, result, nz)
#endif
        call TLabMPI_Trp_ExecJ_Forward(result, wrk3d, tmpi_plan_dy)
        nlines = tmpi_plan_dy%nlines

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        select case (type)
        case (OPR_P2)
!            if (g(2)%der2%need_1der) call FDM_Der1_Solve(nlines, g(2)%der1, g(2)%der1%lu, wrk3d, tmp1, wrk2d, ibc_loc)
            if (.not. y%uniform) call fdm_der1_Y%compute(nlines, wrk3d, tmp1)
            call FDM_Der2_Solve(nlines, g(2)%der2, g(2)%der2%lu, wrk3d, result, tmp1, wrk2d)

        case (OPR_P2_P1)
            ! call FDM_Der1_Solve(nlines, g(2)%der1, g(2)%der1%lu, wrk3d, tmp1, wrk2d, ibc_loc)
            call fdm_der1_Y%compute(nlines, wrk3d, tmp1)
            call FDM_Der2_Solve(nlines, g(2)%der2, g(2)%der2%lu, wrk3d, result, tmp1, wrk2d)

        case (OPR_P1, OPR_P1_ADD, OPR_P1_SUBTRACT)
            ! call FDM_Der1_Solve(nlines, g(2)%der1, g(2)%der1%lu, wrk3d, result, wrk2d, ibc_loc)
            call fdm_der1_Y%compute(nlines, wrk3d, result)

        end select

        ! Put arrays back in the order in which they came in
        select case (type)
        case (OPR_P2_P1)
            call TLabMPI_Trp_ExecJ_Backward(tmp1, wrk3d, tmpi_plan_dy)
            call TLabMPI_Trp_ExecJ_Backward(result, tmp1, tmpi_plan_dy)
            call TLab_Transpose_Real(tmp1, nz, nx*ny, nz, result, nx*ny)
            call TLab_Transpose_Real(wrk3d, nz, nx*ny, nz, tmp1, nx*ny)

        case (OPR_P1_ADD)
            call TLabMPI_Trp_ExecJ_Backward(result, wrk3d, tmpi_plan_dy)
            call TLab_AddTranspose(wrk3d, nz, nx*ny, nz, tmp1, nx*ny)

        case (OPR_P1_SUBTRACT)
            call TLabMPI_Trp_ExecJ_Backward(result, wrk3d, tmpi_plan_dy)
            call TLab_SubtractTranspose(wrk3d, nz, nx*ny, nz, tmp1, nx*ny)

        case default
            call TLabMPI_Trp_ExecJ_Backward(result, wrk3d, tmpi_plan_dy)
            call TLab_Transpose_Real(wrk3d, nz, nx*ny, nz, result, nx*ny)

        end select

        return
    end subroutine OPR_Partial_Y_MPITranspose

    !########################################################################
    !########################################################################
    subroutine OPR_Partial_Y_MPISplit(type, nx, ny, nz, u, result, tmp1, ibc)
        use TLabMPI_PROCS, only: TLabMPI_Halos_Y
        use TLab_Pointers_2D, only: pxz_wrk3d
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
        integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative

        ! -------------------------------------------------------------------
        integer np, np1, np2

        ! ###################################################################
        if (y%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (type == OPR_P2_P1) tmp1 = 0.0_wp
            return
        end if

        ! Transposition: make y-direction the last one
#ifdef USE_ESSL
        call DGETMO(u, nx*ny, nx*ny, nz, result, nz)
#else
        call TLab_Transpose_Real(u, nx*ny, nz, nx*ny, result, nz, locBlock=trans_y_forward)
#endif

        np1 = size(der1_split_y%rhs)/2
        np2 = size(der2_split_y%rhs)/2
        np = max(np1, np2)
        call TLabMPI_Halos_Y(result, nx*nz, np, pxz_halo_m(:, 1), pxz_halo_p(:, 1))

        select case (type)
        case (OPR_P2)
            call FDM_MPISplit_Solve(nx*nz, ny, der2_split_y, result, &
                                    pxz_halo_m(:, np - np2 + 1:np), pxz_halo_p, pxz_wrk3d, wrk2d)

        case (OPR_P2_P1)
            call FDM_MPISplit_Solve(nx*nz, ny, der2_split_y, result, &
                                    pxz_halo_m(:, np - np2 + 1:np), pxz_halo_p, tmp1, wrk2d)
            call FDM_MPISplit_Solve(nx*nz, ny, der1_split_y, result, &
                                    pxz_halo_m(:, np - np1 + 1:np), pxz_halo_p, pxz_wrk3d, wrk2d)

        case (OPR_P1, OPR_P1_ADD, OPR_P1_SUBTRACT)
            call FDM_MPISplit_Solve(nx*nz, ny, der1_split_y, result, &
                                    pxz_halo_m(:, np - np1 + 1:np), pxz_halo_p, pxz_wrk3d, wrk2d)

        end select

        ! Put arrays back in the order in which they came in
        select case (type)
        case (OPR_P2_P1)
            call TLab_Transpose_Real(tmp1, nz, nx*ny, nz, result, nx*ny, locBlock=trans_y_backward)
            call TLab_Transpose_Real(wrk3d, nz, nx*ny, nz, tmp1, nx*ny, locBlock=trans_y_backward)

        case (OPR_P1_ADD)
            call TLab_AddTranspose(wrk3d, nz, nx*ny, nz, tmp1, nx*ny, locBlock=trans_y_backward)

        case (OPR_P1_SUBTRACT)
            call TLab_SubtractTranspose(wrk3d, nz, nx*ny, nz, tmp1, nx*ny, locBlock=trans_y_backward)

        case default
            call TLab_Transpose_Real(wrk3d, nz, nx*ny, nz, result, nx*ny, locBlock=trans_y_backward)

        end select

        return
    end subroutine OPR_Partial_Y_MPISplit

#endif

    !########################################################################
    !########################################################################
    subroutine OPR_Partial_Z(type, nx, ny, nz, u, result, tmp1, ibc)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
        integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative

        ! -------------------------------------------------------------------
        integer ibc_loc

        ! ###################################################################
        if (z%size == 1) then
            result = 0.0_wp
            if (type == OPR_P2_P1) tmp1 = 0.0_wp
            return
        end if

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        select case (type)
        case (OPR_P2)
            ! if (g(3)%der2%need_1der) call FDM_Der1_Solve(nx*ny, g(3)%der1, g(3)%der1%lu, u, wrk3d, wrk2d, ibc_loc)
            if (.not. z%uniform) then
                if (ibc_loc == BCS_NONE) then       ! testing
                    call fdm_der1_Z%compute(nx*ny, u, wrk3d)
                else
                    call FDM_Der1_Solve(nx*ny, g(3)%der1, g(3)%der1%lu, u, wrk3d, wrk2d, ibc_loc)
                end if
            end if
            call FDM_Der2_Solve(nx*ny, g(3)%der2, g(3)%der2%lu, u, result, wrk3d, wrk2d)

        case (OPR_P2_P1)
            if (ibc_loc == BCS_NONE) then       ! testing
                call fdm_der1_Z%compute(nx*ny, u, tmp1)
            else
                call FDM_Der1_Solve(nx*ny, g(3)%der1, g(3)%der1%lu, u, tmp1, wrk2d, ibc_loc)
            end if
            call FDM_Der2_Solve(nx*ny, g(3)%der2, g(3)%der2%lu, u, result, tmp1, wrk2d)

        case (OPR_P1)
            if (ibc_loc == BCS_NONE) then       ! testing
                call fdm_der1_Z%compute(nx*ny, u, result)
            else
                call FDM_Der1_Solve(nx*ny, g(3)%der1, g(3)%der1%lu, u, result, wrk2d, ibc_loc)
            end if

        end select

        return
    end subroutine OPR_Partial_Z

end module OPR_Partial
