#include "tlab_error.h"

module OPR_Partial
    use TLab_Constants, only: wp, wi
    use TLab_Arrays, only: wrk2d, wrk3d
    use TLab_Transpose
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_j
    use TLabMPI_Transpose
    use FDM_Derivative_MPISplit
    use FDM_Derivative_1order
    use FDM_Derivative_2order
#endif
    use TLab_Grid, only: x, y, z
    use FDM, only: fdm_der1_X, fdm_der1_Y, fdm_der1_Z
    use FDM, only: fdm_der2_X, fdm_der2_Y, fdm_der2_Z
    use Thomas_Split
    implicit none
    private

    public :: OPR_Partial_Initialize
    public :: OPR_Partial_X
    public :: OPR_Partial_Y
    public :: OPR_Partial_Z
    public :: OPR_Partial_Z_Bcs

    integer, parameter, public :: OPR_P1 = 1                ! 1. order derivative
    integer, parameter, public :: OPR_P2 = 2                ! 2. order derivative
    integer, parameter, public :: OPR_P2_P1 = 3             ! 2. and 1.order derivatives
    integer, parameter, public :: OPR_P1_ADD = 4
    integer, parameter, public :: OPR_P1_SUBTRACT = 5

    ! -----------------------------------------------------------------------
    procedure(OPR_Partial_interface) :: OPR_Partial_dt
    abstract interface
        subroutine OPR_Partial_interface(type, nx, ny, nz, u, result, tmp1)
            use TLab_Constants, only: wi, wp
            integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
            integer(wi), intent(in) :: nx, ny, nz
            real(wp), intent(in) :: u(nx*ny*nz)
            real(wp), intent(out) :: result(nx*ny*nz)
            real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
        end subroutine
    end interface
    procedure(OPR_Partial_dt), pointer :: OPR_Partial_X, OPR_Partial_Y

#ifdef USE_MPI
    type(der_periodic_mpisplit), public, protected :: fdm_der1_X_split, fdm_der2_X_split
    type(der_periodic_mpisplit), public, protected :: fdm_der1_Y_split, fdm_der2_Y_split
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
        if (trim(adjustl(sRes)) == 'transpose') then; der_mode_i = TYPE_TRANSPOSE
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
                select type (fdm_der1_X)
                type is (der1_periodic)
                    call fdm_der1_X_split%initialize(fdm_der1_X, 'x')
                end select
                np = max(np, size(fdm_der1_X_split%rhs, 2)/2)
                select type (fdm_der2_X)
                type is (der2_extended_periodic)
                    call fdm_der2_X_split%initialize(fdm_der2_X%der2, 'x')
                end select
                np = max(np, size(fdm_der2_X_split%rhs, 2)/2)
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
                select type (fdm_der1_Y)
                type is (der1_periodic)
                    call fdm_der1_Y_split%initialize(fdm_der1_Y, 'y')
                end select
                np = max(np, size(fdm_der1_Y_split%rhs, 2)/2)
                select type (fdm_der2_Y)
                type is (der2_extended_periodic)
                    call fdm_der2_Y_split%initialize(fdm_der2_Y%der2, 'y')
                end select
                np = max(np, size(fdm_der2_Y_split%rhs, 2)/2)
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
    subroutine OPR_Partial_X_Serial(type, nx, ny, nz, u, result, tmp1)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation

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

        select case (type)
        case (OPR_P2)
            if (.not. x%uniform) call fdm_der1_X%compute(ny*nz, result, tmp1)
            call fdm_der2_X%compute(ny*nz, result, wrk3d, tmp1)

        case (OPR_P2_P1)
            call fdm_der1_X%compute(ny*nz, result, wrk3d)
            call fdm_der2_X%compute(ny*nz, result, tmp1, wrk3d)

        case (OPR_P1, OPR_P1_ADD, OPR_P1_SUBTRACT)
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
    subroutine OPR_Partial_X_MPITranspose(type, nx, ny, nz, u, result, tmp1)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation

        ! -------------------------------------------------------------------
        integer(wi) nlines

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
        call DGETMO(result, x%size, x%size, nlines, wrk3d, nlines)
#else
        call TLab_Transpose_Real(result, x%size, nlines, x%size, wrk3d, nlines)
#endif

        select case (type)
        case (OPR_P2)
            if (.not. x%uniform) call fdm_der1_X%compute(nlines, wrk3d, tmp1)
            call fdm_der2_X%compute(ny*nz, wrk3d, result, tmp1)

        case (OPR_P2_P1)
            call fdm_der1_X%compute(nlines, wrk3d, tmp1)
            call fdm_der2_X%compute(ny*nz, wrk3d, result, tmp1)

        case (OPR_P1, OPR_P1_ADD, OPR_P1_SUBTRACT)
            call fdm_der1_X%compute(nlines, wrk3d, result)

        end select

        ! Put arrays back in the order in which they came in
        select case (type)
        case (OPR_P2_P1)
            call TLab_Transpose_Real(tmp1, nlines, x%size, nlines, wrk3d, x%size)
            call TLab_Transpose_Real(result, nlines, x%size, nlines, tmp1, x%size)
            call TLabMPI_Trp_ExecI_Backward(tmp1, result, tmpi_plan_dx)
            call TLabMPI_Trp_ExecI_Backward(wrk3d, tmp1, tmpi_plan_dx)

        case (OPR_P1_ADD)
            call TLab_Transpose_Real(result, nlines, x%size, nlines, wrk3d, x%size)
            call TLabMPI_Trp_ExecI_Backward(wrk3d, result, tmpi_plan_dx)
            tmp1 = tmp1 + result

        case (OPR_P1_SUBTRACT)
            call TLab_Transpose_Real(result, nlines, x%size, nlines, wrk3d, x%size)
            call TLabMPI_Trp_ExecI_Backward(wrk3d, result, tmpi_plan_dx)
            tmp1 = tmp1 - result

        case default
            call TLab_Transpose_Real(result, nlines, x%size, nlines, wrk3d, x%size)
            call TLabMPI_Trp_ExecI_Backward(wrk3d, result, tmpi_plan_dx)

        end select

        return
    end subroutine OPR_Partial_X_MPITranspose

    !########################################################################
    !########################################################################
    subroutine OPR_Partial_X_MPISplit(type, nx, ny, nz, u, result, tmp1)
        use TLabMPI_PROCS, only: TLabMPI_Halos_X
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation

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

        np1 = size(fdm_der1_X_split%rhs, 2)/2
        np2 = size(fdm_der2_X_split%rhs, 2)/2
        np = max(np1, np2)
        call TLabMPI_Halos_X(result, ny*nz, np, pyz_halo_m(:, 1), pyz_halo_p(:, 1))

        select case (type)
        case (OPR_P2)
            call fdm_der2_X_split%compute(ny*nz, result, pyz_halo_m(:, np - np2 + 1:np), pyz_halo_p, wrk3d)

        case (OPR_P2_P1)
            call fdm_der2_X_split%compute(ny*nz, result, pyz_halo_m(:, np - np2 + 1:np), pyz_halo_p, tmp1)
            call fdm_der1_X_split%compute(ny*nz, result, pyz_halo_m(:, np - np1 + 1:np), pyz_halo_p, wrk3d)

        case (OPR_P1, OPR_P1_ADD, OPR_P1_SUBTRACT)
            call fdm_der1_X_split%compute(ny*nz, result, pyz_halo_m(:, np - np1 + 1:np), pyz_halo_p, wrk3d)

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
    subroutine OPR_Partial_Y_Serial(type, nx, ny, nz, u, result, tmp1)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation

        ! -------------------------------------------------------------------
        integer(wi) nlines

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

        select case (type)
        case (OPR_P2)
            if (.not. y%uniform) call fdm_der1_Y%compute(nlines, result, tmp1)
            call fdm_der2_Y%compute(nlines, result, wrk3d, tmp1)

        case (OPR_P2_P1)
            call fdm_der1_Y%compute(nlines, result, wrk3d)
            call fdm_der2_Y%compute(nlines, result, tmp1, wrk3d)

        case (OPR_P1, OPR_P1_ADD, OPR_P1_SUBTRACT)
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
    subroutine OPR_Partial_Y_MPITranspose(type, nx, ny, nz, u, result, tmp1)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation

        ! -------------------------------------------------------------------
        integer(wi) nlines

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

        select case (type)
        case (OPR_P2)
            if (.not. y%uniform) call fdm_der1_Y%compute(nlines, wrk3d, tmp1)
            call fdm_der2_Y%compute(nlines, wrk3d, result, tmp1)

        case (OPR_P2_P1)
            call fdm_der1_Y%compute(nlines, wrk3d, tmp1)
            call fdm_der2_Y%compute(nlines, wrk3d, result, tmp1)

        case (OPR_P1, OPR_P1_ADD, OPR_P1_SUBTRACT)
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
    subroutine OPR_Partial_Y_MPISplit(type, nx, ny, nz, u, result, tmp1)
        use TLabMPI_PROCS, only: TLabMPI_Halos_Y
        use TLab_Pointers_2D, only: pxz_wrk3d
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation

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

        np1 = size(fdm_der1_Y_split%rhs, 2)/2
        np2 = size(fdm_der2_Y_split%rhs, 2)/2
        np = max(np1, np2)
        call TLabMPI_Halos_Y(result, nx*nz, np, pxz_halo_m(:, 1), pxz_halo_p(:, 1))

        select case (type)
        case (OPR_P2)
            call fdm_der2_Y_split%compute(nx*nz, result, pxz_halo_m(:, np - np2 + 1:np), pxz_halo_p, wrk3d)

        case (OPR_P2_P1)
            call fdm_der2_Y_split%compute(nx*nz, result, pxz_halo_m(:, np - np2 + 1:np), pxz_halo_p, tmp1)
            call fdm_der1_Y_split%compute(nx*nz, result, pxz_halo_m(:, np - np1 + 1:np), pxz_halo_p, wrk3d)

        case (OPR_P1, OPR_P1_ADD, OPR_P1_SUBTRACT)
            call fdm_der1_Y_split%compute(nx*nz, result, pxz_halo_m(:, np - np1 + 1:np), pxz_halo_p, wrk3d)

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
    subroutine OPR_Partial_Z(type, nx, ny, nz, u, result, tmp1)
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation

        ! ###################################################################
        if (z%size == 1) then
            result = 0.0_wp
            if (type == OPR_P2_P1) tmp1 = 0.0_wp
            return
        end if

        select case (type)
        case (OPR_P2)
            if (.not. z%uniform) call fdm_der1_Z%compute(nx*ny, u, wrk3d)
            call fdm_der2_Z%compute(nx*ny, u, result, wrk3d)

        case (OPR_P2_P1)
            call fdm_der1_Z%compute(nx*ny, u, tmp1)
            call fdm_der2_Z%compute(nx*ny, u, result, tmp1)

        case (OPR_P1)
            call fdm_der1_Z%compute(nx*ny, u, result)

        end select

        return
    end subroutine OPR_Partial_Z

    !########################################################################
    !########################################################################
    ! Imposing zero derivative at the boundaries
    subroutine OPR_Partial_Z_Bcs(nx, ny, nz, u, result, ibc)
        use TLab_Constants, only: BCS_DD, BCS_ND, BCS_DN, BCS_NN
        use TLab_Arrays, only: wrk2d
        use FDM_Derivative_1order, only: der1_biased
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        integer, intent(in) :: ibc

#define bcs_b(i) wrk2d(i,1)
#define bcs_t(i) wrk2d(i,2)

        select type (fdm_der1_Z)
        type is (der1_biased)
            select case (ibc)
            case (BCS_DD)
                call fdm_der1_Z%bcsDD%compute(nx*ny, u, result)

            case (BCS_ND)
                bcs_b(1:nx*ny) = 0.0_wp
                call fdm_der1_Z%bcsND%compute(nx*ny, u, result, bcs_b(:))

            case (BCS_DN)
                bcs_t(1:nx*ny) = 0.0_wp
                call fdm_der1_Z%bcsDN%compute(nx*ny, u, result, bcs_t(:))

            case (BCS_NN)
                bcs_b(1:nx*ny) = 0.0_wp
                bcs_t(1:nx*ny) = 0.0_wp
                call fdm_der1_Z%bcsNN%compute(nx*ny, u, result, bcs_b(:), bcs_t(:))

            end select
        end select

#undef bcs_b
#undef bcs_t

        return
    end subroutine OPR_Partial_Z_Bcs

end module OPR_Partial
