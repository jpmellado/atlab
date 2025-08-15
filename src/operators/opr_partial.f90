#include "tlab_error.h"

module OPR_Partial
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: BCS_NONE
    use TLab_Arrays, only: wrk2d, wrk3d
#ifdef USE_MPI
    use TLabMPI_VARS, only: ims_npro_i, ims_npro_j
    use TLabMPI_Transpose
    use FDM_Derivative_MPISplit
#endif
    use FDM, only: g
    use FDM_Derivative, only: FDM_Der1_Solve, FDM_Der2_Solve
    use Thomas3_Split
    implicit none
    private

    public :: OPR_Partial_Initialize
    public :: OPR_Partial_X     ! These first 2 could be written in terms of the last one...
    public :: OPR_Partial_Y
    public :: OPR_Partial_Z

    integer, parameter, public :: OPR_P1 = 1                ! 1. order derivative
    integer, parameter, public :: OPR_P2 = 2                ! 2. order derivative
    integer, parameter, public :: OPR_P2_P1 = 3             ! 2. and 1.order derivatives

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
    real(wp), pointer :: pyz_halo_m(:, :) => null(), pyz_halo_p(:, :) => null()
    real(wp), pointer :: pxz_halo_m(:, :) => null(), pxz_halo_p(:, :) => null()
    real(wp), allocatable :: wrk_split(:)
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
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        integer der_mode_i, der_mode_j
        integer, parameter :: TYPE_TRANSPOSE = 1
        integer, parameter :: TYPE_SPLIT = 2

#ifdef USE_MPI
        integer np, idummy
#endif

        ! #######################################################################
        ! Read data
        bakfile = trim(adjustl(inifile))//'.bak'

        block = 'Parallel'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call ScanFile_Char(bakfile, inifile, block, 'DerivativeModeI', 'transpose', sRes)
        if (trim(adjustl(sRes)) == 'transpose') then; der_mode_i = TYPE_TRANSPOSE
        elseif (trim(adjustl(sRes)) == 'split') then; der_mode_i = TYPE_SPLIT
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong DerivativeModeI option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        call ScanFile_Char(bakfile, inifile, block, 'DerivativeModeJ', 'transpose', sRes)
        if (trim(adjustl(sRes)) == 'transpose') then; der_mode_j = TYPE_TRANSPOSE
        elseif (trim(adjustl(sRes)) == 'split') then; der_mode_j = TYPE_SPLIT
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong DerivativeModeJ option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

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

            ! allocate (wrk_split(max(imax*kmax*(ims_npro_j + 1), jmax*kmax*(ims_npro_i + 1))))
            idummy = max(imax*kmax*(ims_npro_j + 1), jmax*kmax*(ims_npro_i + 1))
            call TLab_Allocate_Real(__FILE__, wrk_split, [idummy], 'wrk-split')
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
        if (g(1)%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (present(tmp1)) tmp1 = 0.0_wp
            return
        end if

        ! Transposition: make x-direction the last one
#ifdef USE_ESSL
        call DGETMO(u, g(1)%size, g(1)%size, ny*nz, result, ny*nz)
#else
        call TLab_Transpose(u, g(1)%size, ny*nz, g(1)%size, result, ny*nz)
#endif

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        select case (type)
        case (OPR_P2)
            if (g(1)%der2%need_1der) call FDM_Der1_Solve(ny*nz, ibc_loc, g(1)%der1, g(1)%der1%lu, result, tmp1, wrk2d)
            call FDM_Der2_Solve(ny*nz, g(1)%der2, g(1)%der2%lu, result, wrk3d, tmp1, wrk2d)

        case (OPR_P2_P1)
            call FDM_Der1_Solve(ny*nz, ibc_loc, g(1)%der1, g(1)%der1%lu, result, wrk3d, wrk2d)
            call FDM_Der2_Solve(ny*nz, g(1)%der2, g(1)%der2%lu, result, tmp1, wrk3d, wrk2d)

        case (OPR_P1)
            call FDM_Der1_Solve(ny*nz, ibc_loc, g(1)%der1, g(1)%der1%lu, result, wrk3d, wrk2d)

        end select

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        if (type == OPR_P2_P1) then
            call DGETMO(tmp1, ny*nz, ny*nz, g(1)%size, result, g(1)%size)
            call DGETMO(wrk3d, ny*nz, ny*nz, g(1)%size, tmp1, g(1)%size)
        else
            call DGETMO(wrk3d, ny*nz, ny*nz, g(1)%size, result, g(1)%size)
        end if
#else
        if (type == OPR_P2_P1) then
            call TLab_Transpose(tmp1, ny*nz, g(1)%size, ny*nz, result, g(1)%size)
            call TLab_Transpose(wrk3d, ny*nz, g(1)%size, ny*nz, tmp1, g(1)%size)
        else
            call TLab_Transpose(wrk3d, ny*nz, g(1)%size, ny*nz, result, g(1)%size)
        end if
#endif

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
        if (g(1)%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (present(tmp1)) tmp1 = 0.0_wp
            return
        end if

        nlines = tmpi_plan_dx%nlines

        ! Transposition: make x-direction the last one
        call TLabMPI_Trp_ExecI_Forward(u, result, tmpi_plan_dx)
#ifdef USE_ESSL
        call DGETMO(result, g(1)%size, g(1)%size, nlines, wrk3d, nlines)
#else
        call TLab_Transpose(result, g(1)%size, nlines, g(1)%size, wrk3d, nlines)
#endif

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        select case (type)
        case (OPR_P2)
            if (g(1)%der2%need_1der) call FDM_Der1_Solve(nlines, ibc_loc, g(1)%der1, g(1)%der1%lu, wrk3d, tmp1, wrk2d)
            call FDM_Der2_Solve(nlines, g(1)%der2, g(1)%der2%lu, wrk3d, result, tmp1, wrk2d)

        case (OPR_P2_P1)
            call FDM_Der1_Solve(nlines, ibc_loc, g(1)%der1, g(1)%der1%lu, wrk3d, tmp1, wrk2d)
            call FDM_Der2_Solve(nlines, g(1)%der2, g(1)%der2%lu, wrk3d, result, tmp1, wrk2d)

        case (OPR_P1)
            call FDM_Der1_Solve(nlines, ibc_loc, g(1)%der1, g(1)%der1%lu, wrk3d, result, wrk2d)

        end select

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        if (type == OPR_P2_P1) then
            call DGETMO(tmp1, nlines, nlines, g(1)%size, wrk3d, g(1)%size)
            call DGETMO(result, nlines, nlines, g(1)%size, tmp1, g(1)%size)
            call TLabMPI_Trp_ExecI_Backward(tmp1, result, tmpi_plan_dx)
            call TLabMPI_Trp_ExecI_Backward(wrk3d, tmp1, tmpi_plan_dx)
        else
            call DGETMO(result, nlines, nlines, g(1)%size, wrk3d, g(1)%size)
            call TLabMPI_Trp_ExecI_Backward(wrk3d, result, tmpi_plan_dx)
        end if
#else
        if (type == OPR_P2_P1) then
            call TLab_Transpose(tmp1, nlines, g(1)%size, nlines, wrk3d, g(1)%size)
            call TLab_Transpose(result, nlines, g(1)%size, nlines, tmp1, g(1)%size)
            call TLabMPI_Trp_ExecI_Backward(tmp1, result, tmpi_plan_dx)
            call TLabMPI_Trp_ExecI_Backward(wrk3d, tmp1, tmpi_plan_dx)
        else
            call TLab_Transpose(result, nlines, g(1)%size, nlines, wrk3d, g(1)%size)
            call TLabMPI_Trp_ExecI_Backward(wrk3d, result, tmpi_plan_dx)
        end if
#endif

        return
    end subroutine OPR_Partial_X_MPITranspose

    !########################################################################
    !########################################################################
    subroutine OPR_Partial_X_MPISplit(type, nx, ny, nz, u, result, tmp1, ibc)
        use TLabMPI_PROCS, only: TLabMPI_Halos_X
        use TLab_Pointers_2D, only: pyz_wrk3d
        integer(wi), intent(in) :: type                         ! OPR_P1, OPR_P2, OPR_P2_P1
        integer(wi), intent(in) :: nx, ny, nz
        real(wp), intent(in) :: u(nx*ny*nz)
        real(wp), intent(out) :: result(nx*ny*nz)
        real(wp), intent(inout), optional :: tmp1(nx*ny*nz)     ! 1. order derivative in 2. order calculation
        integer, intent(in), optional :: ibc                    ! boundary conditions 1. order derivative

        ! -------------------------------------------------------------------
        integer np, np1, np2

        ! ###################################################################
        if (g(1)%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (present(tmp1)) tmp1 = 0.0_wp
            return
        end if

        ! Transposition: make x-direction the last one
#ifdef USE_ESSL
        call DGETMO(u, nx, nx, ny*nz, result, ny*nz)
#else
        call TLab_Transpose(u, nx, ny*nz, nx, result, ny*nz)
#endif

        ! -------------------------------------------------------------------
        np1 = size(der1_split_x%rhs)/2
        np2 = size(der2_split_x%rhs)/2
        np = max(np1, np2)
        call TLabMPI_Halos_X(result, ny*nz, np, pyz_halo_m(:, 1), pyz_halo_p(:, 1))

        select case (type)
        case (OPR_P2)
            call FDM_MPISplit_Solve(ny*nz, nx, der2_split_x, result, &
                                    pyz_halo_m(:, np - np2 + 1:np), pyz_halo_p, pyz_wrk3d, wrk_split)

        case (OPR_P2_P1)
            call FDM_MPISplit_Solve(ny*nz, nx, der2_split_x, result, &
                                    pyz_halo_m(:, np - np2 + 1:np), pyz_halo_p, tmp1, wrk_split)
            call FDM_MPISplit_Solve(ny*nz, nx, der1_split_x, result, &
                                    pyz_halo_m(:, np - np1 + 1:np), pyz_halo_p, pyz_wrk3d, wrk_split)

        case (OPR_P1)
            call FDM_MPISplit_Solve(ny*nz, nx, der1_split_x, result, &
                                    pyz_halo_m(:, np - np1 + 1:np), pyz_halo_p, pyz_wrk3d, wrk_split)

        end select

        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        if (type == OPR_P2_P1) then
            call DGETMO(tmp1, ny*nz, ny*nz, nx, result, nx)
            call DGETMO(wrk3d, ny*nz, ny*nz, nx, tmp1, nx)
        else
            call DGETMO(wrk3d, ny*nz, ny*nz, nx, result, nx)
        end if
#else
        if (type == OPR_P2_P1) then
            call TLab_Transpose(tmp1, ny*nz, nx, ny*nz, result, nx)
            call TLab_Transpose(wrk3d, ny*nz, nx, ny*nz, tmp1, nx)
        else
            call TLab_Transpose(wrk3d, ny*nz, nx, ny*nz, result, nx)
        end if
#endif

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
        if (g(2)%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (present(tmp1)) tmp1 = 0.0_wp
            return
        end if

        ! Transposition: make y-direction the last one
#ifdef USE_ESSL
        call DGETMO(u, nx*ny, nx*ny, nz, result, nz)
#else
        call TLab_Transpose(u, nx*ny, nz, nx*ny, result, nz)
#endif
        nlines = nx*nz

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        select case (type)
        case (OPR_P2)
            if (g(2)%der2%need_1der) call FDM_Der1_Solve(nlines, ibc_loc, g(2)%der1, g(2)%der1%lu, result, tmp1, wrk2d)
            call FDM_Der2_Solve(nlines, g(2)%der2, g(2)%der2%lu, result, wrk3d, tmp1, wrk2d)

        case (OPR_P2_P1)
            call FDM_Der1_Solve(nlines, ibc_loc, g(2)%der1, g(2)%der1%lu, result, wrk3d, wrk2d)
            call FDM_Der2_Solve(nlines, g(2)%der2, g(2)%der2%lu, result, tmp1, wrk3d, wrk2d)

        case (OPR_P1)
            call FDM_Der1_Solve(nlines, ibc_loc, g(2)%der1, g(2)%der1%lu, result, wrk3d, wrk2d)

        end select

        ! ###################################################################
        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        if (type == OPR_P2_P1) then
            call DGETMO(tmp1, nz, nz, nx*ny, result, nx*ny)
            call DGETMO(wrk3d, nz, nz, nx*ny, tmp1, nx*ny)
        else
            call DGETMO(wrk3d, nz, nz, nx*ny, result, nx*ny)
        end if
#else
        if (type == OPR_P2_P1) then
            call TLab_Transpose(tmp1, nz, nx*ny, nz, result, nx*ny)
            call TLab_Transpose(wrk3d, nz, nx*ny, nz, tmp1, nx*ny)
        else
            call TLab_Transpose(wrk3d, nz, nx*ny, nz, result, nx*ny)
        end if
#endif

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
        if (g(2)%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (present(tmp1)) tmp1 = 0.0_wp
            return
        end if

        ! Transposition: make y-direction the last one
#ifdef USE_ESSL
        call DGETMO(u, nx*ny, nx*ny, nz, result, nz)
#else
        call TLab_Transpose(u, nx*ny, nz, nx*ny, result, nz)
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
            if (g(2)%der2%need_1der) call FDM_Der1_Solve(nlines, ibc_loc, g(2)%der1, g(2)%der1%lu, wrk3d, tmp1, wrk2d)
            call FDM_Der2_Solve(nlines, g(2)%der2, g(2)%der2%lu, wrk3d, result, tmp1, wrk2d)

        case (OPR_P2_P1)
            call FDM_Der1_Solve(nlines, ibc_loc, g(2)%der1, g(2)%der1%lu, wrk3d, tmp1, wrk2d)
            call FDM_Der2_Solve(nlines, g(2)%der2, g(2)%der2%lu, wrk3d, result, tmp1, wrk2d)

        case (OPR_P1)
            call FDM_Der1_Solve(nlines, ibc_loc, g(2)%der1, g(2)%der1%lu, wrk3d, result, wrk2d)

        end select

        ! ###################################################################
        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        if (type == OPR_P2_P1) then
            call TLabMPI_Trp_ExecJ_Backward(tmp1, wrk3d, tmpi_plan_dy)
            call TLabMPI_Trp_ExecJ_Backward(result, tmp1, tmpi_plan_dy)
            call DGETMO(tmp1, nz, nz, nx*ny, result, nx*ny)
            call DGETMO(wrk3d, nz, nz, nx*ny, tmp1, nx*ny)
        else
            call TLabMPI_Trp_ExecJ_Backward(result, wrk3d, tmpi_plan_dy)
            call DGETMO(wrk3d, nz, nz, nx*ny, result, nx*ny)
        end if
#else
        if (type == OPR_P2_P1) then
            call TLabMPI_Trp_ExecJ_Backward(tmp1, wrk3d, tmpi_plan_dy)
            call TLabMPI_Trp_ExecJ_Backward(result, tmp1, tmpi_plan_dy)
            call TLab_Transpose(tmp1, nz, nx*ny, nz, result, nx*ny)
            call TLab_Transpose(wrk3d, nz, nx*ny, nz, tmp1, nx*ny)
        else
            call TLabMPI_Trp_ExecJ_Backward(result, wrk3d, tmpi_plan_dy)
            call TLab_Transpose(wrk3d, nz, nx*ny, nz, result, nx*ny)
        end if
#endif

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
        if (g(2)%size == 1) then ! Set to zero in 2D case
            result = 0.0_wp
            if (present(tmp1)) tmp1 = 0.0_wp
            return
        end if

        ! Transposition: make y-direction the last one
#ifdef USE_ESSL
        call DGETMO(u, nx*ny, nx*ny, nz, result, nz)
#else
        call TLab_Transpose(u, nx*ny, nz, nx*ny, result, nz)
#endif

        ! -------------------------------------------------------------------
        np1 = size(der1_split_y%rhs)/2
        np2 = size(der2_split_y%rhs)/2
        np = max(np1, np2)
        call TLabMPI_Halos_Y(result, nx*nz, np, pxz_halo_m(:, 1), pxz_halo_p(:, 1))

        select case (type)
        case (OPR_P2)
            call FDM_MPISplit_Solve(nx*nz, ny, der2_split_y, result, &
                                    pxz_halo_m(:, np - np2 + 1:np), pxz_halo_p, pxz_wrk3d, wrk_split)

        case (OPR_P2_P1)
            call FDM_MPISplit_Solve(nx*nz, ny, der2_split_y, result, &
                                    pxz_halo_m(:, np - np2 + 1:np), pxz_halo_p, tmp1, wrk_split)
            call FDM_MPISplit_Solve(nx*nz, ny, der1_split_y, result, &
                                    pxz_halo_m(:, np - np1 + 1:np), pxz_halo_p, pxz_wrk3d, wrk_split)

        case (OPR_P1)
            call FDM_MPISplit_Solve(nx*nz, ny, der1_split_y, result, &
                                    pxz_halo_m(:, np - np1 + 1:np), pxz_halo_p, pxz_wrk3d, wrk_split)

        end select

        ! ###################################################################
        ! Put arrays back in the order in which they came in
#ifdef USE_ESSL
        if (type == OPR_P2_P1) then
            call DGETMO(tmp1, nz, nz, nx*ny, result, nx*ny)
            call DGETMO(wrk3d, nz, nz, nx*ny, tmp1, nx*ny)
        else
            call DGETMO(wrk3d, nz, nz, nx*ny, result, nx*ny)
        end if
#else
        if (type == OPR_P2_P1) then
            call TLab_Transpose(tmp1, nz, nx*ny, nz, result, nx*ny)
            call TLab_Transpose(wrk3d, nz, nx*ny, nz, tmp1, nx*ny)
        else
            call TLab_Transpose(wrk3d, nz, nx*ny, nz, result, nx*ny)
        end if
#endif

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
        if (g(3)%size == 1) then
            result = 0.0_wp
            if (present(tmp1)) tmp1 = 0.0_wp
            return
        end if

        if (present(ibc)) then
            ibc_loc = ibc
        else
            ibc_loc = BCS_NONE
        end if

        select case (type)
        case (OPR_P2)
            if (g(3)%der2%need_1der) call FDM_Der1_Solve(nx*ny, ibc_loc, g(3)%der1, g(3)%der1%lu, u, wrk3d, wrk2d)
            call FDM_Der2_Solve(nx*ny, g(3)%der2, g(3)%der2%lu, u, result, wrk3d, wrk2d)

        case (OPR_P2_P1)
            call FDM_Der1_Solve(nx*ny, ibc_loc, g(3)%der1, g(3)%der1%lu, u, tmp1, wrk2d)
            call FDM_Der2_Solve(nx*ny, g(3)%der2, g(3)%der2%lu, u, result, tmp1, wrk2d)

        case (OPR_P1)
            call FDM_Der1_Solve(nx*ny, ibc_loc, g(3)%der1, g(3)%der1%lu, u, result, wrk2d)

        end select

        return
    end subroutine OPR_Partial_Z

end module OPR_Partial
