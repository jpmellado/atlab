#include "tlab_error.h"

module TLabMPI_PROCS
    use mpi_f08
    use TLab_Constants, only: wp, dp, sp, wi, lfile, efile
    use TLab_Memory, only: imax, jmax, kmax
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    use TLabMPI_VARS
    implicit none
    private

    public :: TLabMPI_Initialize
    public :: TLabMPI_Halos_X, TLabMPI_Halos_Y
    public :: TLabMPI_Panic

contains
    ! ######################################################################
    ! ######################################################################
    subroutine TLabMPI_Initialize(inifile)
        use TLab_Grid, only: xSubgrid, ySubgrid, zSubgrid
        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        integer(wi) dims(2), coord(2)
        logical period(2), remain_dims(2), reorder
        type(MPI_Comm) :: ims_comm_xy                               ! Plane communicators

        character(len=512) line
        character*64 lstr

        ! #######################################################################
        xMpi%num_processors = xSubgrid%parent%size/xSubgrid%size
        yMpi%num_processors = ySubgrid%parent%size/ySubgrid%size
        zMpi%num_processors = zSubgrid%parent%size/zSubgrid%size

        ! consistency check
        if (xMpi%num_processors*yMpi%num_processors == ims_npro) then
            line = 'Initializing MPI domain partition'
            write (lstr, *) xMpi%num_processors
            line = trim(adjustl(line))//' '//trim(adjustl(lstr))
            write (lstr, *) yMpi%num_processors
            line = trim(adjustl(line))//'x'//trim(adjustl(lstr))
            call TLab_Write_ASCII(lfile, trim(adjustl(line)))
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Inconsistency in total number of PEs')
            call TLab_Stop(DNS_ERROR_DIMGRID)
        end if

        ! #######################################################################
        call TLab_Write_ASCII(lfile, 'Creating MPI communicators.')

        ! the first index in the grid corresponds to j, the second to i
        ! we tried transposing the mpi axes, but similar time
        dims(1) = yMpi%num_processors
        dims(2) = xMpi%num_processors
        period = .true.; reorder = .false.
        call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dims, period, reorder, ims_comm_xy, ims_err)

        remain_dims(1) = .false.; remain_dims(2) = .true.
        call MPI_CART_SUB(ims_comm_xy, remain_dims, xMpi%comm, ims_err)

        remain_dims(1) = .true.; remain_dims(2) = .false.
        call MPI_CART_SUB(ims_comm_xy, remain_dims, yMpi%comm, ims_err)

        call MPI_CART_COORDS(ims_comm_xy, ims_pro, 2, coord, ims_err)
        xMpi%rank = coord(2)
        yMpi%rank = coord(1)

        ! to be removed
        ims_npro_i = xMpi%num_processors !xSubgrid%parent%size/xSubgrid%size
        ims_npro_j = yMpi%num_processors !ySubgrid%parent%size/ySubgrid%size
        ims_npro_k = zMpi%num_processors !1
        !
        ims_comm_x = xMpi%comm
        ims_comm_y = yMpi%comm
        ims_pro_i = xMpi%rank
        ims_pro_j = yMpi%rank

        ! #######################################################################
        ! local offset in grid points##############
        xSubgrid%offset = xSubgrid%size*xMpi%rank
        ySubgrid%offset = ySubgrid%size*yMpi%rank

        ! #######################################################################
        ! Control of MPI type
        select case (wp)
        case (dp)
            TLAB_MPI_REAL_TYPE = MPI_REAL8
        case (sp)
            TLAB_MPI_REAL_TYPE = MPI_REAL4
        end select

        return
    end subroutine TLabMPI_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine TLabMPI_Panic(location, mpi_error_code)
        character(len=*), intent(in) :: location
        integer, intent(in) :: mpi_error_code

        !##############################
        character error_string*1024
        integer error_local, error_len

        call MPI_Error_String(mpi_error_code, error_string, error_len, error_local)
        call TLab_Write_ASCII(efile, 'MPI-ERROR: Source file'//trim(adjustl(LOCATION)), .true.)
        call TLab_Write_ASCII(efile, error_string, .true.)

        call TLab_Stop(mpi_error_code)
        ! Not supposed to return from this subroutine

    end subroutine TLabMPI_Panic

    ! ###################################################################
    ! ###################################################################
    subroutine TLabMPI_Halos_X(a, size_plane, n_halo_planes, halo_m, halo_p)
        real(wp), intent(in) :: a(:)
        integer(wi), intent(in) :: size_plane
        integer(wi), intent(in) :: n_halo_planes
        real(wp), intent(out) :: halo_m(:)      ! minus, coming from left/west processor
        real(wp), intent(out) :: halo_p(:)      ! plus, coming from right/east processor

        integer(wi) :: counts, disp
        integer source, dest

        ! ###################################################################
        counts = size_plane*n_halo_planes

        ! pass to previous processor
        dest = mod(ims_pro_i - 1 + ims_npro_i, ims_npro_i)
        source = mod(ims_pro_i + 1, ims_npro_i)
        disp = 1
        call MPI_Sendrecv(a(disp), counts, MPI_REAL8, dest, 0, &
                          halo_p, counts, MPI_REAL8, source, 0, &
                          ims_comm_x, MPI_STATUS_IGNORE, ims_err)

        ! pass to following processor
        dest = mod(ims_pro_i + 1, ims_npro_i)
        source = mod(ims_pro_i - 1 + ims_npro_i, ims_npro_i)
        disp = size(a) - counts + 1
        call MPI_Sendrecv(a(disp), counts, MPI_REAL8, dest, 1, &
                          halo_m, counts, MPI_REAL8, source, 1, &
                          ims_comm_x, MPI_STATUS_IGNORE, ims_err)

        return
    end subroutine TLabMPI_Halos_X

    ! ###################################################################
    ! ###################################################################
    subroutine TLabMPI_Halos_Y(a, size_plane, n_halo_planes, halo_m, halo_p)
        real(wp), intent(in) :: a(:)
        integer(wi), intent(in) :: size_plane
        integer(wi), intent(in) :: n_halo_planes
        real(wp), intent(out) :: halo_m(:)      ! minus, coming from left/west processor
        real(wp), intent(out) :: halo_p(:)      ! plus, coming from right/east processor

        integer(wi) :: counts, disp
        integer source, dest

        ! ###################################################################
        counts = size_plane*n_halo_planes

        ! pass to previous processor
        dest = mod(ims_pro_j - 1 + ims_npro_j, ims_npro_j)
        source = mod(ims_pro_j + 1, ims_npro_j)
        disp = 1
        call MPI_Sendrecv(a(disp), counts, MPI_REAL8, dest, 0, &
                          halo_p, counts, MPI_REAL8, source, 0, &
                          ims_comm_y, MPI_STATUS_IGNORE, ims_err)

        ! pass to following processor
        dest = mod(ims_pro_j + 1, ims_npro_j)
        source = mod(ims_pro_j - 1 + ims_npro_j, ims_npro_j)
        disp = size(a) - counts + 1
        call MPI_Sendrecv(a(disp), counts, MPI_REAL8, dest, 1, &
                          halo_m, counts, MPI_REAL8, source, 1, &
                          ims_comm_y, MPI_STATUS_IGNORE, ims_err)

        return
    end subroutine TLabMPI_Halos_Y

end module TLabMPI_PROCS
