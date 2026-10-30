#include "tlab_error.h"

module TLabMPI_Grid
    use mpi_f08
    use TLab_Constants, only: wp, wi
    use TLab_Constants, only: efile
    use TLab_WorkFlow, only: num_processors, rank
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Log
    implicit none
    private

    public :: TLabMPI_Grid_Initialize
    public :: TLabMPI_Halos_X, TLabMPI_Halos_Y

    type, public :: mpi_axis_dt
        type(MPI_Comm) :: comm
        integer :: num_processors       ! we could name it size, in analogy to spatial grid, but this might be clearer
        integer :: rank
    end type

    type, extends(mpi_axis_dt) :: mpi_grid_dt
        type(mpi_axis_dt) axes(3)
    end type

    type(mpi_grid_dt), target, public :: mpiGrid
    type(mpi_axis_dt), pointer, public :: xMpi => null()    ! For readability in the code
    type(mpi_axis_dt), pointer, public :: yMpi => null()
    type(mpi_axis_dt), pointer, public :: zMpi => null()

    ! -----------------------------------------------------------------------
    integer ierr
    
contains
    ! ######################################################################
    ! ######################################################################
    subroutine TLabMPI_Grid_Initialize(inifile)
        use TLab_Grid, only: xSubgrid, ySubgrid, zSubgrid, subgrid
        character(len=*), intent(in) :: inifile

        ! -----------------------------------------------------------------------
        integer(wi) dims(2), coord(2)
        logical period(2), remain_dims(2), reorder
        type(MPI_Comm) :: mpi_comm_xy                             ! Plane communicators

        character(len=512) line
        character*64 lstr
        integer ig

        ! #######################################################################
        mpiGrid%comm = MPI_COMM_WORLD               ! Pass global information to mpiGrid structure
        mpiGrid%num_processors = num_processors
        mpiGrid%rank = rank

        do ig = 1, 3
            mpiGrid%axes(ig)%num_processors = subgrid%axes(ig)%parent%size/subgrid%axes(ig)%size
        end do

        xMpi => mpiGrid%axes(1)     ! for readability
        yMpi => mpiGrid%axes(2)
        zMpi => mpiGrid%axes(3)

        ! consistency check
        if (xMpi%num_processors*yMpi%num_processors == mpiGrid%num_processors) then
            line = 'Initializing MPI domain partition'
            write (lstr, *) xMpi%num_processors
            line = trim(adjustl(line))//' '//trim(adjustl(lstr))
            write (lstr, *) yMpi%num_processors
            line = trim(adjustl(line))//'x'//trim(adjustl(lstr))
            call TLab_Log(trim(adjustl(line)))
        else
            call TLab_Write_ASCII(efile, __FILE__//'. Inconsistency in total number of PEs')
            call TLab_Stop(DNS_ERROR_DIMGRID)
        end if

        ! #######################################################################
        call TLab_Log('Creating MPI communicators.')

        ! the first index in the grid corresponds to j, the second to i
        ! we tried transposing the mpi axes, but similar time
        dims(1) = yMpi%num_processors
        dims(2) = xMpi%num_processors
        period = .true.; reorder = .false.
        call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dims, period, reorder, mpi_comm_xy, ierr)

        remain_dims(1) = .false.; remain_dims(2) = .true.
        call MPI_CART_SUB(mpi_comm_xy, remain_dims, xMpi%comm, ierr)

        remain_dims(1) = .true.; remain_dims(2) = .false.
        call MPI_CART_SUB(mpi_comm_xy, remain_dims, yMpi%comm, ierr)

        call MPI_CART_COORDS(mpi_comm_xy, mpiGrid%rank, 2, coord, ierr)
        xMpi%rank = coord(2)
        yMpi%rank = coord(1)

        ! #######################################################################
        ! local offset in grid points
        xSubgrid%offset = xSubgrid%size*xMpi%rank
        ySubgrid%offset = ySubgrid%size*yMpi%rank

        return
    end subroutine TLabMPI_Grid_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine TLabMPI_Halos_X(a, size_plane, n_halo_planes, halo_m, halo_p)
        real(wp), intent(in) :: a(:)
        integer(wi), intent(in) :: size_plane
        integer(wi), intent(in) :: n_halo_planes
        real(wp), intent(out) :: halo_m(*)      ! minus, coming from left/west processor
        real(wp), intent(out) :: halo_p(*)      ! plus, coming from right/east processor

        integer(wi) :: counts, disp
        integer source, dest

        ! ###################################################################
        counts = size_plane*n_halo_planes

        ! pass to previous processor
        dest = mod(xMpi%rank - 1 + xMpi%num_processors, xMpi%num_processors)
        source = mod(xMpi%rank + 1, xMpi%num_processors)
        disp = 1
        call MPI_Sendrecv(a(disp), counts, MPI_REAL8, dest, 0, &
                          halo_p, counts, MPI_REAL8, source, 0, &
                          xMpi%comm, MPI_STATUS_IGNORE, ierr)

        ! pass to following processor
        dest = mod(xMpi%rank + 1, xMpi%num_processors)
        source = mod(xMpi%rank - 1 + xMpi%num_processors, xMpi%num_processors)
        disp = size(a) - counts + 1
        call MPI_Sendrecv(a(disp), counts, MPI_REAL8, dest, 1, &
                          halo_m, counts, MPI_REAL8, source, 1, &
                          xMpi%comm, MPI_STATUS_IGNORE, ierr)

        return
    end subroutine TLabMPI_Halos_X

    ! ###################################################################
    ! ###################################################################
    subroutine TLabMPI_Halos_Y(a, size_plane, n_halo_planes, halo_m, halo_p)
        real(wp), intent(in) :: a(:)
        integer(wi), intent(in) :: size_plane
        integer(wi), intent(in) :: n_halo_planes
        real(wp), intent(out) :: halo_m(*)      ! minus, coming from left/west processor
        real(wp), intent(out) :: halo_p(*)      ! plus, coming from right/east processor

        integer(wi) :: counts, disp
        integer source, dest

        ! ###################################################################
        counts = size_plane*n_halo_planes

        ! pass to previous processor
        dest = mod(yMpi%rank - 1 + yMpi%num_processors, yMpi%num_processors)
        source = mod(yMpi%rank + 1, yMpi%num_processors)
        disp = 1
        call MPI_Sendrecv(a(disp), counts, MPI_REAL8, dest, 0, &
                          halo_p, counts, MPI_REAL8, source, 0, &
                          yMpi%comm, MPI_STATUS_IGNORE, ierr)

        ! pass to following processor
        dest = mod(yMpi%rank + 1, yMpi%num_processors)
        source = mod(yMpi%rank - 1 + yMpi%num_processors, yMpi%num_processors)
        disp = size(a) - counts + 1
        call MPI_Sendrecv(a(disp), counts, MPI_REAL8, dest, 1, &
                          halo_m, counts, MPI_REAL8, source, 1, &
                          yMpi%comm, MPI_STATUS_IGNORE, ierr)

        return
    end subroutine TLabMPI_Halos_Y

end module TLabMPI_Grid
