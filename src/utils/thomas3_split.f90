#include "tlab_error.h"

! Matrix splitting

module Thomas3_Split
    use TLab_Constants, only: wp, wi, small_wp, efile
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
#ifdef USE_MPI
    use mpi_f08, only: MPI_Comm
#endif
    use Thomas3
    implicit none
    private

    public :: Thomas3_Split_Initialize
    public :: Thomas3_Split_Solve_Serial
    public :: Thomas3_Split_Solve_Serial2
#ifdef USE_MPI
    public :: Thomas3_Split_Solve_MPI
    public :: Thomas3_Split_Solve_MPI2
#endif

    type, public :: thomas3_split_dt
#ifdef USE_MPI
        type(MPI_Comm) communicator
        integer rank, n_ranks
#endif
        integer :: block_id
        integer(wi) :: nmin, nmax
        logical :: circulant = .true.
        real(wp), allocatable :: lhs(:, :)
        real(wp), allocatable :: z(:, :)
        real(wp), allocatable :: y(:, :)
        real(wp) :: alpha(2) = [0.0_wp, 0.0_wp]
        real(wp) :: beta = 0.0_wp
        real(wp) :: gamma = 0.0_wp
    end type thomas3_split_dt

    type, public :: data_dt                     ! for testing algorithm in serial
        real(wp), pointer :: p(:, :)
    end type data_dt

contains
    !########################################################################
    !########################################################################
    subroutine Thomas3_Split_Initialize(a, b, c, points, split)
        real(wp), intent(inout) :: a(:), b(:), c(:)
        integer(wi), intent(in) :: points(:)   ! sequence of splitting points in ascending order
        type(thomas3_split_dt), intent(inout) :: split

        ! -------------------------------------------------------------------
        integer nblocks, m, k
        integer(wi) nsize
        integer(wi) p, p_plus_1

        real(wp), allocatable :: aloc(:), bloc(:), cloc(:), zloc(:)
        real(wp) alpha_0(2), alpha(2)
        real(wp) delta
        real(wp) alpha_previous(2), beta_previous, gamma_loc

        !########################################################################
        nblocks = size(points)
        ! Number of coefficients are nblocks-1 for the tridiagonal case
        ! and we add one for the circulant case, which is managed by last block

        nsize = size(a)

        k = split%block_id
        split%nmin = points(mod(k - 2 + nblocks, nblocks) + 1)
        split%nmin = mod(split%nmin, nsize) + 1
        split%nmax = points(k)

        ! -------------------------------------------------------------------
        ! memory management
        nsize = split%nmax - split%nmin + 1

        if (allocated(split%lhs)) deallocate (split%lhs)
        allocate (split%lhs(1:nsize, 3))
        split%lhs(:, :) = 0.0_wp

        if (allocated(split%z)) deallocate (split%z)
        allocate (split%z(1:nsize, nblocks))
        split%z(:, :) = 0.0_wp

        if (allocated(split%y)) deallocate (split%y)
        allocate (split%y(1:nsize, nblocks))
        split%y(:, :) = 0.0_wp

        ! -------------------------------------------------------------------
        ! temporary arrays to calculate z_j
        nsize = size(a)
        allocate (aloc(nsize), bloc(nsize), cloc(nsize), zloc(nsize))

        if (split%circulant) then
            m = nblocks             ! last block has the circulant coefficient
            p = points(m)           ! index of splitting point
            p_plus_1 = 1

            call Splitting(a, b, c, p, p_plus_1, alpha_0, zloc)

            if (split%block_id == m) then
                split%alpha(1:2) = alpha_0(1:2)
            end if
            split%z(:, m) = zloc(split%nmin:split%nmax)
            ! second formulation
            split%y(:, m) = zloc(split%nmin:split%nmax)
            !
        end if

        do m = 1, nblocks - 1       ! loop over remaining coefficients
            p = points(m)           ! index of splitting point
            p_plus_1 = p + 1

            call Splitting(a, b, c, p, p_plus_1, alpha, zloc)

            if (split%block_id == m) then
                split%alpha(1:2) = alpha(1:2)
                split%gamma = alpha_0(1)*zloc(nsize) + alpha_0(2)*zloc(1)
            end if
            split%z(:, m) = zloc(split%nmin:split%nmax)

            ! second formulation
            split%y(:, m) = zloc(split%nmin:split%nmax)
            if (split%circulant) then
                gamma_loc = alpha_0(1)*zloc(nsize) + alpha_0(2)*zloc(1)
                split%y(:, m) = split%y(:, m) + gamma_loc*split%y(:, nblocks)
            end if
            if (m > 1) then
                p = points(m - 1)
                p_plus_1 = p + 1
                beta_previous = alpha_previous(1)*zloc(p) + alpha_previous(2)*zloc(p_plus_1)
                split%y(:, m) = split%y(:, m) + beta_previous*split%y(:, m - 1)
            end if
            alpha_previous(:) = alpha(:)
            !

            if (split%block_id == m - 1) then
                p = points(split%block_id)
                p_plus_1 = p + 1
                split%beta = split%alpha(1)*zloc(p) + split%alpha(2)*zloc(p_plus_1)
            end if

        end do

        ! -------------------------------------------------------------------
        ! block matrix Am and LU decomposition
        split%lhs(:, 1) = aloc(split%nmin:split%nmax)
        split%lhs(:, 2) = bloc(split%nmin:split%nmax)
        split%lhs(:, 3) = cloc(split%nmin:split%nmax)

        deallocate (aloc, bloc, cloc, zloc)

        return

    contains
        subroutine Splitting(a, b, c, p, p_plus_1, alpha, z)
            real(wp), intent(inout) :: a(:), b(:), c(:)
            integer(wi), intent(in) :: p, p_plus_1
            real(wp), intent(inout) :: alpha(2)
            real(wp), intent(inout) :: z(:)

            ! Start definition of alpha
            alpha(1) = a(p_plus_1)
            alpha(2) = c(p)

            ! Generate matrix Ak
            b(p) = b(p) - a(p_plus_1)
            a(p_plus_1) = 0.0_wp
            b(p_plus_1) = b(p_plus_1) - c(p)
            c(p) = 0.0_wp

            ! Generate vector zk
            aloc(:) = a(:); bloc(:) = b(:); cloc(:) = c(:)
            call Thomas3_LU(size(aloc), aloc, bloc, cloc)

            z(:) = 0.0_wp
            z(p) = 1.0_wp
            z(p_plus_1) = 1.0_wp
            call Thomas3_Solve(size(aloc), 1, aloc, bloc, cloc, z(:))

            ! Complete definition of alpha
            delta = 1.0_wp + alpha(1)*z(p) + alpha(2)*z(p_plus_1)
            if (abs(delta) < small_wp) then
                call TLab_Write_ASCII(efile, __FILE__//'. Singular matrix M.')
                call TLab_Stop(DNS_ERROR_THOMAS)
            end if
            alpha(1) = -alpha(1)/delta
            alpha(2) = -alpha(2)/delta

            return
        end subroutine Splitting

    end subroutine Thomas3_Split_Initialize

    !########################################################################
    !########################################################################
    subroutine Thomas3_Split_Solve_Serial(split, f)
        type(thomas3_split_dt), intent(in) :: split(:)
        type(data_dt), intent(inout) :: f(:)

        integer(wi) n, nsize, nlines
        integer(wi) m
        integer(wi) k, nblocks, k_plus_1
        real(wp), allocatable :: tmp(:), xp(:, :), alpha(:, :), alpha_0(:, :)

        !########################################################################
        nblocks = size(split)
        nlines = size(f(1)%p, 1)

        ! Solving block system Am in each block
        ! this part involves lhs and could be out of this routine.
        ! this routine is then called Thomas3_Split_Reduction
        do k = 1, nblocks                       ! loop over blocks
            ! nsize = split(k)%nmax - split(k)%nmin + 1
            nsize = size(f(k)%p, 2)
            call Thomas3_Solve(nsize, nlines, &
                               split(k)%lhs(:, 1), split(k)%lhs(:, 2), split(k)%lhs(:, 3), f(k)%p(:, :))
        end do

        ! -------------------------------------------------------------------
        ! Reduction step

        ! pass x(:,1) to previous block
        allocate (xp(nlines, nblocks))
        do k = nblocks, 1, -1                   ! loop over blocks
            k_plus_1 = mod(k + nblocks, nblocks) + 1
            xp(:, k) = f(k_plus_1)%p(:, 1)      ! solution at left boundary of block k+1
        end do

        allocate (tmp(nlines))
        allocate (alpha(nlines, nblocks))       ! the idea is that each block needs only one alpha
        allocate (alpha_0(nlines, nblocks))     ! the idea is that each block needs only one alpha

        do k = nblocks, 1, -1                   ! loop over blocks
            m = k
            ! nsize = split(k)%nmax - split(k)%nmin + 1
            nsize = size(f(k)%p, 2)
            alpha(:, k) = split(m)%alpha(1)*f(k)%p(:, nsize) + split(m)%alpha(2)*xp(:, k)

            if (k <= nblocks - 2) then          ! the last 2 have no beta correction
                tmp(:) = alpha(:, m + 1)        ! this info comes from block k+1
                alpha(:, k) = alpha(:, k) + split(k)%beta*tmp(:)
            end if

            if (k <= nblocks - 1) then          ! update solution
                ! send alpha_k to all PEs with larger rank
                do m = k, nblocks
                    ! nsize = split(m)%nmax - split(m)%nmin + 1
                    nsize = size(f(m)%p, 2)
                    do n = 1, nsize
                        f(m)%p(:, n) = f(m)%p(:, n) + alpha(:, k)*split(m)%z(n, k)
                    end do
                end do
            end if

            if (split(k)%circulant) then        ! accumulate alpha_0
                if (k <= nblocks - 1) then
                    tmp(:) = alpha_0(:, k + 1)  ! this info comes from block k+1
                    alpha_0(:, k) = tmp(:) + split(k)%gamma*alpha(:, k)
                else
                    alpha_0(:, k) = alpha(:, k)
                end if
            end if

        end do

        if (split(nblocks)%circulant) then
            ! broadcast alpha_0(:) from block 1 into tmp of all ranks
            tmp(:) = alpha_0(:, 1)
            do k = nblocks, 1, -1               ! loop over blocks
                ! nsize = split(k)%nmax - split(k)%nmin + 1
                nsize = size(f(k)%p, 2)
                do n = 1, nsize
                    f(k)%p(:, n) = f(k)%p(:, n) + tmp(:)*split(k)%z(n, nblocks)
                end do
            end do
        end if

        deallocate (tmp, xp, alpha, alpha_0)

        return
    end subroutine Thomas3_Split_Solve_Serial

    !########################################################################
    !########################################################################
    subroutine Thomas3_Split_Solve_Serial2(split, f)
        type(thomas3_split_dt), intent(in) :: split(:)
        type(data_dt), intent(inout) :: f(:)

        integer(wi) n, nsize, nlines
        integer(wi) m, mmax
        integer(wi) k, nblocks, k_plus_1
        real(wp), allocatable :: xp(:, :), alpha(:, :)

        !########################################################################
        nblocks = size(split)
        nlines = size(f(1)%p, 1)

        ! Solving block system Am in each block
        ! this part involves lhs and could be out of this routine.
        ! this routine is then called Thomas3_Split_Reduction
        do k = 1, nblocks                       ! loop over blocks
            ! nsize = split(k)%nmax - split(k)%nmin + 1
            nsize = size(f(k)%p, 2)
            call Thomas3_Solve(nsize, nlines, &
                               split(k)%lhs(:, 1), split(k)%lhs(:, 2), split(k)%lhs(:, 3), f(k)%p(:, :))
        end do

        ! -------------------------------------------------------------------
        ! Reduction step

        ! pass x(:,1) to previous block
        allocate (xp(nlines, nblocks))
        do k = nblocks, 1, -1                   ! loop over blocks
            k_plus_1 = mod(k + nblocks, nblocks) + 1
            xp(:, k) = f(k_plus_1)%p(:, 1)      ! solution at left boundary of block k+1
        end do

        !  calculate local alpha
        allocate (alpha(nlines, nblocks))       ! the idea is that each block needs only one alpha
        do k = nblocks, 1, -1                   ! loop over blocks
            nsize = size(f(k)%p, 2)
            alpha(:, k) = split(k)%alpha(1)*f(k)%p(:, nsize) + split(k)%alpha(2)*xp(:, k)
        end do

        ! send alpha to all blocks

        ! calculate solution
        do k = nblocks, 1, -1                   ! loop over blocks
            nsize = size(f(k)%p, 2)
            if (split(k)%circulant) then
                mmax = nblocks
            else
                mmax = nblocks - 1
            end if

            do m = 1, mmax
                do n = 1, nsize
                    f(k)%p(:, n) = f(k)%p(:, n) + alpha(:, m)*split(k)%y(n, m)
                end do
            end do

        end do

        deallocate (xp, alpha)

        return
    end subroutine Thomas3_Split_Solve_Serial2

#ifdef USE_MPI
    !########################################################################
    !########################################################################
    subroutine Thomas3_Split_Solve_MPI(split, f, alphas, tmp)
        use mpi_f08

        type(thomas3_split_dt), intent(in) :: split
        real(wp), intent(inout) :: f(:, :)
        real(wp), intent(inout) :: alphas(:, :)    ! auxiliary memory space for (alpha_j, alpha_0)
        real(wp), intent(inout) :: tmp(:, :)       ! auxiliary memory space

        integer(wi) n, nsize, nlines
        integer(wi) m
        integer(wi) k, nblocks

        type(MPI_Status) status
        integer ims_err
        integer source, dest, tag
        type(MPI_Request), allocatable :: request(:)
        integer l

        !########################################################################
        nblocks = split%n_ranks
        nlines = size(f, 1)
        nsize = size(f, 2)              ! Assume all blocks have same size

        if (allocated(request)) deallocate (request)
        allocate (request(2*split%n_ranks))

        !########################################################################
        ! Solving block system Am in each block
        ! this part involves lhs and probably should be out of this routine.
        ! this routine is then called Thomas3_Split_Reduction
        call Thomas3_Solve(nsize, nlines, &
                           split%lhs(:, 1), split%lhs(:, 2), split%lhs(:, 3), f(:, :))

        !########################################################################
        ! Reduction step
        ! Assume circulant matrix and need alpha_0

        ! -------------------------------------------------------------------
        ! pass x(:,1) to previous block
#define xp(j) tmp(:,1)

        dest = mod(split%rank - 1 + split%n_ranks, split%n_ranks)
        source = mod(split%rank + 1, split%n_ranks)
        tag = 0
        call MPI_Sendrecv(f(:, 1), nlines, MPI_REAL8, dest, tag, &
                          xp(:), nlines, MPI_REAL8, source, tag, &
                          split%communicator, status, ims_err)

        ! -------------------------------------------------------------------
        ! Calculate coefficients

        nsize = size(f, 2)
        alphas(:, 1) = split%alpha(1)*f(:, nsize) + split%alpha(2)*xp(:)    ! alpha_j, first contribution
        ! if (split%circulant) then
        alphas(:, 2) = alphas(:, 1)                                         ! alpha_0, initial value
        ! end if
#undef xp

        tag = 1
        if (split%rank < nblocks - 1) then
            source = split%rank + 1
            call MPI_Recv(tmp, 2*nlines, MPI_REAL8, source, tag, split%communicator, status, ims_err)
            if (split%block_id <= nblocks - 2) then
                alphas(:, 1) = alphas(:, 1) + split%beta*tmp(:, 1)          ! alpha_j, second contribution
            end if
            ! if (split%circulant) then
            alphas(:, 2) = tmp(:, 2) + split%gamma*alphas(:, 1)             ! alpha_0, accumulation
            !endif
        end if

        if (split%rank /= 0) then
            dest = split%rank - 1
            call MPI_Send(alphas, 2*nlines, MPI_REAL8, dest, tag, split%communicator, ims_err)
        end if

        l = 0                                               ! send coefficients to processors with higher rank
        if (split%block_id <= nblocks - 1) then
            do dest = split%rank + 1, split%n_ranks - 1
                l = l + 1
                tag = split%block_id
                call MPI_ISend(alphas(:, 1), nlines, MPI_REAL8, dest, tag, split%communicator, request(l), ims_err)
            end do
        end if

        ! -------------------------------------------------------------------
        ! update solution
        if (split%block_id <= nblocks - 1) then             ! local coefficient
            m = split%block_id
            do n = 1, nsize
                f(:, n) = f(:, n) + alphas(:, 1)*split%z(n, m)
            end do
        end if

        if (split%rank /= 0) then                              ! coefficients from processors with lower rank
            do source = split%rank - 1, 0, -1
                l = l + 1
                tag = source + 1
                call MPI_IRecv(tmp(:, 1), nlines, MPI_REAL8, source, tag, split%communicator, request(l), ims_err)
                call MPI_Wait(request(l), status, ims_err)
                ! call MPI_Recv(tmp(:, 1), nlines, MPI_REAL8, source, tag, split%communicator, status, ims_err)
                m = tag
                do n = 1, nsize
                    f(:, n) = f(:, n) + tmp(:, 1)*split%z(n, m)
                end do
            end do
        end if

        ! if (split%circulant) then                         ! circulant coefficient
        call MPI_Bcast(alphas(:, 2), nlines, MPI_REAL8, 0, split%communicator, ims_err)
        do n = 1, nsize
            f(:, n) = f(:, n) + alphas(:, 2)*split%z(n, nblocks)
        end do
        ! end if

        return
    end subroutine Thomas3_Split_Solve_MPI

    !########################################################################
    !########################################################################
    subroutine Thomas3_Split_Solve_MPI2(split, f, alpha, tmp)
        use mpi_f08

        type(thomas3_split_dt), intent(in) :: split
        real(wp), intent(inout) :: f(:, :)
        real(wp), intent(inout) :: alpha(:)         ! auxiliary memory space for local alpha
        real(wp), intent(inout) :: tmp(:, :)        ! auxiliary memory space for all alphas

        integer(wi) n, nsize, nlines
        integer(wi) m
        integer(wi) k, nblocks

        type(MPI_Status) status
        integer ims_err
        integer source, dest, tag
        type(MPI_Request), allocatable :: request(:)
        integer l

        !########################################################################
        nblocks = split%n_ranks
        nlines = size(f, 1)
        nsize = size(f, 2)              ! Assume all blocks have same size

        if (allocated(request)) deallocate (request)
        allocate (request(2*split%n_ranks))

        !########################################################################
        ! Solving block system Am in each block
        ! this part involves lhs and probably should be out of this routine.
        ! this routine is then called Thomas3_Split_Reduction
        call Thomas3_Solve(nsize, nlines, &
                           split%lhs(:, 1), split%lhs(:, 2), split%lhs(:, 3), f(:, :))

        !########################################################################
        ! Reduction step
        ! Assume circulant matrix and need alpha_0

        ! -------------------------------------------------------------------
        ! pass x(:,1) to previous block and calculate local coefficient
#define xp(j) alpha(j)

        dest = mod(split%rank - 1 + split%n_ranks, split%n_ranks)
        source = mod(split%rank + 1, split%n_ranks)
        tag = 0
        call MPI_ISend(f(:, 1), nlines, MPI_REAL8, dest, tag, &
                       split%communicator, request(1), ims_err)
        call MPI_IRecv(xp(:), nlines, MPI_REAL8, source, tag, &
                       split%communicator, request(2), ims_err)
        call MPI_Wait(request(2), status, ims_err)
        alpha(:) = split%alpha(1)*f(:, nsize) + split%alpha(2)*xp(:)

#undef xp

        ! Distribute coefficients
        call MPI_Allgather(alpha, nlines, MPI_REAL8, &
                           tmp, nlines, MPI_REAL8, split%communicator, ims_err)

        ! -------------------------------------------------------------------
        ! update solution as the coefficient is received
        do m = 1, nblocks
            do n = 1, nsize
                f(:, n) = f(:, n) + tmp(:, m)*split%y(n, m)
            end do
        end do

        return
    end subroutine Thomas3_Split_Solve_MPI2
#endif

end module Thomas3_Split
