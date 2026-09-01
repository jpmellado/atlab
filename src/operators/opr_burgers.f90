module OPR_Burgers
    use TLab_Constants, only: wp, wi
    implicit none
    private

    public :: burgers1d
    public :: burgers1d_boussinesq
    public :: burgers1d_anelastic

    ! -----------------------------------------------------------------------
    type, abstract :: burgers1d
        real(wp) :: diffusivity
    contains
        procedure :: initialize => burgers1d_initialize
        procedure :: compute => burgers1d_compute
        procedure :: compute_setrhou => burgers1d_compute_setrhou
    end type

    type, extends(burgers1d) :: burgers1d_boussinesq
    contains
    end type

    type, extends(burgers1d) :: burgers1d_anelastic
        real(wp), allocatable :: rho(:)
    contains
    end type

contains
    !########################################################################
    !########################################################################
    subroutine burgers1d_initialize(self, diffusivity, axis, rbackground)
        class(burgers1d), intent(out) :: self
        real(wp), intent(in) :: diffusivity
        character(len=*), intent(in) :: axis
        real(wp), intent(in), optional :: rbackground(:)

        self%diffusivity = diffusivity

        select type (self)
        type is (burgers1d_anelastic)
            call anelastic_initialize_rho(self%rho, axis, rbackground)
        end select

        return

    end subroutine burgers1d_initialize

!########################################################################
!########################################################################
    subroutine burgers1d_compute(self, nlines, nsize, der1, der2, rhou)
        class(burgers1d) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(inout) :: der2(nlines, nsize)
        real(wp), intent(in) :: rhou(nlines, nsize)

#define result(i,j) der2(i,j)

        result(:, :) = der2(:, :)*self%diffusivity - rhou(:, :)*der1(:, :)

#undef result

        return
    end subroutine burgers1d_compute

!########################################################################
!########################################################################
    subroutine burgers1d_compute_setrhou(self, nlines, nsize, der1, der2, rhou)
        class(burgers1d) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(inout) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)

        select type (self)
        type is (burgers1d_anelastic)
            call anelastic_compute_setrho(self, nlines, nsize, der1, der2, rhou)

        type is (burgers1d_boussinesq)
            call self%compute(nlines, nsize, der1, der2, rhou)

        end select

        return
    end subroutine burgers1d_compute_setrhou

!########################################################################
!########################################################################
    subroutine anelastic_compute_setrho(self, nlines, nsize, der1, der2, rhou)
        class(burgers1d_anelastic) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(inout) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)

        integer ij

#define result(i,j) der2(i,j)

        do ij = 1, nsize
            rhou(:, ij) = rhou(:, ij)*self%rho(:)
            result(:, ij) = der2(:, ij)*self%diffusivity - rhou(:, ij)*der1(:, ij)
        end do

#undef result

        return
    end subroutine anelastic_compute_setrho

!########################################################################
!########################################################################
    subroutine anelastic_initialize_rho(rho, axis, rbackground)
        use TLab_Memory, only: imax, jmax, kmax
#ifdef USE_MPI
        use TLabMPI_VARS, only: xMpi, yMpi
        use TLabMPI_Transpose_DerivedTypes, only: tmpi_trp_X, tmpi_trp_Y
#endif
        use TLab_Grid, only: z
#ifdef USE_MPI
        use OPR_Partial, only: der_mode_i, der_mode_j, TYPE_TRANSPOSE
#endif

        real(wp), intent(in) :: rbackground(:)
        real(wp), allocatable, intent(out) :: rho(:)
        character(len=*), intent(in) :: axis

        integer(wi) ip, j
        integer(wi) nlines, offset

        !########################################################################
        select case (trim(adjustl(axis)))
            ! -----------------------------------------------------------------------
            ! Density correction term in the burgers operator along X
        case ('x')
#ifdef USE_MPI
            if (xMpi%num_processors > 1 .and. der_mode_i == TYPE_TRANSPOSE) then
                ! nlines = tmpi_plan_dx%nlines
                nlines = tmpi_trp_X%nlines
                offset = nlines*xMpi%rank
            else
#endif
                nlines = jmax*kmax
                offset = 0
#ifdef USE_MPI
            end if
#endif
            allocate (rho(nlines))
            do j = 1, nlines
                ip = (offset + j - 1)/jmax + 1
                rho(j) = rbackground(ip)
            end do

            ! -----------------------------------------------------------------------
            ! Density correction term in the burgers operator along Y
        case ('y')
#ifdef USE_MPI
            if (yMpi%num_processors > 1 .and. der_mode_j == TYPE_TRANSPOSE) then
                ! nlines = tmpi_plan_dy%nlines
                nlines = tmpi_trp_Y%nlines
                offset = nlines*yMpi%rank
            else
#endif
                nlines = imax*kmax
                offset = 0
#ifdef USE_MPI
            end if
#endif
            allocate (rho(nlines))
            do j = 1, nlines
                ip = mod(offset + j - 1, z%size) + 1
                rho(j) = rbackground(ip)
            end do

        end select

        return
    end subroutine anelastic_initialize_rho

end module
