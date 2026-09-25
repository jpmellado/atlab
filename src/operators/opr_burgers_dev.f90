! Operator \mu d^2 s - [ rho u - (rho u)_background ] d s
module OPR_Burgers_Dev
    use TLab_Constants, only: wp, wi
    implicit none
    private

    public :: burgers1d         ! polymorphic

    public :: burgers1d_XY
    public :: burgers1d_Z

    public :: burgers1d_subsidence_XY
    public :: burgers1d_subsidence_Z

    ! -----------------------------------------------------------------------
    type :: burgers1d_base
        real(wp) :: diffusivity                         ! coefficient mu
    contains
        procedure :: initialize => burgers1d_initialize
        procedure :: compute => burgers1d_compute
        procedure :: add => burgers1d_add               ! compute the operator and add it to a rhs array
    end type

    type, abstract, extends(burgers1d_base) :: burgers1d  ! handle form of advection field
        real(wp), allocatable :: rho(:)
        real(wp), allocatable :: rhou_background(:)
    contains
        procedure :: initialize_setrho => burgers1d_initialize_setrho
        procedure :: compute_setrhou => burgers1s_compute_setrhou
        procedure :: add_setrhou => burgers1s_add_setrhou
    end type

    ! -----------------------------------------------------------------------
    ! Subroutines that include a rho term
    type, extends(burgers1d) :: burgers1d_XY
    contains
        procedure :: compute_setrhou => compute_setrhou_XY
        procedure :: add_setrhou => add_setrhou_XY
    end type

    type, extends(burgers1d) :: burgers1d_Z
    contains
        procedure :: compute_setrhou => compute_setrhou_Z
        procedure :: add_setrhou => add_setrhou_Z
    end type

    ! -----------------------------------------------------------------------
    ! Subroutines that include a subsidence term to reduce memory calls
    type, extends(burgers1d) :: burgers1d_subsidence_XY
    contains
        procedure :: compute_setrhou => compute_setrhou_subsidence_XY
        procedure :: add_setrhou => add_setrhou_subsidence_XY
    end type

    type, extends(burgers1d) :: burgers1d_subsidence_Z
    contains
        procedure :: compute_setrhou => compute_setrhou_subsidence_Z
        procedure :: add_setrhou => add_setrhou_subsidence_Z
    end type

contains
    !########################################################################
    !########################################################################
    subroutine burgers1d_initialize(self, diffusivity)
        class(burgers1d_base), intent(out) :: self
        real(wp), intent(in) :: diffusivity

        self%diffusivity = diffusivity

        return
    end subroutine

    subroutine burgers1d_compute(self, nlines, nsize, der1, der2, rhou)
        class(burgers1d_base) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(inout) :: der2(nlines, nsize)
        real(wp), intent(in) :: rhou(nlines, nsize)

#define result(i,j) der2(i,j)

        result(:, :) = der2(:, :)*self%diffusivity - rhou(:, :)*der1(:, :)

#undef result

        return
    end subroutine

    subroutine burgers1d_add(self, nlines, nsize, der1, der2, rhou, result)
        class(burgers1d_base) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(in) :: der2(nlines, nsize)
        real(wp), intent(in) :: rhou(nlines, nsize)
        real(wp), intent(out) :: result(nlines, nsize)

        result(:, :) = result(:, :) + der2(:, :)*self%diffusivity - rhou(:, :)*der1(:, :)

        return
    end subroutine

    !########################################################################
    !########################################################################
    ! Handle different forms of advection term
    subroutine burgers1d_initialize_setrho(self, diffusivity, axis, rbackground, wbackground)
        class(burgers1d), intent(out) :: self
        real(wp), intent(in) :: diffusivity
        character(len=*), intent(in), optional :: axis
        real(wp), intent(in), optional :: rbackground(:)
        real(wp), intent(in), optional :: wbackground(:)

        call self%initialize(diffusivity)

        if (present(rbackground)) call anelastic_initialize_rho(self%rho, axis, rbackground)
        if (present(wbackground)) allocate (self%rhou_background, source=wbackground)
        if (allocated(self%rho) .and. allocated(self%rhou_background)) then
            self%rhou_background(:) = self%rhou_background(:)*self%rho(:)
        end if

        return
    end subroutine

    ! -----------------------------------------------------------------------
    ! Wrappers for the case in which the advection field is simply u
    subroutine burgers1s_compute_setrhou(self, nlines, nsize, der1, der2, rhou)
        class(burgers1d) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(inout) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)

        call burgers1d_compute(self, nlines, nsize, der1, der2, rhou)

        return
    end subroutine

    subroutine burgers1s_add_setrhou(self, nlines, nsize, der1, der2, rhou, result)
        class(burgers1d) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(in) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)
        real(wp), intent(out) :: result(nlines, nsize)

        call burgers1d_add(self, nlines, nsize, der1, der2, rhou, result)

        return
    end subroutine

    ! -----------------------------------------------------------------------
    subroutine compute_setrhou_XY(self, nlines, nsize, der1, der2, rhou)
        class(burgers1d_XY) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(inout) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)

        integer n

#define result(i,j) der2(i,j)
        do n = 1, nsize
            rhou(:, n) = rhou(:, n)*self%rho(:)
            result(:, n) = der2(:, n)*self%diffusivity - rhou(:, n)*der1(:, n)
        end do
#undef result

        return
    end subroutine

    subroutine add_setrhou_XY(self, nlines, nsize, der1, der2, rhou, result)
        class(burgers1d_XY) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(in) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)
        real(wp), intent(out) :: result(nlines, nsize)

        integer n

        do n = 1, nsize
            rhou(:, n) = rhou(:, n)*self%rho(:)
            result(:, n) = result(:, n) + der2(:, n)*self%diffusivity - rhou(:, n)*der1(:, n)
        end do

        return
    end subroutine

    subroutine compute_setrhou_Z(self, nlines, nsize, der1, der2, rhou)
        class(burgers1d_Z) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(inout) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)

        integer n

#define result(i,j) der2(i,j)
        do n = 1, nsize
            rhou(:, n) = rhou(:, n)*self%rho(n)
            result(:, n) = der2(:, n)*self%diffusivity - rhou(:, n)*der1(:, n)
        end do
#undef result

        return
    end subroutine

    subroutine add_setrhou_Z(self, nlines, nsize, der1, der2, rhou, result)
        class(burgers1d_Z) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(in) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)
        real(wp), intent(out) :: result(nlines, nsize)

        integer n

        do n = 1, nsize
            rhou(:, n) = rhou(:, n)*self%rho(n)
            result(:, n) = result(:, n) + der2(:, n)*self%diffusivity - rhou(:, n)*der1(:, n)
        end do

        return
    end subroutine

    ! -----------------------------------------------------------------------
    subroutine compute_setrhou_subsidence_XY(self, nlines, nsize, der1, der2, rhou)
        class(burgers1d_subsidence_XY) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(inout) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)

        integer n

#define result(i,j) der2(i,j)
        do n = 1, nsize
            rhou(:, n) = rhou(:, n)*self%rho(:) - self%rhou_background(:)
            result(:, n) = der2(:, n)*self%diffusivity - rhou(:, n)*der1(:, n)
        end do
#undef result

        return
    end subroutine

    subroutine add_setrhou_subsidence_XY(self, nlines, nsize, der1, der2, rhou, result)
        class(burgers1d_subsidence_XY) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(in) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)
        real(wp), intent(out) :: result(nlines, nsize)

        integer n

        do n = 1, nsize
            rhou(:, n) = rhou(:, n)*self%rho(:) - self%rhou_background(:)
            result(:, n) = result(:, n) + der2(:, n)*self%diffusivity - rhou(:, n)*der1(:, n)
        end do

        return
    end subroutine

    subroutine compute_setrhou_subsidence_Z(self, nlines, nsize, der1, der2, rhou)
        class(burgers1d_subsidence_Z) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(inout) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)

        integer n

#define result(i,j) der2(i,j)
        do n = 1, nsize
            rhou(:, n) = rhou(:, n)*self%rho(n) - self%rhou_background(n)
            result(:, n) = der2(:, n)*self%diffusivity - rhou(:, n)*der1(:, n)
        end do
#undef result

        return
    end subroutine

    subroutine add_setrhou_subsidence_Z(self, nlines, nsize, der1, der2, rhou, result)
        class(burgers1d_subsidence_Z) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(in) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)
        real(wp), intent(out) :: result(nlines, nsize)

        integer n

        do n = 1, nsize
            rhou(:, n) = rhou(:, n)*self%rho(n) - self%rhou_background(n)
            result(:, n) = result(:, n) + der2(:, n)*self%diffusivity - rhou(:, n)*der1(:, n)
        end do

        return
    end subroutine

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

            ! -----------------------------------------------------------------------
            ! Density correction term in the burgers operator along Z
        case ('z')
            allocate (rho, source=rbackground)

        end select

        return
    end subroutine anelastic_initialize_rho

end module
