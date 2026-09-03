module OPR_Burgers
    use TLab_Constants, only: wp, wi
    implicit none
    private

    public :: burgers1d         ! polymorphic

    public :: burgers1d_boussinesq
    public :: burgers1d_anelastic
    ! public :: burgers1d_compressible

    public :: burgers1d_subsidence_boussinesq
    public :: burgers1d_subsidence_anelastic
    ! public :: burgers1d_subsidence_compressible

    ! -----------------------------------------------------------------------
    type, abstract :: burgers1d
        real(wp) :: diffusivity
    contains
        procedure :: initialize => burgers1d_initialize
        procedure :: compute => burgers1d_compute
        procedure :: add => burgers1d_add
        !
        procedure(compute_setrhou_ice), deferred :: compute_setrhou     ! handle density contribution depending on type
    end type

    type, extends(burgers1d) :: burgers1d_boussinesq
    contains
        procedure :: compute_setrhou => boussinesq_compute_setrho
    end type

    type, extends(burgers1d) :: burgers1d_anelastic
        real(wp), allocatable :: rho(:)
    contains
        procedure :: initialize => burgers1d_anelastic_initialize
        procedure :: compute_setrhou => anelastic_compute_setrho
    end type

    abstract interface
        subroutine compute_setrhou_ice(self, nlines, nsize, der1, der2, rhou)
            import wp, wi, burgers1d
            class(burgers1d) self
            integer(wi), intent(in) :: nlines, nsize
            real(wp), intent(in) :: der1(nlines, nsize)
            real(wp), intent(inout) :: der2(nlines, nsize)
            real(wp), intent(inout) :: rhou(nlines, nsize)
        end subroutine compute_setrhou_ice
    end interface

    ! -----------------------------------------------------------------------
    type, abstract, extends(burgers1d) :: burgers1d_subsidence
        real(wp), allocatable :: rhou_background(:)
    contains
        procedure :: initialize => burgers1d_subsidence_initialize
        procedure :: add => burgers1d_subsidence_add
    end type

    type, extends(burgers1d_subsidence) :: burgers1d_subsidence_boussinesq
    contains
        procedure :: compute_setrhou => boussinesq_subsidence_compute_setrho
    end type

    type, extends(burgers1d_subsidence) :: burgers1d_subsidence_anelastic
        real(wp), allocatable :: rho(:)
    contains
        procedure :: initialize => burgers1d_subsidence_anelastic_initialize
        procedure :: compute_setrhou => anelastic_subsidence_compute_setrho
    end type

contains
    !########################################################################
    !########################################################################
    subroutine burgers1d_initialize(self, diffusivity, axis, rbackground, wbackground)
        class(burgers1d), intent(out) :: self
        real(wp), intent(in) :: diffusivity
        character(len=*), intent(in), optional :: axis
        real(wp), intent(in), optional :: rbackground(:)
        real(wp), intent(in), optional :: wbackground(:)

        self%diffusivity = diffusivity

        return
    end subroutine

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
    end subroutine

    subroutine burgers1d_add(self, nlines, nsize, der1, der2, rhou, result)
        class(burgers1d) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(in) :: der2(nlines, nsize)
        real(wp), intent(in) :: rhou(nlines, nsize)
        real(wp), intent(out) :: result(nlines, nsize)

        result(:, :) = result(:, :) + der2(:, :)*self%diffusivity - rhou(:, :)*der1(:, :)

        return
    end subroutine

    subroutine burgers1d_subsidence_initialize(self, diffusivity, axis, rbackground, wbackground)
        class(burgers1d_subsidence), intent(out) :: self
        real(wp), intent(in) :: diffusivity
        character(len=*), intent(in), optional :: axis
        real(wp), intent(in), optional :: rbackground(:)
        real(wp), intent(in), optional :: wbackground(:)

        self%diffusivity = diffusivity
        allocate (self%rhou_background, source=wbackground)

        return
    end subroutine

    subroutine burgers1d_subsidence_add(self, nlines, nsize, der1, der2, rhou, result)
        class(burgers1d_subsidence) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(in) :: der2(nlines, nsize)
        real(wp), intent(in) :: rhou(nlines, nsize)
        real(wp), intent(out) :: result(nlines, nsize)

        integer n

        do n = 1, nsize
            result(:, n) = result(:, n) + der2(:, n)*self%diffusivity + (self%rhou_background(n) - rhou(:, n))*der1(:, n)
        end do

        return
    end subroutine

    !########################################################################
    !########################################################################
    ! Handle rho contribution depending on type
    subroutine boussinesq_compute_setrho(self, nlines, nsize, der1, der2, rhou)
        ! same as compute but with inout type for rhou
        class(burgers1d_boussinesq) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(inout) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)

#define result(i,j) der2(i,j)

        result(:, :) = der2(:, :)*self%diffusivity - rhou(:, :)*der1(:, :)

#undef result

        return
    end subroutine

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
    end subroutine

    subroutine boussinesq_subsidence_compute_setrho(self, nlines, nsize, der1, der2, rhou)
        ! same as compute but with inout type for rhou
        class(burgers1d_subsidence_boussinesq) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(inout) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)

        ! TBD

        return
    end subroutine

    subroutine anelastic_subsidence_compute_setrho(self, nlines, nsize, der1, der2, rhou)
        class(burgers1d_subsidence_anelastic) self
        integer(wi), intent(in) :: nlines, nsize
        real(wp), intent(in) :: der1(nlines, nsize)
        real(wp), intent(inout) :: der2(nlines, nsize)
        real(wp), intent(inout) :: rhou(nlines, nsize)

        ! TBD

        return
    end subroutine

    !########################################################################
    !########################################################################
    subroutine burgers1d_anelastic_initialize(self, diffusivity, axis, rbackground, wbackground)
        class(burgers1d_anelastic), intent(out) :: self
        real(wp), intent(in) :: diffusivity
        character(len=*), intent(in), optional :: axis
        real(wp), intent(in), optional :: rbackground(:)
        real(wp), intent(in), optional :: wbackground(:)

        self%diffusivity = diffusivity

        call anelastic_initialize_rho(self%rho, axis, rbackground)

        return
    end subroutine burgers1d_anelastic_initialize

    subroutine burgers1d_subsidence_anelastic_initialize(self, diffusivity, axis, rbackground, wbackground)
        class(burgers1d_subsidence_anelastic), intent(out) :: self
        real(wp), intent(in) :: diffusivity
        character(len=*), intent(in), optional :: axis
        real(wp), intent(in), optional :: rbackground(:)
        real(wp), intent(in), optional :: wbackground(:)

        self%diffusivity = diffusivity

        call anelastic_initialize_rho(self%rho, axis, rbackground)

        allocate (self%rhou_background, source=wbackground)
        self%rhou_background(:) = self%rhou_background(:)*rbackground(:)
        ! call anelastic_initialize_rho(self%rhou_background, axis, rbackground)

        return
    end subroutine burgers1d_subsidence_anelastic_initialize

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
