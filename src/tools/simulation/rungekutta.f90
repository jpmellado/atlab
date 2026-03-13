! Integrate d_t a = b(a,t) from t to t + delta t

module RungeKutta
    use TLab_Constants, only: wp, wi
    implicit none
    private

    public :: rungekutta_lowstorage_dt
    public :: rk3_dt
    public :: rk45_dt

    ! -----------------------------------------------------------------------
    type, abstract :: rungekutta_lowstorage_dt
        integer :: num_substep
        real(wp), allocatable :: coef_a(:)
        real(wp), allocatable :: coef_b(:)
        real(wp), allocatable :: coef_t(:)
        integer :: substep
        real(wp) :: time                    ! local time at each substep of the RK step
    contains
        procedure(initialize_ice), deferred :: initialize
        procedure :: AdvanceSubstep_Boussinesq
        procedure :: AdvanceSubstep_Anelastic
    end type
    abstract interface
        subroutine initialize_ice(self)
            import rungekutta_lowstorage_dt
            class(rungekutta_lowstorage_dt), intent(out) :: self
        end subroutine
    end interface

    type, extends(rungekutta_lowstorage_dt) :: rk3_dt
    contains
        procedure :: initialize => rk3_initialize
    end type

    type, extends(rungekutta_lowstorage_dt) :: rk45_dt
    contains
        procedure :: initialize => rk45_initialize
    end type

contains
    ! ###################################################################
    ! ###################################################################
    subroutine rk3_initialize(self)
        class(rk3_dt), intent(out) :: self

        integer i

        self%num_substep = 3
        allocate (self%coef_a(1:self%num_substep))
        allocate (self%coef_b(1:self%num_substep))
        allocate (self%coef_t(1:self%num_substep))

        self%coef_a(1:3) = [1.0_wp/3.0_wp, &
                            15.0_wp/16.0_wp, &
                            8.0_wp/15.0_wp]

        self%coef_b(1:3) = [-5.0_wp/9.0_wp, &
                            -153.0_wp/128.0_wp, &
                            0.0_wp]

        self%coef_t(1:3) = [1.0_wp/3.0_wp, &
                            3.0_wp/4.0_wp, &
                            1.0_wp]
        ! go to increments
        do i = self%num_substep, 2, -1
            self%coef_t(i) = self%coef_t(i) - self%coef_t(i - 1)
        end do

        return
    end subroutine

    ! ###################################################################
    ! ###################################################################
    subroutine rk45_initialize(self)
        class(rk45_dt), intent(out) :: self

        integer i

        self%num_substep = 5
        allocate (self%coef_a(1:self%num_substep))
        allocate (self%coef_b(1:self%num_substep))
        allocate (self%coef_t(1:self%num_substep))

        self%coef_a(1) = 1432997174477.0_wp/9575080441755.0_wp
        self%coef_a(2) = 5161836677717.0_wp/13612068292357.0_wp
        self%coef_a(3) = 1720146321549.0_wp/2090206949498.0_wp
        self%coef_a(4) = 3134564353537.0_wp/4481467310338.0_wp
        self%coef_a(5) = 2277821191437.0_wp/14882151754819.0_wp

        self%coef_b(1) = -567301805773.0_wp/1357537059087.0_wp
        self%coef_b(2) = -2404267990393.0_wp/2016746695238.0_wp
        self%coef_b(3) = -3550918686646.0_wp/2091501179385.0_wp
        self%coef_b(4) = -1275806237668.0_wp/842570457699.0_wp
        self%coef_b(5) = 0.0_wp

        self%coef_t(1) = self%coef_a(1)
        self%coef_t(2) = 2526269341429.0_wp/6820363962896.0_wp
        self%coef_t(3) = 2006345519317.0_wp/3224310063776.0_wp
        self%coef_t(4) = 2802321613138.0_wp/2924317926251.0_wp
        self%coef_t(5) = 1.0_wp
        ! go to increments
        do i = self%num_substep, 2, -1
            self%coef_t(i) = self%coef_t(i) - self%coef_t(i - 1)
        end do

        return
    end subroutine

    ! ###################################################################
    ! ###################################################################
    ! Memory intensive, cache-optimized
    ! create size of memory blocks during initialization
    subroutine AdvanceSubstep_Boussinesq(self, a, b, delta_t)
        class(rungekutta_lowstorage_dt) self
        real(wp), intent(inout) :: a(:, :, :)
        real(wp), intent(inout) :: b(:, :, :)
        real(wp), intent(in) :: delta_t

        integer k, kmax, iq

        kmax = size(b, 2)

        if (self%substep < self%num_substep) then

            do iq = 1, size(b, 3)
                do k = 1, kmax
                    a(:, k, iq) = a(:, k, iq) + self%coef_a(self%substep)*delta_t*b(:, k, iq)
                    b(:, k, iq) = b(:, k, iq)*self%coef_b(self%substep)
                end do
            end do

        else
            do iq = 1, size(b, 3)
                do k = 1, kmax
                    a(:, k, iq) = a(:, k, iq) + self%coef_a(self%substep)*delta_t*b(:, k, iq)
                end do
            end do

        end if

        return
    end subroutine

    ! ###################################################################
    ! ###################################################################
    subroutine AdvanceSubstep_Anelastic(self, a, b, delta_t, weight)
        class(rungekutta_lowstorage_dt) self
        real(wp), intent(inout) :: a(:, :, :)
        real(wp), intent(inout) :: b(:, :, :)
        real(wp), intent(in) :: delta_t
        real(wp), intent(in) :: weight(:)

        integer k, kmax, iq

        kmax = size(b, 2)

        if (self%substep < self%num_substep) then

            do iq = 1, size(b, 3)
                do k = 1, kmax
                    a(:, k, iq) = a(:, k, iq) + self%coef_a(self%substep)*delta_t*b(:, k, iq)*weight(k)
                    b(:, k, iq) = b(:, k, iq)*self%coef_b(self%substep)
                end do
            end do

        else
            do iq = 1, size(b, 3)
                do k = 1, kmax
                    a(:, k, iq) = a(:, k, iq) + self%coef_a(self%substep)*delta_t*b(:, k, iq)*weight(k)
                end do
            end do

        end if

        return
    end subroutine

end module
