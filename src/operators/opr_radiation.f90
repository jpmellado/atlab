! Solve radiative transfer equation along 1 direction

module OPR_Radiation
    use TLab_Constants, only: wp, wi, pi_wp
    use TLab_Constants, only: BCS_MAX, BCS_MIN, BCS_BOTH
    use TLab_Arrays, only: wrk2d
    use TLab_Grid, only: z
    use OPR_ODES
    use FDM_Integral_1
    use Integration
    implicit none
    private

    public :: IR_RTE1_OnlyLiquid
    public :: IR_RTE1_OnlyLiquid_UpwardFirst
    public :: IR_RTE1_Global
    public :: IR_RTE1_Incremental
    public :: IR_RTE1_Local

    real(wp), parameter, public :: sigma = &    ! Stefan-Boltzmann constant, W /m^2 /K^4
                                   5.67037442e-8_wp

    real(wp), parameter, public :: mu = &       ! mean direction parameter
                                   0.5_wp*(1.0_wp/sqrt(3.0_wp) + 1.0_wp/sqrt(2.0_wp))  ! mean of (1/sqrt{3}, 1/sqrt{2})
    !                               1.0_wp/sqrt(2.0_wp)
    !                               0.5_wp       ! testing

    real(wp), public :: epsilon                 ! surface emissivity at the bottom of the domain

contains
    !########################################################################
    !########################################################################
    subroutine IR_RTE1_OnlyLiquid(nlines, ny, fdmi, a_source, tau, bcs_hb, bcs_ht, ibc, flux_down, flux_up)
        integer(wi), intent(in) :: nlines, ny
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(inout) :: a_source(nlines, ny)     ! input as bulk absorption coefficient, output as source
        real(wp), intent(inout) :: tau(nlines, ny)          ! optical depth and transmission functions
        real(wp), intent(inout) :: bcs_hb(nlines)           ! flux boundary condition at the bottom
        real(wp), intent(inout) :: bcs_ht(nlines)           ! flux boundary condition at the top
        integer :: ibc                                      ! boundary condition: top, down, both
        real(wp), intent(out), optional :: flux_down(nlines, ny), flux_up(nlines, ny)

        ! -----------------------------------------------------------------------
        integer(wi) j

        ! #######################################################################
        ! calculate f_j = exp(-tau(z, zmax)/\mu)
        tau(:, ny) = 0.0_wp                                   ! boundary condition
        call FDM_Int1_Solve(nlines, fdmi(BCS_MAX), fdmi(BCS_MAX)%rhs, a_source, tau, wrk2d)       ! recall this gives the negative of the integral
        ! call Int_Trapezoidal_f(a_source, z%nodes(:), tau, BCS_MAX)
        ! call Int_Simpson_Biased_f(a_source, z%nodes(:), tau, BCS_MAX)
        do j = ny, 1, -1
            tau(:, j) = exp(tau(:, j))
        end do
        !  tau = dexp(tau)         seg-fault; need ulimit -u unlimited

        ! Calculate heating rate
        if (ibc == BCS_BOTH) then
            do j = ny, 1, -1
                a_source(:, j) = a_source(:, j)*(tau(:, j)*bcs_ht(1:nlines) &                     ! downward flux
                                                 + tau(:, 1)/tau(:, j)*bcs_hb(1:nlines))        ! upward flux
            end do
        else
            do j = ny, 1, -1
                a_source(:, j) = a_source(:, j)*tau(:, j)*bcs_ht(1:nlines)                        ! only downward
            end do
        end if

        ! Calculate flux, if necessary
        if (present(flux_up)) then
            do j = ny, 1, -1
                flux_down(:, j) = bcs_ht(1:nlines)*tau(:, j)                  ! downward flux
                flux_up(:, j) = bcs_hb(1:nlines)*tau(:, 1)/tau(:, j)        ! upward flux
            end do

        end if

        return
    end subroutine IR_RTE1_OnlyLiquid

    !########################################################################
    !########################################################################
    subroutine IR_RTE1_OnlyLiquid_UpwardFirst(nlines, ny, fdmi, a_source, tau, bcs_hb, bcs_ht, ibc, flux_down, flux_up)
        integer(wi), intent(in) :: nlines, ny
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(inout) :: a_source(nlines, ny)     ! input as bulk absorption coefficent, output as source
        real(wp), intent(inout) :: tau(nlines, ny)          ! optical depth and transmission functions
        real(wp), intent(inout) :: bcs_hb(nlines)           ! flux boundary condition at the bottom
        real(wp), intent(inout) :: bcs_ht(nlines)           ! flux boundary condition at the top
        integer :: ibc                                      ! boundary condition: top, down, both
        real(wp), intent(out), optional :: flux_down(nlines, ny), flux_up(nlines, ny)

        ! -----------------------------------------------------------------------
        integer(wi) j

        ! #######################################################################
        ! calculate f_j = exp(-tau(0, z)/\mu)
        tau(:, 1) = 0.0_wp                                   ! boundary condition
        call FDM_Int1_Solve(nlines, fdmi(BCS_MIN), fdmi(BCS_MIN)%rhs, a_source, tau, wrk2d)
        do j = 1, ny
            tau(:, j) = exp(-tau(:, j))
        end do

        ! Calculate heating rate
        if (ibc == BCS_BOTH) then
            do j = ny, 1, -1
                a_source(:, j) = a_source(:, j)*(tau(:, ny)/tau(:, j)*bcs_ht(1:nlines) &    ! downward flux
                                                 + tau(:, j)*bcs_hb(1:nlines))                ! upward flux
            end do
        else
            do j = ny, 1, -1
                a_source(:, j) = a_source(:, j)*tau(:, ny)/tau(:, j)*bcs_ht(1:nlines)       ! only downward
            end do
        end if

        ! Calculate flux, if necessary
        if (present(flux_up)) then
            do j = ny, 1, -1
                flux_down(:, j) = bcs_ht(1:nlines)*tau(:, ny)/tau(:, j)                     ! downward flux
                flux_up(:, j) = bcs_hb(1:nlines)*tau(:, j)                                    ! upward flux
            end do

        end if

        return
    end subroutine IR_RTE1_OnlyLiquid_UpwardFirst

    !########################################################################
    !########################################################################
    subroutine IR_RTE1_Incremental(nlines, ny, fdmi, a_source, b, tau, bcs_hb, bcs_ht, flux_down, flux_up)
        integer(wi), intent(in) :: nlines, ny
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(inout) :: a_source(nlines, ny)     ! input as bulk absorption coefficent, output as source
        real(wp), intent(inout) :: b(nlines, ny)            ! input as emission function, output as upward flux, if flux is to be return
        real(wp), intent(inout) :: tau(nlines, ny)          ! optical depth and transmission functions
        real(wp), intent(inout) :: bcs_hb(nlines)           ! flux boundary condition at the bottom
        real(wp), intent(inout) :: bcs_ht(nlines)           ! flux boundary condition at the top
        real(wp), intent(inout) :: flux_down(nlines, ny)    ! flux_down for intermediate calculations and net flux as output
        real(wp), intent(out), optional :: flux_up(nlines, ny)

        ! -----------------------------------------------------------------------
        integer(wi) j
        real(wp) dummy
        real(wp), pointer :: p_wrk2d_1(:) => null()
        real(wp), pointer :: p_wrk2d_2(:) => null()

        ! #######################################################################
        p_wrk2d_1(1:nlines) => wrk2d(1:nlines, 1)
        p_wrk2d_2(1:nlines) => wrk2d(1:nlines, 2)

        ! ###################################################################
        ! absorption coefficient; divide by mean direction
        dummy = 1.0_wp/mu
        a_source = a_source*dummy

        ! emission function
        bcs_hb(1:nlines) = b(1:nlines, 1)   ! save for calculation of surface flux
        b = b*a_source                      ! absorption coefficient times emission function

        ! transmission function I_{j-1,j}  = exp(-tau(z_{j-1}, z_j)/\mu)
        tau(:, 1) = 0.0_wp                                    ! boundary condition
        ! call FDM_Int1_Solve(nxz, fdmi(BCS_MIN), fdmi(BCS_MIN)%rhs, a_source, tau, wrk2d)
        ! call Int_Trapezoidal_f(a_source, z%nodes(:), tau, BCS_MIN)
        call Int_Simpson_Biased_f(a_source, z%nodes(:), tau, BCS_MIN)
        do j = ny, 2, -1
            tau(:, j) = exp(tau(:, j - 1) - tau(:, j))
        end do
        tau(:, 1) = 1.0_wp        ! this value is not used

        ! ###################################################################
        ! downward flux; positive going down
        j = ny
        flux_down(1:nlines, j) = bcs_ht(1:nlines)

        do j = ny - 1, 1, -1
            ! Integral contribution from emission function using a trapezoidal rule
            p_wrk2d_2 = b(:, j + 1)
            p_wrk2d_1 = b(:, j)/tau(:, j + 1)
            p_wrk2d_1 = 0.5_wp*(p_wrk2d_1 + p_wrk2d_2)*(z%nodes(j + 1) - z%nodes(j))

            flux_down(:, j) = tau(:, j + 1)*(flux_down(:, j + 1) + p_wrk2d_1)
        end do

        ! ###################################################################
        ! bottom boundary condition; calculate upward flux at the bottom
        bcs_hb(1:nlines) = epsilon*bcs_hb(1:nlines) + (1.0_wp - epsilon)*flux_down(:, 1)

        ! ###################################################################
        ! upward flux and source
        j = 1
        a_source(:, j) = a_source(:, j)*(bcs_hb(1:nlines) + flux_down(:, j)) - 2.0_wp*b(:, j)

        if (present(flux_up)) then                                          ! Save fluxes
            flux_up(:, j) = bcs_hb(1:nlines)                                ! upward flux

            do j = 2, ny
                ! Integral contribution from emission function using a trapezoidal rule
                p_wrk2d_1 = b(:, j - 1)
                p_wrk2d_2 = b(:, j)/tau(:, j)
                p_wrk2d_1 = 0.5_wp*(p_wrk2d_1 + p_wrk2d_2)*(z%nodes(j) - z%nodes(j - 1))

                bcs_hb = tau(:, j)*(bcs_hb(1:nlines) + p_wrk2d_1)
                a_source(:, j) = a_source(:, j)*(bcs_hb(1:nlines) + flux_down(:, j)) - 2.0_wp*b(:, j)

                flux_up(:, j) = bcs_hb(1:nlines)
            end do

        else                                                                ! Do not save fluxes
            do j = 2, ny
                ! Integral contribution from emission function using a trapezoidal rule
                p_wrk2d_1 = b(:, j - 1)
                p_wrk2d_2 = b(:, j)/tau(:, j)
                p_wrk2d_1 = 0.5_wp*(p_wrk2d_1 + p_wrk2d_2)*(z%nodes(j) - z%nodes(j - 1))

                bcs_hb(1:nlines) = tau(:, j)*(bcs_hb(1:nlines) + p_wrk2d_1(1:nlines))
                a_source(:, j) = a_source(:, j)*(bcs_hb(1:nlines) + flux_down(:, j)) - 2.0_wp*b(:, j)
            end do

        end if

        ! ###################################################################
        nullify (p_wrk2d_1, p_wrk2d_2)

        return
    end subroutine IR_RTE1_Incremental

    !########################################################################
    !########################################################################
    subroutine IR_RTE1_Local(nlines, ny, fdmi, a_source, b, tau, bcs_hb, bcs_ht, flux_down, tmp2, flux_up)
        integer(wi), intent(in) :: nlines, ny
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(inout) :: a_source(nlines, ny)     ! input as bulk absorption coefficent, output as source
        real(wp), intent(inout) :: b(nlines, ny)            ! input as emission function, output as upward flux, if flux is to be return
        real(wp), intent(inout) :: tau(nlines, ny)          ! optical depth and transmission functions
        real(wp), intent(inout) :: bcs_hb(nlines)           ! flux boundary condition at the bottom
        real(wp), intent(inout) :: bcs_ht(nlines)           ! flux boundary condition at the top
        real(wp), intent(inout) :: flux_down(nlines, ny)    ! flux_down for intermediate calculations and net flux as output
        real(wp), intent(inout) :: tmp2(nlines, ny)
        real(wp), intent(out), optional :: flux_up(nlines, ny)

        ! -----------------------------------------------------------------------
        integer(wi) j, k
        real(wp) dummy
        real(wp), pointer :: p_flux_1(:) => null()
        real(wp), pointer :: p_flux_2(:) => null()

        ! #######################################################################
        p_flux_1(1:nlines) => wrk2d(1:nlines, 1)
        p_flux_2(1:nlines) => wrk2d(1:nlines, 2)

        ! ###################################################################
        ! absorption coefficient; divide by mean direction
        dummy = 1.0_wp/mu
        a_source = a_source*dummy

        ! emission function
        bcs_hb(1:nlines) = b(1:nlines, 1)   ! save for calculation of surface flux
        b = b*a_source                      ! absorption coefficient times emission function

        ! transmission function I_{j-1,j} = exp(-tau(z_{j-1}, z_j)/\mu)
        tau(:, 1) = 0.0_wp                                    ! boundary condition
        ! call FDM_Int1_Solve(nxz, fdmi(BCS_MIN), fdmi(BCS_MIN)%rhs, a_source, tau, wrk2d)
        ! call Int_Trapezoidal_f(a_source, z%nodes(:), tau, BCS_MIN)
        call Int_Simpson_Biased_f(a_source, z%nodes(:), tau, BCS_MIN)
        do j = ny, 2, -1
            tau(:, j) = exp(tau(:, j - 1) - tau(:, j))
        end do
        tau(:, 1) = 1.0_wp        ! this value is not used

        ! ###################################################################
        ! downward flux; positive going down
        j = ny
        p_flux_1(1:nlines) = bcs_ht(1:nlines)
        flux_down(:, j) = p_flux_1

        do j = ny - 1, 1, -1
            p_flux_1 = p_flux_1*tau(:, j + 1)                     ! accumulate transmission

            p_flux_2 = 1.0_wp                                       ! calculate emission; p_flux_2 is aux array in the next loop
            tmp2(:, j) = b(:, j)
            do k = j + 1, ny
                p_flux_2 = p_flux_2*tau(:, k)
                tmp2(:, k) = b(:, k)*p_flux_2
            end do
            call Int_Simpson_v(tmp2(:, j:ny), z%nodes(j:ny), flux_down(:, j))
            ! call Int_Trapezoidal_v(tmp2(:, j:ny), z%nodes(j:ny), flux_down(:, j))
            flux_down(:, j) = flux_down(:, j) + p_flux_1
        end do

        ! ###################################################################
        ! bottom boundary condition; calculate upward flux at the bottom
        bcs_hb(1:nlines) = epsilon*bcs_hb(1:nlines) + (1.0_wp - epsilon)*flux_down(:, 1)

        ! ###################################################################
        ! upward flux and net terms
        j = 1
        p_flux_1(1:nlines) = bcs_hb(1:nlines)
        a_source(:, j) = a_source(:, j)*(p_flux_1 + flux_down(:, j)) - 2.0_wp*b(:, j)

        if (present(flux_up)) then                                  ! I need additional space to store fluxes
            flux_up(:, j) = p_flux_1                                ! upward flux

            do j = 2, ny
                p_flux_1 = p_flux_1*tau(:, j)

                p_flux_2 = 1.0_wp       ! used as auxiliary array in the next loop
                tmp2(:, j) = b(:, j)
                do k = j - 1, 1, -1
                    p_flux_2 = p_flux_2*tau(:, k + 1)
                    tmp2(:, k) = b(:, k)*p_flux_2
                end do
                call Int_Simpson_v(tmp2(:, 1:j), z%nodes(1:j), p_flux_2)
                ! call Int_Trapezoidal_v(tmp2(:, 1:j), z%nodes(1:j), p_flux_2)
                a_source(:, j) = a_source(:, j)*(p_flux_2 + p_flux_1 + flux_down(:, j)) - 2.0_wp*b(:, j)

                flux_up(:, j) = p_flux_2 + p_flux_1                 ! upward flux
            end do

        else                ! we only need the source term
            do j = 2, ny
                p_flux_1 = p_flux_1*tau(:, j)                     ! accumulate transmission

                p_flux_2 = 1.0_wp                                   ! calculate emission; p_flux_2 is aux array in the next loop
                tmp2(:, j) = b(:, j)
                do k = j - 1, 1, -1
                    p_flux_2 = p_flux_2*tau(:, k + 1)
                    tmp2(:, k) = b(:, k)*p_flux_2
                end do
                call Int_Simpson_v(tmp2(:, 1:j), z%nodes(1:j), p_flux_2)
                ! call Int_Trapezoidal_v(tmp2(:, 1:j), z%nodes(1:j), p_flux_2)
                a_source(:, j) = a_source(:, j)*(p_flux_2 + p_flux_1 + flux_down(:, j)) - 2.0_wp*b(:, j)

            end do

        end if

        ! #######################################################################
        nullify (p_flux_1, p_flux_2)

        return
    end subroutine IR_RTE1_Local

    !########################################################################
    !########################################################################
    subroutine IR_RTE1_Global(nlines, ny, fdmi, a_source, b, tau, bcs_hb, bcs_ht, flux_down, flux_up)
        integer(wi), intent(in) :: nlines, ny
        type(fdm_integral_dt), intent(in) :: fdmi(2)
        real(wp), intent(inout) :: a_source(nlines, ny)     ! input as bulk absorption coefficent, output as source
        real(wp), intent(inout) :: b(nlines, ny)            ! input as emission function, output as upward flux, if flux is to be return
        real(wp), intent(inout) :: tau(nlines, ny)          ! optical depth and transmission functions
        real(wp), intent(inout) :: bcs_hb(nlines)           ! flux boundary condition at the bottom
        real(wp), intent(inout) :: bcs_ht(nlines)           ! flux boundary condition at the top
        real(wp), intent(inout) :: flux_down(nlines, ny)    ! flux_down for intermediate calculations and net flux as output
        real(wp), intent(inout) :: flux_up(nlines, ny)

        ! -----------------------------------------------------------------------
        integer(wi) j
        real(wp) dummy

        ! ###################################################################
        ! absorption coefficient; divide by mean direction
        dummy = 1.0_wp/mu
        a_source = a_source*dummy

        ! emission function
        bcs_hb(1:nlines) = b(1:nlines, 1)   ! save for calculation of surface flux
        b = b*a_source                      ! absorption coefficient times emission function

        ! ###################################################################
        ! downward flux; positive going down

        ! transmission function I_j = exp(-tau(z_j, zmax)/\mu)
        tau(:, ny) = 0.0_wp                                   ! boundary condition
        ! call FDM_Int1_Solve(nlines, fdmi(BCS_MAX), fdmi(BCS_MIN)%rhs, a_source, tau, wrk2d)         ! recall this gives the negative of the integral
        ! call Int_Trapezoidal_f(a_source, z%nodes(:), tau, BCS_MAX)
        call Int_Simpson_Biased_f(a_source, z%nodes(:), tau, BCS_MAX)
        do j = ny, 1, -1
            tau(:, j) = exp(-tau(:, j))
        end do
        !  tau = dexp(tau)         seg-fault; need ulimit -u unlimited

        flux_down = b/tau
        ! call Int_Trapezoidal_Increments_InPlace(flux_down, z%nodes(:), BCS_MAX)                   ! Calculate I_j = int_{x_{j}}^{x_{j+1}}
        call Int_Simpson_Biased_Increments_InPlace(flux_down, z%nodes(:), wrk2d(:, 1), BCS_MAX)   ! Calculate I_j = int_{x_{j}}^{x_{j+1}}
        j = ny
        wrk2d(1:nlines, 1) = 0.0_wp                                         ! accumulate emission; using wrk2d as aux array
        flux_down(:, j) = tau(:, j)*(bcs_ht(1:nlines) + wrk2d(1:nlines, 1))
        do j = ny - 1, 1, -1
            wrk2d(1:nlines, 1) = wrk2d(1:nlines, 1) + flux_down(:, j)            ! accumulate emission
            flux_down(:, j) = tau(:, j)*(bcs_ht(1:nlines) + wrk2d(1:nlines, 1))
        end do

        ! ###################################################################
        ! bottom boundary condition; calculate upward flux at the bottom
        bcs_hb(1:nlines) = epsilon*bcs_hb(1:nlines) + (1.0_wp - epsilon)*flux_down(:, 1)

        ! ###################################################################
        ! upward flux and net terms

        ! transmission function I_j = exp(-tau(zmin, z)/\mu)
        tau(:, 1) = 0.0_wp                                                ! boundary condition
        ! call FDM_Int1_Solve(nlines, fdmi(BCS_MIN), fdmi(BCS_MIN)%rhs, a_source, tau, wrk2d)
        ! call Int_Trapezoidal_f(a_source, z%nodes(:), tau, BCS_MIN)
        call Int_Simpson_Biased_f(a_source, z%nodes(:), tau, BCS_MIN)
        do j = 1, ny
            tau(:, j) = exp(-tau(:, j))
        end do

        flux_up = b/tau
        ! call Int_Trapezoidal_Increments_InPlace(flux_up, z%nodes(:), BCS_MIN)                     ! Calculate I_j = int_{x_{j-1}}^{x_{j}}
        call Int_Simpson_Biased_Increments_InPlace(flux_up, z%nodes(:), wrk2d(:, 1), BCS_MIN)     ! Calculate I_j = int_{x_{j-1}}^{x_{j}}
        j = 1
        wrk2d(1:nlines, 1) = 0.0_wp                                         ! accumulate emission; using wrk2d as aux array
        flux_up(:, j) = tau(:, j)*(bcs_hb(1:nlines) + wrk2d(1:nlines, 1))
        !
        a_source(:, j) = a_source(:, j)*(flux_up(:, j) + flux_down(:, j)) - 2.0_wp*b(:, j)
        do j = 2, ny
            wrk2d(1:nlines, 1) = wrk2d(1:nlines, 1) + flux_up(:, j)         ! accumulate emission
            flux_up(:, j) = tau(:, j)*(bcs_hb(1:nlines) + wrk2d(1:nlines, 1))
            !
            a_source(:, j) = a_source(:, j)*(flux_up(:, j) + flux_down(:, j)) - 2.0_wp*b(:, j)
        end do

        return
    end subroutine IR_RTE1_Global

end module
