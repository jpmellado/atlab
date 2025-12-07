!########################################################################
!# Sources (processes) in the evolution equations.
!########################################################################
module TLab_Sources
    use TLab_Constants, only: wp, wi, small_wp
    use TLab_Memory, only: imax, jmax, kmax, isize_field, inb_scal, inb_scal_array
    use NavierStokes, only: nse_eqns, DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC
    use Thermo_Anelastic, only: ribackground, Thermo_Anelastic_AddBuoyancy, Thermo_Anelastic_Weight_Add
    use Thermo_Anelastic, only: Thermo_Anelastic_AddBuoyancy_PerVolume
    use Gravity, only: gravityProps, Gravity_AddSource
    use Rotation, only: coriolisProps, Rotation_AddCoriolis, Rotation_AddCoriolis_PerVolume
    use Microphysics, only: sedimentationProps, Microphysics_Sedimentation_Z
    use Microphysics, only: evaporationProps, Microphysics_Evaporation
    use Radiation, only: infraredProps, Radiation_Infrared_Z
    use SpecialForcing
    ! use LargeScaleForcing
    implicit none
    private

    public :: TLab_Sources_Flow
    public :: TLab_Sources_Scal

contains
! #######################################################################
! #######################################################################
    subroutine TLab_Sources_Flow(q, s, hq, tmp1)
        use TLab_Time, only: rtime
        real(wp), intent(in) :: q(isize_field, *), s(isize_field, *)
        real(wp), intent(out) :: hq(isize_field, *)
        real(wp), intent(inout) :: tmp1(isize_field)

        ! -----------------------------------------------------------------------
        integer iq

        ! #######################################################################
        select case (nse_eqns)
        case (DNS_EQNS_BOUSSINESQ)
            call Rotation_AddCoriolis(coriolisProps, imax, jmax, kmax, q, hq)
        case (DNS_EQNS_ANELASTIC)
            call Rotation_AddCoriolis_PerVolume(coriolisProps, imax, jmax, kmax, q, hq)

        end select

        do iq = 1, 3

            if (gravityProps%active(iq)) then
                select case (nse_eqns)
                case (DNS_EQNS_BOUSSINESQ)
                    call Gravity_AddSource(gravityProps, imax, jmax, kmax, s, hq(:, iq), gravityProps%vector(iq))

                case (DNS_EQNS_ANELASTIC)
                    call Thermo_Anelastic_AddBuoyancy_PerVolume(imax, jmax, kmax, s, hq(:, iq), gravityProps%vector(iq))

                end select

            end if

            ! -----------------------------------------------------------------------
            ! Subsidence
            ! Implemented directly in burgers

            ! -----------------------------------------------------------------------
            if (forcingProps%active(iq)) then
                call SpecialForcing_Source(forcingProps, imax, jmax, kmax, iq, rtime, q(:, iq), hq(:, iq), tmp1)

                hq(:, iq) = hq(:, iq) + tmp1(:)

            end if

        end do

        return
    end subroutine TLab_Sources_Flow

    ! #######################################################################
    ! #######################################################################
    subroutine TLab_Sources_Scal(s, hs, tmp1, tmp2, tmp3, tmp4)
        real(wp), intent(in) :: s(isize_field, *)
        real(wp), intent(out) :: hs(isize_field, *)
        real(wp), intent(inout) :: tmp1(isize_field), tmp2(isize_field), tmp3(isize_field), tmp4(isize_field)

        ! -----------------------------------------------------------------------
        integer is

        ! #######################################################################
        do is = 1, inb_scal

            if (infraredProps%active(is)) then
                call Radiation_Infrared_Z(infraredProps, imax, jmax, kmax, s, tmp1, tmp2, tmp3, tmp4)
                hs(:, is) = hs(:, is) + tmp1(:)
            end if

            if (sedimentationProps%active(is)) then
                call Microphysics_Sedimentation_Z(sedimentationProps, imax, jmax, kmax, is, s, tmp1, tmp2)
                hs(:, is) = hs(:, is) + tmp1(:)
            end if

            if (evaporationProps%active(is)) then
                call Microphysics_Evaporation(evaporationProps, imax, jmax, kmax, is, s, tmp1)
                hs(:, is) = hs(:, is) + tmp1(:)
            end if

            ! -----------------------------------------------------------------------
            ! Subsidence
            ! Implemented directly in burgers

        end do

        return
    end subroutine TLab_Sources_Scal

end module TLab_Sources
