program vThermo
    use TLab_Constants, only: wp, wi, fmt_r
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use Thermodynamics
    use Thermo_Anelastic

    implicit none

    real(wp) h(1), qt(1), ql(1)
    real(wp) e(1), p(1), ps(1), t(1), r(1), qs(1), qv(1)
    real(wp) dummy(1), dqldqt(1), ep(1), theta(1), theta_v(1), theta_l(1), theta_e(1), Td(1), wrk(1)
    real(wp) heat1(1), heat2(1), cp1(1), cp2(1), alpha(1), as(1), bs(1)
    real(wp) h1(1), s(3), r1(1)
    integer(wi) iopt

! ###################################################################
    call TLab_Start()

    imixture = MIXT_TYPE_AIRWATER
    nondimensional = .false.
    call Thermo_Initialize()
    ep = 0.0_wp
    dsmooth = 1.0_wp
    scaleheightinv = 1.0_wp

    write (*, *) 'Case p-t (1) or d-e (2) or p-h (3)?'
    read (*, *) iopt

    if (iopt == 1) then
        write (*, *) 'temperature (C) ?'
        read (*, *) t
        t = t + 273.15
        write (*, *) 'pressure (hPa) ?'
        read (*, *) p
        p = p*100.0_wp

    else if (iopt == 2) then
        write (*, *) 'density ?'
        read (*, *) r
        write (*, *) 'energy ?'
        read (*, *) e

    else if (iopt == 3) then
        write (*, *) 'enthalpy (kJ/kg)?'
        read (*, *) h
        write (*, *) 'pressure (hPa) ?'
        read (*, *) p
        p = p*100.0_wp

    end if

    allocate (pbackground(1), epbackground(1), rbackground(1))
    epbackground(1) = ep(1)
    pbackground(1) = p(1)
    rbackground(1) = r(1)

    write (*, *) 'water specific humidity (g/kg) ?'
    read (*, *) qt
    qt = qt*0.001_wp

! ###################################################################
    if (iopt == 1) then
        call Thermo_Psat_Polynomial(1, t, ps)
        qs = 1.0_wp/(p/ps - 1.0_wp)*rd_ov_rv
        qs = qs/(1.0_wp + qs)
        if (qt(1) > qs(1)) then
            qv = qs*(1 - qt)
            ql = qt - qv
        else
            qv = qt
            ql = 0.0_wp
        end if
        ! z1(1) = qt(1)
        ! z1(2) = ql(1)
        ! call THERMO_CALORIC_ENTHALPY(1, z1, t, h)
        ! call THERMO_CALORIC_ENERGY(1, z1, t, e)
        ! call THERMO_THERMAL_DENSITY(1, z1, p, t, r)

        s(1) = h(1)
        s(2) = qt(1)
        s(3) = ql(1)
        call Thermo_Anelastic_Theta(1, 1, 1, s, theta, wrk)
        call Thermo_Anelastic_ThetaV(1, 1, 1, s, theta_v, wrk)
        call Thermo_Anelastic_ThetaL(1, 1, 1, s, theta_l, wrk)
        call Thermo_Anelastic_ThetaE(1, 1, 1, s, theta_e, wrk)
        ! call Thermo_Anelastic_DEWPOINT(1, 1, 1, s, dummy, Td, wrk)

    else if (iopt == 2) then
        z1(1) = qt(1)
        call THERMO_CALORIC_TEMPERATURE(1, z1, e, r, T, dqldqt)
        ql = z1(2)
        qv = qt - ql
        qs = qv ! initial condition for next routine
        call THERMO_THERMAL_PRESSURE(1, z1, r, t, p)
        call Thermo_Psat_Polynomial(1, t, ps)
        qs = 1.0_wp/(p/ps - 1.0_wp)*rd_ov_rv
        qs = qs/(1.0_wp + qs)
        call THERMO_CALORIC_ENTHALPY(1, z1, t, h)

        s(1) = h(1); s(2:3) = z1(1:2)
        call Thermo_Anelastic_THETA_V(1, 1, 1, s, theta_v, wrk)
        call Thermo_Anelastic_THETA_L(1, 1, 1, s, theta_l, wrk)
        call Thermo_Anelastic_THETA_E(1, 1, 1, s, theta_e, wrk)
        call Thermo_Anelastic_DEWPOINT(1, 1, 1, s, dummy, Td, wrk)

    else if (iopt == 3) then
        ! h = h/TREF/1.007
        z1(1) = qt(1)
        call Thermo_Anelastic_PH(1, 1, 1, z1, h)
        s(1) = h(1); s(2:3) = z1(1:2)
        call Thermo_Anelastic_TEMPERATURE(1, 1, 1, s, T)
        ! CALL THERMO_AIRWATER_PH_RE(1, z1, p, h, T)
        ql(1) = z1(2)
        qv = qt - ql

        call Thermo_Psat_Polynomial(1, T, ps)
        qs = 1.0_wp/(p/ps - 1.0_wp)*rd_ov_rv
        qs = qs/(1.0_wp + qs)
        call THERMO_THERMAL_DENSITY(1, z1, p, T, r)
        call THERMO_CALORIC_ENERGY(1, z1, T, e)

        call Thermo_Anelastic_THETA_V(1, 1, 1, s, theta_v, wrk)
        call Thermo_Anelastic_THETA_L(1, 1, 1, s, theta_l, wrk)
        call Thermo_Anelastic_THETA_E(1, 1, 1, s, theta_e, wrk)
        call Thermo_Anelastic_DEWPOINT(1, 1, 1, s, dummy, Td, wrk)

! check
        call Thermo_Anelastic_DENSITY(1, 1, 1, s, r1, wrk)
!     r2 = p/(T*(1- qt +qv/rd_ov_rv ) )
        call THERMO_CALORIC_ENTHALPY(1, z1, T, h1)

    end if

    write (*, fmt_r) 'Saturation specific humidity (g/kg):', qs*1.0e3_wp
    write (*, fmt_r) 'Vapor specific humidity (g/kg).....:', qv*1.0e3_wp
    write (*, fmt_r) 'Liquid specific humidity (g/kg)....:', ql*1.0e3_wp
    write (*, fmt_r) 'Density ...........................:', r
    write (*, fmt_r) 'Pressure (hPa) ....................:', p/100.0_wp
    write (*, fmt_r) 'Saturation pressure (hPa) .........:', ps/100.0_wp
    write (*, fmt_r) 'Temperature (K) ...................:', t !- 273.15
    write (*, fmt_r) 'Dewpoint temperature (K) ..........:', Td
    write (*, fmt_r) 'Specific heat capacity (J/kg/T)....:', Cd + qt*Cdv + ql*Cvl
    write (*, fmt_r) 'Specific energy  (J/kg)............:', e
    write (*, fmt_r) 'Specific enthalpy (J/kg)...........:', h
    write (*, fmt_r) 'Reference latent heat (J/kg) ......:', -THERMO_AI(6, 1, 3)
    write (*, fmt_r) 'Latent heat (J/kg) ................:', (-Cl - t*Lvl)
    write (*, fmt_r) 'Virtual potential T (K) ...........:', theta_v
    write (*, fmt_r) 'Liquid-water potential T (K) ......:', theta_l
    write (*, fmt_r) 'Equivalent potential T (K) ........:', theta_e
    if (iopt == 3) then
        write (*, fmt_r) 'Density ...........................:', r1
        write (*, fmt_r) 'Specific enthalpy .................:', h1
    end if

    stop

end program vThermo
