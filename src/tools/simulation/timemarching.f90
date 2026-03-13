#include "tlab_error.h"

module TimeMarching
    use TLab_Constants, only: wp, wi, big_wp
    use TLab_Constants, only: efile, wfile
    use TLab_WorkFlow, only: flow_on, scal_on
    use TLab_Memory, only: imax, jmax, kmax, isize_field
    use TLab_Time, only: rtime
    use TLab_Memory, only: inb_flow, inb_scal
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop
    ! use PARTICLE_VARS
#ifdef USE_MPI
    use mpi_f08
    use TLabMPI_VARS
#endif
    use NavierStokes, only: nse_eqns, DNS_EQNS_COMPRESSIBLE, DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC
    use NavierStokes, only: visc, schmidt, prandtl
    use DNS_Arrays
    use RungeKutta

    implicit none
    private

    public :: TMarch_Initialize
    public :: TMarch_RungeKutta
    public :: TMarch_Courant

    class(rungekutta_lowstorage_dt), allocatable, protected :: TMarch

    real(wp), public :: dtime                       ! time step
    logical, public :: use_variable_timestep = .true.

    ! real(wp) :: dtrkm                       ! time step of each substep
    ! real(wp) :: dte                         ! explicit scheme coefficient at each RK substep; to be removed from public
    logical :: remove_divergence            ! Remove residual divergence every time step

    ! ! -------------------------------------------------------------------
    ! integer(wi) :: rkm_mode                     ! Type of Runge-Kutta scheme
    ! integer, parameter :: RKM_EXP3 = 3
    ! integer, parameter :: RKM_EXP4 = 4
    ! integer, parameter :: RKM_IMP3_DIFFUSION = 5

    ! integer(wi) :: rkm_endstep                  ! number of substeps
    ! integer(wi) :: rkm_substep                  ! substep counter

    real(wp) :: cfla, cfld, cflr                ! CFL numbers

    ! real(wp) kdt(5), kco(4), ktime(6)           ! explicit scheme coefficients
    ! real(wp) kex(3), kim(3)                     ! implicit scheme coefficients
    ! real(wp) :: etime                           ! time at each RK substep; rtime is time at each iteration, not at each RK substep

    real(wp) schmidtfactor, dx2i
    real(wp) dummy

    type :: ds_dt
        real(wp), allocatable :: one_ov_ds1(:)
        real(wp), allocatable :: one_ov_ds2(:)
    end type
    type(ds_dt) :: ds(3)

contains

    ! ###################################################################
    ! ###################################################################
    subroutine TMarch_Initialize(inifile)
        use TLab_Memory, only: TLab_Allocate_Real
        use TLab_Arrays, only: wrk1d
        use TLab_Grid, only: globalGrid, x, y, z, xSubgrid, ySubgrid, zSubgrid
        use FDM, only: fdm_der1_X, fdm_der1_Y, fdm_der1_Z

        character*(*) inifile

        ! -------------------------------------------------------------------
        character(len=32) bakfile, block, lstr
        character(len=128) eStr
        character(len=512) sRes
        integer ig
        integer(wi) i, j, k

        ! ###################################################################
        ! read
        bakfile = trim(adjustl(inifile))//'.bak'

        block = 'Time'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        call TLab_Write_ASCII(bakfile, '#')
        call TLab_Write_ASCII(bakfile, '#['//trim(adjustl(block))//']')
        call TLab_Write_ASCII(bakfile, '#Scheme=<RungeKuttaExplicit3/RungeKuttaExplicit4>')
        call TLab_Write_ASCII(bakfile, '#TimeStep=<value>')
        call TLab_Write_ASCII(bakfile, '#MaxCFL=<value>')
        call TLab_Write_ASCII(bakfile, '#MaxDiffusiveCFL=<value>')
        call TLab_Write_ASCII(bakfile, '#MaxReactiveCFL=<value>')
        call TLab_Write_ASCII(bakfile, '#RemoveDivergence=<none/remove>')

        call ScanFile_Char(bakfile, inifile, block, 'Scheme', 'dummy', sRes)
        select case (trim(adjustl(sRes)))
        case ('rungekuttaexplicit3')
            allocate (rk3_dt :: TMarch)
            ! rkm_mode = RKM_EXP3
            lstr = '0.6'; 
        case ('rungekuttaexplicit4')
            allocate (rk45_dt :: TMarch)
            ! rkm_mode = RKM_EXP4
            lstr = '1.2'; 
        case default
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong Scheme option.')
            call TLab_Stop(DNS_ERROR_RKORDER)
        end select

        ! Default cfla value set in lstr while reading Scheme
        call ScanFile_Real(bakfile, inifile, block, 'MaxCFL', trim(adjustl(lstr)), cfla)
        write (lstr, *) 0.25_wp*cfla ! Default value for diffusive CFL
        call ScanFile_Real(bakfile, inifile, block, 'MaxDiffusiveCFL', trim(adjustl(lstr)), cfld)
        write (lstr, *) 0.5_wp*cfla ! Default value for reactive CFL
        call ScanFile_Real(bakfile, inifile, block, 'MaxReactiveCFL', trim(adjustl(lstr)), cflr)

        call ScanFile_Char(bakfile, inifile, block, 'TimeStep', 'void', sRes)
        if (trim(adjustl(sRes)) /= 'void') then
            use_variable_timestep = .false.
            read (sRes, *) dtime

            call ScanFile_Char(bakfile, inifile, block, 'MaxCFL', 'void', sRes)
            if (trim(adjustl(sRes)) /= 'void') then
                call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Cannot impose both time step and max CFL.')
                call TLab_Stop(DNS_ERROR_OPTION)
            end if

        end if

        call ScanFile_Char(bakfile, inifile, block, 'RemoveDivergence', 'yes', sRes)
        if (trim(adjustl(sRes)) == 'no') then; remove_divergence = .false.
        else if (trim(adjustl(sRes)) == 'yes') then; remove_divergence = .true.
        else
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Wrong RemoveDivergence option.')
            call TLab_Stop(DNS_ERROR_OPTION)
        end if

        ! ###################################################################
        ! initialize
        call TMarch%initialize()

        ! ! RK coefficients
        ! select case (rkm_mode)
        ! case (RKM_EXP3)             ! Runge-Kutta explicit 3th order from Williamson 1980
        !     rkm_endstep = 3

        !     kdt(1:3) = [1.0_wp/3.0_wp, 15.0_wp/16.0_wp, 8.0_wp/15.0_wp]
        !     ktime(1:3) = [0.0_wp, 1.0_wp/3.0_wp, 3.0_wp/4.0_wp]
        !     kco(1:2) = [-5.0_wp/9.0_wp, -153.0_wp/128.0_wp]

        ! case (RKM_EXP4)             ! Runge-Kutta explicit 4th order 5 stages from Carpenter & Kennedy 1994
        !     rkm_endstep = 5

            ! kdt(1) = 1432997174477.0_wp/9575080441755.0_wp
            ! kdt(2) = 5161836677717.0_wp/13612068292357.0_wp
            ! kdt(3) = 1720146321549.0_wp/2090206949498.0_wp
            ! kdt(4) = 3134564353537.0_wp/4481467310338.0_wp
            ! kdt(5) = 2277821191437.0_wp/14882151754819.0_wp

        ! ktime(1) = 0.0_wp
        ! ktime(2) = kdt(1)
        ! ktime(3) = 2526269341429.0_wp/6820363962896.0_wp
        ! ktime(4) = 2006345519317.0_wp/3224310063776.0_wp
        ! ktime(5) = 2802321613138.0_wp/2924317926251.0_wp
        ! ktime(6) = 1.0_wp

        !     kco(1) = -567301805773.0_wp/1357537059087.0_wp
        !     kco(2) = -2404267990393.0_wp/2016746695238.0_wp
        !     kco(3) = -3550918686646.0_wp/2091501179385.0_wp
        !     kco(4) = -1275806237668.0_wp/842570457699.0_wp

        ! end select

        ! ###################################################################
        ! Memory management
        call TLab_Allocate_Real(__FILE__, hq, [isize_field, inb_flow], 'flow-rhs')
        call TLab_Allocate_Real(__FILE__, hs, [isize_field, inb_scal], 'scal-rhs')

        p_hq(1:imax, 1:jmax, 1:kmax, 1:inb_flow) => hq(1:imax*jmax*kmax*inb_flow, 1)
        p_hs(1:imax, 1:jmax, 1:kmax, 1:inb_scal) => hs(1:imax*jmax*kmax*inb_scal, 1)

        pxy_hq(1:imax*jmax, 1:kmax, 1:inb_flow) => hq(1:imax*jmax*kmax*inb_flow, 1)
        pxy_hs(1:imax*jmax, 1:kmax, 1:inb_scal) => hs(1:imax*jmax*kmax*inb_scal, 1)

        ! ###################################################################
        ! maximum diffusivities for TMarch_Courant
        schmidtfactor = 1.0_wp
        dummy = 1.0_wp/prandtl
        schmidtfactor = max(schmidtfactor, dummy)
        dummy = 1.0_wp/minval(schmidt(1:inb_scal))
        schmidtfactor = max(schmidtfactor, dummy)

        ! ###################################################################
        ! Calculate inverse of Jacobian (inverse of grid spacing)
        do ig = 1, 3
            allocate (ds(ig)%one_ov_ds1(globalGrid%axes(ig)%size))
            if (globalGrid%axes(ig)%uniform) then
                if (globalGrid%axes(ig)%size > 1) then
                    ds(ig)%one_ov_ds1(:) = 1.0_wp/(globalGrid%axes(ig)%nodes(2) - globalGrid%axes(ig)%nodes(1))
                else
                    ds(ig)%one_ov_ds1(:) = 1.0_wp   ! 2d case
                end if
            else
                wrk1d(1:globalGrid%axes(ig)%size, 1) = [(real(i - 1, wp), i=1, globalGrid%axes(ig)%size)]
                select case (ig)
                case (1)
                    call fdm_der1_X%compute(1, wrk1d(:, 1), ds(ig)%one_ov_ds1(:))
                case (2)
                    call fdm_der1_Y%compute(1, wrk1d(:, 1), ds(ig)%one_ov_ds1(:))
                case (3)
                    call fdm_der1_Z%compute(1, wrk1d(:, 1), ds(ig)%one_ov_ds1(:))
                end select
            end if

            allocate (ds(ig)%one_ov_ds2(globalGrid%axes(ig)%size))
            ds(ig)%one_ov_ds2(:) = ds(ig)%one_ov_ds1(:)*ds(ig)%one_ov_ds1(:)

        end do

        ! Maximum of (1/dx^2 + 1/dy^2 + 1/dz^2) for TMarch_Courant
        dx2i = 0.0_wp
        do k = 1, kmax
            do j = 1, jmax
                do i = 1, imax
                    dummy = 0.0_wp
                    if (x%size > 1) dummy = dummy + ds(1)%one_ov_ds2(i + xSubgrid%offset)
                    if (y%size > 1) dummy = dummy + ds(2)%one_ov_ds2(j + ySubgrid%offset)
                    if (z%size > 1) dummy = dummy + ds(3)%one_ov_ds2(k + zSubgrid%offset)
                    dx2i = max(dx2i, dummy)
                end do
            end do
        end do

#ifdef USE_MPI
        call MPI_ALLREDUCE(dx2i, dummy, 1, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
        dx2i = dummy
#endif

        call TMarch_Courant()

        return
    end subroutine TMarch_Initialize

    ! ###################################################################
    ! ###################################################################
    subroutine TMarch_RungeKutta()
        use TLab_Arrays, only: q, s, txc
        use TLab_Sources, only: TLab_Sources_Flow, TLab_Sources_Scal, TLab_Sources_Scal_Implicit
        use TLab_Pointers_2D, only: pxy_q, pxy_s
        use DNS_Control, only: DNS_Limit_Bounds
        use Thermo_Anelastic, only: ribackground
        use Buffer

        ! -------------------------------------------------------------------
        integer rkm_substep
#ifdef USE_PROFILE
        integer(wi) t_srt, t_end, t_dif, idummy, PROC_CYCLES, MAX_CYCLES
        character*256 time_string
#endif

        !########################################################################
        ! Initialize arrays to zero for the explicit low-storage algorithm
        if (flow_on) hq = 0.0_wp
        if (scal_on) hs = 0.0_wp
        TMarch%time = rtime

        !########################################################################
        ! Loop over the sub-stages
        !########################################################################
        do rkm_substep = 1, TMarch%num_substep !rkm_endstep
#ifdef USE_PROFILE
            call system_clock(t_srt, PROC_CYCLES, MAX_CYCLES)
#endif
            TMarch%substep = rkm_substep

            ! dte = dtime*kdt(rkm_substep)
            ! etime = rtime + dtime*ktime(rkm_substep)
            ! dtrkm = dtime*(ktime(rkm_substep + 1) - ktime(rkm_substep))

            ! -------------------------------------------------------------------
            ! Explicit part
            ! call TLab_Sources_Flow(q, s, etime, hq, txc(:, 1))
            ! call TLab_Sources_Scal(s, hs, etime, txc(:, 1), txc(:, 2), txc(:, 3), txc(:, 4))
            call TLab_Sources_Flow(q, s, TMarch%time, hq, txc(:, 1))
            call TLab_Sources_Scal(s, hs, TMarch%time, txc(:, 1), txc(:, 2), txc(:, 3), txc(:, 4))

            if (bufferType == BUFFER_TYPE_NUDGE) call Buffer_Nudge()

            select case (nse_eqns)
            case (DNS_EQNS_BOUSSINESQ)
                ! call NSE_Boussinesq(dte, remove_divergence)
                call NSE_Boussinesq(TMarch%coef_a(TMarch%substep)*dtime, remove_divergence)
                call NSE_Boussinesq_BscFlow()
                call NSE_Boussinesq_BscScal()
                call TMarch%AdvanceSubstep_Boussinesq(pxy_q, pxy_hq, dtime)
                call TMarch%AdvanceSubstep_Boussinesq(pxy_s, pxy_hs, dtime)
                ! call TMarch_Substep_Boussinesq_Explicit()

            case (DNS_EQNS_ANELASTIC)
                ! call NSE_Anelastic_PerVolume(dte, remove_divergence)
                call NSE_Anelastic_PerVolume(TMarch%coef_a(TMarch%substep)*dtime, remove_divergence)
                call NSE_Anelastic_PerVolume_BscFlow()
                call NSE_Anelastic_PerVolume_BscScal()
                call TMarch%AdvanceSubstep_Anelastic(pxy_q, pxy_hq, dtime, ribackground)
                call TMarch%AdvanceSubstep_Anelastic(pxy_s, pxy_hs, dtime, ribackground)
                ! call TMarch_Substep_Anelastic_Explicit()

            case (DNS_EQNS_COMPRESSIBLE)

            end select

            ! Implicit part
            ! call TLab_Sources_Scal_Implicit(time_step=dtrkm, s=s)
            call TLab_Sources_Scal_Implicit(time_step=TMarch%coef_t(TMarch%substep)*dtime, s=s)
            ! Update of boundary condition still missing

            ! Limiters
            call DNS_Limit_Bounds()

            ! Diagnostics
            call TLab_Diagnostic(imax, jmax, kmax, s)

            ! TMarch%substep = TMarch%substep + 1
            TMarch%time = TMarch%time + TMarch%coef_t(TMarch%substep)*dtime

            ! -------------------------------------------------------------------
            ! Profiling data
#ifdef USE_PROFILE
            call system_clock(t_end, PROC_CYCLES, MAX_CYCLES)
            idummy = t_end - t_srt

#ifdef USE_MPI
            call MPI_REDUCE(idummy, t_dif, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)
            if (mpiGrid%rank == 0) then
         write (time_string, 999) mpiGrid%num_processors, xMpi%num_processors, zMpi%num_processors, rkm_substep, t_dif/1.0_wp/PROC_CYCLES/mpiGrid%num_processors
999             format(I5.5, ' (xMpi%num_processors X zMpi%num_processors:', I4.4, 'x', I4.4, 1x, ') RK-Substep', I1, ':', E13.5, 's')
                call TLab_Write_ASCII(lfile, time_string)
            end if
#else
            t_dif = idummy
            write (time_string, 999) rkm_substep, t_dif/1.0_wp/PROC_CYCLES/mpiGrid%num_processors
999         format('RK-Substep', I1, ':', E13.5, 's')
            call TLab_Write_ASCII(lfile, time_string)
#endif

#endif
        end do

        return
    end subroutine TMarch_RungeKutta

    ! !########################################################################
    ! !########################################################################
    ! ! Perform the time stepping
    ! ! Memory intensive, cache-optimized
    ! ! create size of memory blocks during initialization
    ! subroutine TMarch_Substep_Boussinesq_Explicit()
    !     use TLab_Pointers_2D, only: pxy_q, pxy_s

    !     real(wp) alpha
    !     integer iq, is, k

    !     ! #######################################################################
    !     if (rkm_substep < rkm_endstep) then
    !         alpha = kco(rkm_substep)

    !         do iq = 1, inb_flow
    !             do k = 1, kmax
    !                 pxy_q(:, k, iq) = pxy_q(:, k, iq) + dte*pxy_hq(:, k, iq)
    !                 pxy_hq(:, k, iq) = pxy_hq(:, k, iq)*alpha
    !             end do
    !         end do

    !         do iq = 1, inb_scal
    !             do k = 1, kmax
    !                 pxy_s(:, k, iq) = pxy_s(:, k, iq) + dte*pxy_hs(:, k, iq)
    !                 pxy_hs(:, k, iq) = pxy_hs(:, k, iq)*alpha
    !             end do
    !         end do

    !     else
    !         do is = 1, inb_flow
    !             do k = 1, kmax
    !                 pxy_q(:, k, is) = pxy_q(:, k, is) + dte*pxy_hq(:, k, is)
    !             end do
    !         end do

    !         do is = 1, inb_scal
    !             do k = 1, kmax
    !                 pxy_s(:, k, is) = pxy_s(:, k, is) + dte*pxy_hs(:, k, is)
    !             end do
    !         end do
    !     end if

    !     return
    ! end subroutine TMarch_Substep_Boussinesq_Explicit

    ! !########################################################################
    ! !########################################################################
    ! subroutine TMarch_Substep_Anelastic_Explicit()
    !     use TLab_Pointers_2D, only: pxy_q, pxy_s
    !     use Thermo_Anelastic, only: ribackground

    !     real(wp) alpha
    !     integer iq, is, k

    !     ! #######################################################################
    !     if (rkm_substep < rkm_endstep) then
    !         alpha = kco(rkm_substep)

    !         do iq = 1, inb_flow
    !             do k = 1, kmax
    !                 pxy_q(:, k, iq) = pxy_q(:, k, iq) + dte*pxy_hq(:, k, iq)*ribackground(k)
    !                 pxy_hq(:, k, iq) = pxy_hq(:, k, iq)*alpha
    !             end do
    !         end do

    !         do iq = 1, inb_scal
    !             do k = 1, kmax
    !                 pxy_s(:, k, iq) = pxy_s(:, k, iq) + dte*pxy_hs(:, k, iq)*ribackground(k)
    !                 pxy_hs(:, k, iq) = pxy_hs(:, k, iq)*alpha
    !             end do
    !         end do

    !     else
    !         do is = 1, inb_flow
    !             do k = 1, kmax
    !                 pxy_q(:, k, is) = pxy_q(:, k, is) + dte*pxy_hq(:, k, is)*ribackground(k)
    !             end do
    !         end do

    !         do is = 1, inb_scal
    !             do k = 1, kmax
    !                 pxy_s(:, k, is) = pxy_s(:, k, is) + dte*pxy_hs(:, k, is)*ribackground(k)
    !             end do
    !         end do
    !     end if

    !     return
    ! end subroutine TMarch_Substep_Anelastic_Explicit

    !########################################################################
    !#
    !# Determine the variable time step.
    !# For constant time step, this routine
    !# calculates CFL and diffustion numbers for log files
    !#
    !# The diffusion number is fixed in terms of the CFL, which is the input.
    !# This depends on the scheme used. From Lele (1992), page 32, we have that, if
    !# the sixth order tridiagonal scheme is used, then the maximum CFL number
    !# for a 4RK is 2.9/1.989, about 1.43. For the (5)4RK from CarpenterKennedy1994
    !# used here we have 3.36/1.989, about 1.69.
    !# This holds for periodic case. A safety margin leads to the common value of 1.2.
    !#
    !# If second order finite different operator is used, then the maximum
    !# diffusion number is 2.9/6.857, about 0.42.
    !# For the (5)4RK from CarpenterKennedy1994 it is 4.639/6.857 = 0.68
    !# If the extension by Lamballais et al is used, then the maximum
    !# diffusion number is 2.9/pi^2, about 0.29.
    !# For the (5)4RK from CarpenterKennedy1994 it is 4.639/pi^2 = 0.47.
    !#
    !# If twice the first order finite difference operator is used, then the
    !# maximum diffusion number is 2.9/1.989^2, about 0.73.
    !# For the (5)4RK from CarpenterKennedy1994 it is 4.639/1.989^2 = 1.17
    !#
    !########################################################################
    subroutine TMarch_Courant()
        use TLab_Grid, only: y, xSubgrid, ySubgrid
        use DNS_Control, only: logs_data, logs_dtime
        use TLab_Pointers_3D, only: u, v, w, p_wrk3d

        ! -------------------------------------------------------------------
        integer(wi) i, j, k
        integer(wi) ipmax, j_glo
        real(wp) pmax(3), dtc, dtd
#ifdef USE_MPI
        real(wp) pmax_aux(3)
#endif

        ! ###################################################################
        dtc = big_wp    ! So that the minimum non-zero determines dt at the end
        dtd = big_wp

        ipmax = 0       ! Initialize counter of time constraints

        ! ###################################################################
        ! CFL number condition
        ! ###################################################################
        ipmax = ipmax + 1

        ! -------------------------------------------------------------------
        ! Incompressible: Calculate global maximum of u/dx + v/dy + w/dz
        ! -------------------------------------------------------------------
        select case (nse_eqns)
        case (DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC)
            if (y%size > 1) then
                do k = 1, kmax
                    do j = 1, jmax
                        j_glo = j + ySubgrid%offset
                        do i = 1, imax
                            p_wrk3d(i, j, k) = abs(u(i, j, k))*ds(1)%one_ov_ds1(i + xSubgrid%offset) &
                                               + abs(v(i, j, k))*ds(2)%one_ov_ds1(j_glo) &
                                               + abs(w(i, j, k))*ds(3)%one_ov_ds1(k)
                        end do
                    end do
                end do
            else    ! do I need this?
                do k = 1, kmax
                    do j = 1, jmax
                        do i = 1, imax
                            p_wrk3d(i, j, k) = abs(u(i, j, k))*ds(1)%one_ov_ds1(i + xSubgrid%offset) &
                                               + abs(w(i, j, k))*ds(3)%one_ov_ds1(k)
                        end do
                    end do
                end do
            end if

        end select

        pmax(1) = maxval(p_wrk3d)

        ! ###################################################################
        ! Diffusion number condition
        ! ###################################################################
        ipmax = ipmax + 1

        ! -------------------------------------------------------------------
        ! Incompressible: Calculate global maximum of \mu*(1/dx^2 + 1/dy^2 + 1/dz^2)
        ! -------------------------------------------------------------------
        select case (nse_eqns)
        case (DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC)
            pmax(2) = schmidtfactor*visc*dx2i

        end select

        ! ###################################################################
        ! Final operations
        ! ###################################################################
#ifdef USE_MPI
        call MPI_ALLREDUCE(pmax, pmax_aux, ipmax, MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ims_err)
        pmax(1:ipmax) = pmax_aux(1:ipmax)
#endif

        if (use_variable_timestep) then
            if (pmax(1) > 0.0_wp) dtc = cfla/pmax(1) ! Set time step for the given CFL number
            if (pmax(2) > 0.0_wp) dtd = cfld/pmax(2) ! Set time step for the given diffusion number

            dtime = min(dtc, big_wp)
            ! select case (rkm_mode)
            ! case (RKM_EXP3, RKM_EXP4)       ! Explicit diffusion
            dtime = min(dtd, dtime)

            ! end select

        end if

        ! Real CFL and diffusion numbers being used, for the logfile
        logs_dtime = dtime
        logs_data(2) = dtime*pmax(1)
        logs_data(3) = dtime*pmax(2)

        return

    end subroutine TMarch_Courant

end module TimeMarching
