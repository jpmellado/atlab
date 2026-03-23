#include "tlab_error.h"

program PDFS
    use TLab_Constants, only: wp, wi, small_wp
    use TLab_Constants, only: ifile, efile, lfile, gfile, tag_flow, tag_scal
    use TLab_Pointers, only: pointers_dt
    use TLab_Time, only: itime, rtime
    use TLab_Memory, only: imax, jmax, kmax, inb_scal_array, inb_txc, inb_flow, inb_scal, isize_field, inb_wrk2d
    use TLab_Arrays
    use TLab_WorkFlow, only: TLab_Write_ASCII, TLab_Stop, TLab_Start
    use TLab_Memory, only: TLab_Initialize_Memory
#ifdef USE_MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
#endif
    use IO_Fields
    use TLab_Grid
    use FDM, only: FDM_Initialize
    use NavierStokes
    use Thermodynamics, only: Thermo_Initialize
    use TLab_Background, only: TLab_Initialize_Background
    use Gravity, only: Gravity_Initialize, gravityProps, Gravity_AddSource, bbackground
    use SpecialForcing, only: SpecialForcing_Initialize
    use Rotation, only: Rotation_Initialize
    use Microphysics, only: Microphysics_Initialize
    use Radiation, only: Radiation_Initialize
    use LargeScaleForcing, only: LargeScaleForcing_Initialize
    use OPR_Partial
    use OPR_Fourier
    use OPR_Elliptic
    use NSE_Burgers, only: NSE_Burgers_Initialize
    use NSE_Pressure
    use FI_VECTORCALCULUS
    use FI_STRAIN_EQN
    use FI_GRADIENT_EQN
    use FI_VORTICITY_EQN
    use Tensor
    use Reductions, only: Reduce_Block_InPlace

    implicit none

    integer(wi), parameter :: itime_size_max = 3000 ! iterations to be processed
    integer(wi) itime_size, it
    integer(wi) itime_vec(itime_size_max)

    integer(wi), parameter :: iopt_size_max = 30    ! options to be processed
    integer(wi) iopt_size
    integer(wi) opt_vec(iopt_size_max)
    ! real(wp) opt_vec2(iopt_size_max)

    integer(wi), parameter :: params_size_max = 2

    ! -------------------------------------------------------------------
    ! Additional local arrays
    real(wp), allocatable :: pdf(:), z_aux(:)

    integer(1), allocatable :: gate(:)
    type(pointers_dt) :: vars(16)

    character*32 fname!, bakfile
    character*64 str

    integer opt_main, opt_block, opt_bins(2)
    integer nfield, ifield, is, ij, kmax_aux
    integer isize_pdf
    real(wp) dummy, eloc1, eloc2, eloc3, cos1, cos2, cos3
    logical iread_flow, iread_scal
    real(wp) params(2)
    integer(wi) ibc(16)
    real(wp) vmin(16), vmax(16)
    logical reduce_data

    ! ! Gates for the definition of the intermittency function (partition of the fields)
    ! integer(wi) opt_cond, opt_cond_scal, opt_cond_relative
    ! integer(wi), parameter :: igate_size_max = 8
    ! integer(wi) igate_size
    ! real(wp) gate_threshold(igate_size_max)
    integer(1) gate_level

    !########################################################################
    !########################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
    call IO_Initialize()

    call TLab_Grid_Initialize()

#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    call TLabMPI_Trp_Initialize(ifile)
#endif

    call NavierStokes_Initialize_Parameters(ifile)
    call Thermo_Initialize(ifile)

    call Gravity_Initialize(ifile)
    call Rotation_Initialize(ifile)
    call SpecialForcing_Initialize(ifile)
    call Microphysics_Initialize(ifile)
    call Radiation_Initialize(ifile)
    call LargeScaleForcing_Initialize(ifile)

    call TLab_Consistency_Check()

    call Pdfs_Initialize()

    ! #######################################################################
    call TLab_Initialize_Memory(__FILE__)

    call FDM_Initialize(ifile)

    call OPR_Partial_Initialize(ifile)
    call OPR_Fourier_Initialize()
    call OPR_Elliptic_Initialize(ifile)
    call OPR_Check()

    call TLab_Initialize_Background(ifile)  ! Initialize thermodynamic quantities
    call NSE_Burgers_Initialize(ifile)

    ! do ig = 1, 3
    !     call OPR_FILTER_INITIALIZE(g(ig), PressureFilter(ig))
    ! end do

    allocate (gate(isize_field))

    ! in case z%size is not divisible by opt_block, drop the upper most planes
    kmax_aux = z%size/opt_block
    allocate (z_aux(kmax_aux))          ! Reduced vertical grid
    z_aux(:) = 0.0_wp
    do ij = 1, kmax_aux*opt_block
        is = (ij - 1)/opt_block + 1
        z_aux(is) = z_aux(is) + z%nodes(ij)/real(opt_block, wp)
    end do

    ! Space for the min and max of sampling variable at opt_bins+1,opt_bins+2
    ! Space for the 3D pdf at kmax_aux+1
    allocate (pdf(isize_pdf*(kmax_aux + 1)*nfield))

    ibc(1:nfield) = 1

    ! ###################################################################
    ! Postprocess given list of files
    ! ###################################################################
    do it = 1, itime_size
        itime = itime_vec(it)

        write (str, *) itime; str = 'Processing iteration It'//trim(adjustl(str))
        call TLab_Write_ASCII(lfile, str)

        if (iread_scal) then
            write (fname, *) itime; fname = trim(adjustl(tag_scal))//trim(adjustl(fname))
            call IO_Read_Fields(fname, imax, jmax, kmax, itime, inb_scal, 0, s, params(1:1))
            rtime = params(1)
        end if

        if (iread_flow) then
            write (fname, *) itime; fname = trim(adjustl(tag_flow))//trim(adjustl(fname))
            call IO_Read_Fields(fname, imax, jmax, kmax, itime, inb_flow, 0, q, params(1:1))
            rtime = params(1)
        end if

        call TLab_Diagnostic(imax, jmax, kmax, s)  ! Initialize diagnostic thermodynamic quantities

        ! ! -------------------------------------------------------------------
        ! ! Calculate intermittency
        ! ! -------------------------------------------------------------------
        ! if (opt_cond == 1) then ! External file
        !     write (fname, *) itime; fname = 'gate.'//trim(adjustl(fname))
        !     call IO_Read_Field_INT1(fname, imax, jmax, kmax, itime, gate, params(1:2))
        !     igate_size = int(params(2))

        ! else if (opt_cond > 1) then
        !     opt_cond_scal = 1 ! Scalar field to use for the conditioning
        !     if (imixture == MIXT_TYPE_AIRWATER .or. imixture == MIXT_TYPE_AIRWATER_LINEAR) then
        !         opt_cond_scal = inb_scal_array
        !     end if

        !     call TLab_Write_ASCII(lfile, 'Calculating partition...')
        !     call FI_GATE(opt_cond, opt_cond_relative, opt_cond_scal, &
        !                  imax, jmax, kmax, igate_size, gate_threshold, q, s, txc, gate)

        !     if (kmax_aux*opt_block /= z%size) then
                !   call REDUCE_BLOCK_INPLACE_INT1(imax, jmax, kmax, 1, 1, 1, imax, jmax*opt_block, kmax_aux, gate, wrk1d)
        !     end if

        ! end if

        ! -------------------------------------------------------------------
        ! Type of PDFs
        ! -------------------------------------------------------------------
        ifield = 0
        reduce_data = .true.

        select case (opt_main)

            ! ###################################################################
            ! Main variable 2D-PDF
            ! ###################################################################
        case (1)
            ifield = ifield + 1; vars(ifield)%field => q(:, 1); vars(ifield)%tag = 'u'
            ifield = ifield + 1; vars(ifield)%field => q(:, 2); vars(ifield)%tag = 'v'
            ifield = ifield + 1; vars(ifield)%field => q(:, 3); vars(ifield)%tag = 'w'
            if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
                call NSE_Pressure_Incompressible(q, s, txc(:, 1), txc(:, 2), txc(:, 5), txc(:, 6))
                ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'p'
            else
                ifield = ifield + 1; vars(ifield)%field => q(:, 6); vars(ifield)%tag = 'p'
                ifield = ifield + 1; vars(ifield)%field => q(:, 5); vars(ifield)%tag = 'r'
                ifield = ifield + 1; vars(ifield)%field => q(:, 7); vars(ifield)%tag = 't'
            end if

            do is = 1, inb_scal_array
                ifield = ifield + 1; vars(ifield)%field => s(:, is); vars(ifield)%tag = 's'
                write (str, *) is; vars(ifield)%tag = trim(adjustl(vars(ifield)%tag))//trim(adjustl(str))
            end do

            do is = 1, ifield ! In case we want same interval for all heights
                if (ibc(is) == 0) call MINMAX(imax, jmax, kmax, vars(is)%field, vmin(is), vmax(is))
            end do

            ! ###################################################################
            ! Scalar gradient equation
            ! ###################################################################
        case (2)
            call TLab_Write_ASCII(lfile, 'Computing scalar gradient equation...')

            call FI_GRADIENT_PRODUCTION(imax, jmax, kmax, s, q(1, 1), q(1, 2), q(1, 3), &
                                        txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            call FI_GRADIENT_DIFFUSION(imax, jmax, kmax, s, & ! array q used as auxiliar
                                       txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), q(1, 1))
            txc(1:isize_field, 2) = txc(1:isize_field, 2)*visc/schmidt(inb_scal)
            call FI_GRADIENT(imax, jmax, kmax, s, txc(1, 3), txc(1, 4))
            txc(1:isize_field, 5) = txc(1:isize_field, 1)/txc(1:isize_field, 3)
            txc(1:isize_field, 4) = log(txc(1:isize_field, 3))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'GiGi'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'LnGiGi'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'ProductionMsGiGjSij'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'DiffusionNuGiLapGi'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 5); vars(ifield)%tag = 'StrainAMsNiNjSij'; ibc(ifield) = 2

            ! ###################################################################
            ! Enstrophy equation
            ! ###################################################################
        case (3)
            call TLab_Write_ASCII(lfile, 'Computing enstrophy equation...')

            if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
                ! txc(:, 4) = 0.0_wp; txc(:, 5) = 0.0_wp; txc(:, 6) = 0.0_wp

                select case (nse_eqns)
                case (DNS_EQNS_ANELASTIC)
                    ! call Thermo_Anelastic_BUOYANCY(imax, jmax, kmax, s, wrk3d)

                case (DNS_EQNS_BOUSSINESQ)
                    wrk1d(1:kmax, 1) = bbackground(1:kmax)
                    bbackground(1:kmax) = 0.0_wp
                    wrk3d(1:isize_field) = 0.0_wp
                    call Gravity_AddSource(gravityProps, imax, jmax, kmax, s, wrk3d, gravityProps%vector(3))
                    bbackground(1:kmax) = wrk1d(1:kmax, 1)
                end select
                s(1:isize_field, 1) = wrk3d(1:isize_field)

                call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, s, txc(1, 4))
                txc(:, 4) = -txc(:, 4)
                txc(:, 5) = 0.0_wp
                call OPR_Partial_X(OPR_P1, imax, jmax, kmax, s, txc(1, 6))

            else
                call FI_VORTICITY_BAROCLINIC(imax, jmax, kmax, q(1, 5), q(1, 6), txc(1, 4), txc(1, 3), txc(1, 7))
            end if
            ! result vector in txc1, txc2, txc3
            call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 7))
            ! scalar product, store in txc8
            txc(1:isize_field, 8) = txc(1:isize_field, 1)*txc(1:isize_field, 4) &
                                    + txc(1:isize_field, 2)*txc(1:isize_field, 5) + txc(1:isize_field, 3)*txc(1:isize_field, 6)

            call FI_VORTICITY_PRODUCTION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                                         txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))

            call FI_VORTICITY_DIFFUSION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 2), &
                                        txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7))
            txc(1:isize_field, 2) = visc*txc(1:isize_field, 2)

            call FI_VORTICITY(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 3), txc(1, 4), txc(1, 5))

            call FI_INVARIANT_P(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 4), txc(1, 5))

            txc(1:isize_field, 5) = txc(1:isize_field, 4)*txc(1:isize_field, 3) ! -w^2 div(u)
            txc(1:isize_field, 4) = txc(1:isize_field, 1)/txc(1:isize_field, 3) ! production rate
            txc(1:isize_field, 6) = log(txc(1:isize_field, 3))                  ! ln(w^2)

            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'WiWi'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 6); vars(ifield)%tag = 'LnWiWi'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'ProductionWiWjSij'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'DiffusionNuWiLapWi'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 5); vars(ifield)%tag = 'DilatationMsWiWiDivU'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'Baroclinic'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'RateANiNjSij'; ibc(ifield) = 2

            ! ###################################################################
            ! Strain equation
            ! ###################################################################
        case (4)
            call TLab_Write_ASCII(lfile, 'Computing strain equation...')

            if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
                call NSE_Pressure_Incompressible(q, s, txc(:, 1), txc(:, 2), txc(:, 5), txc(:, 6))
                call FI_STRAIN_PRESSURE(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                                        txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            else
                call FI_STRAIN_PRESSURE(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), q(1, 6), &
                                        txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            end if
            txc(1:isize_field, 1) = 2.0_wp*txc(1:isize_field, 2)

            call FI_STRAIN_PRODUCTION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), &
                                      txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7))
            txc(1:isize_field, 2) = 2.0_wp*txc(1:isize_field, 2)

            call FI_STRAIN_DIFFUSION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), &
                                     txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6), txc(1, 7), txc(1, 8))
            txc(1:isize_field, 3) = 2.0_wp*visc*txc(1:isize_field, 3)

            call FI_STRAIN(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            txc(1:isize_field, 4) = 2.0_wp*txc(1:isize_field, 4)
            txc(1:isize_field, 5) = log(txc(1:isize_field, 4))

            ifield = ifield + 1; vars(1)%field => txc(:, 4); vars(ifield)%tag = '2SijSij'; ibc(ifield) = 2
            ifield = ifield + 1; vars(2)%field => txc(:, 5); vars(ifield)%tag = 'Ln2SijSij'; ibc(ifield) = 2
            ifield = ifield + 1; vars(3)%field => txc(:, 2); vars(ifield)%tag = 'ProductionMs2SijSjkS_ki'; ibc(ifield) = 2
            ifield = ifield + 1; vars(4)%field => txc(:, 3); vars(ifield)%tag = 'DiffusionNuSijLapSij'; ibc(ifield) = 2
            ifield = ifield + 1; vars(5)%field => txc(:, 1); vars(ifield)%tag = 'Pressure2SijPij'; ibc(ifield) = 2

            ! ###################################################################
            ! Velocity gradient invariants
            ! ###################################################################
        case (5)
            call TLab_Write_ASCII(lfile, 'Computing velocity gradient invariants...')

            call FI_INVARIANT_R(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            call FI_INVARIANT_Q(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
            call FI_INVARIANT_P(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 3), txc(1, 4))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'InvP'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'InvQ'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'InvR'; ibc(ifield) = 2

            if (kmax_aux*opt_block /= z%size .and. reduce_data) then ! I already need it here
                do is = 1, ifield
                    call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, 1, 1, 1, imax, jmax*opt_block, kmax_aux, vars(is)%field)
                end do
                reduce_data = .false.
            end if

            write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))//'.RQ'
            call PDF2V(fname, rtime, imax, jmax*opt_block, kmax_aux, opt_bins, z_aux, txc(1, 1), txc(1, 2), pdf)

            ! ###################################################################
            ! Joint PDF W^2 and 2S^2
            ! ###################################################################
        case (6)
            call TLab_Write_ASCII(lfile, 'Computing enstrophy-strain pdf...')

            call FI_VORTICITY(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3))
            call FI_STRAIN(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 2), txc(1, 3), txc(1, 4))
            txc(1:isize_field, 2) = 2.0_wp*txc(1:isize_field, 2)

            if (kmax_aux*opt_block /= z%size .and. reduce_data) then ! I already need it here
                do is = 1, ifield
                    call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, 1, 1, 1, imax, jmax*opt_block, kmax_aux, vars(is)%field)
                end do
                reduce_data = .false.
            end if

            write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))//'.WS'
            call PDF2V(fname, rtime, imax, jmax*opt_block, kmax_aux, opt_bins, z_aux, txc(1, 1), txc(1, 2), pdf)

            ! ###################################################################
            ! Joint PDF Scalar and Scalar Gradient
            ! ###################################################################
        case (7)
            call TLab_Write_ASCII(lfile, 'Computing scalar-scalar--gradient pdf...')

            call FI_GRADIENT(imax, jmax, kmax, s, txc(1, 1), txc(1, 2))
            txc(1:isize_field, 2) = log(txc(1:isize_field, 1))

            ifield = ifield + 1; vars(1)%field => s(:, 1); vars(ifield)%tag = 's'; ibc(ifield) = 1
            ifield = ifield + 1; vars(2)%field => txc(:, 1); vars(ifield)%tag = 'GiGi'; ibc(ifield) = 2
            ifield = ifield + 1; vars(3)%field => txc(:, 2); vars(ifield)%tag = 'LnGiGi'; ibc(ifield) = 3

            if (kmax_aux*opt_block /= z%size .and. reduce_data) then ! I already need it here
                do is = 1, ifield
                    call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, 1, 1, 1, imax, jmax*opt_block, kmax_aux, vars(is)%field)
                end do
                reduce_data = .false.
            end if

            write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))//'.SLnG'
            call PDF2V(fname, rtime, imax, jmax*opt_block, kmax_aux, opt_bins, s(1, 1), txc(1, 2), z_aux, pdf)

            ! write (fname, *) itime; fname = 'cavgGiGi'//trim(adjustl(fname))
            ! call CAVG1V_N(fname, rtime, imax*opt_block, kmax_aux, kmax, &
            !               1, opt_bins(1), ibc, vmin, vmax, vars, gate_level, gate, txc(1, 1), z_aux, pdf)

            ! write (fname, *) itime; fname = 'cavgLnGiGi'//trim(adjustl(fname))
            ! call CAVG1V_N(fname, rtime, imax*opt_block, kmax_aux, kmax, &
            !               1, opt_bins(1), ibc, vmin, vmax, vars, gate_level, gate, txc(1, 2), z_aux, pdf)

            ! ###################################################################
            ! Scalar gradient components
            ! ###################################################################
            ! Check angles for new frame of reference
        case (8)
            call TLab_Write_ASCII(lfile, 'Computing scalar gradient components...')

            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, s, txc(1, 1))
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, s, txc(1, 2))
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, s, txc(1, 3))
            ! Angles; s array is overwritten to save space
            do ij = 1, isize_field
                dummy = txc(ij, 3)/sqrt(txc(ij, 1)*txc(ij, 1) + txc(ij, 2)*txc(ij, 2) + txc(ij, 3)*txc(ij, 3))
                txc(ij, 4) = asin(dummy)                 ! with Oz
                s(ij, 1) = atan2(txc(ij, 2), txc(ij, 1))  ! with Ox in plane xOy
            end do

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'Gx'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'Gy'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'Gz'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => s(:, 1); vars(ifield)%tag = 'Gtheta'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'Gphi'; ibc(ifield) = 2

            write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))//'.GphiS'
            call PDF2V(fname, rtime, imax, jmax*opt_block, kmax_aux, opt_bins, s(1, 1), txc(1, 4), z_aux, pdf)

            ! ###################################################################
            ! eigenvalues of rate-of-strain tensor
            ! ###################################################################
        case (9)
            call TLab_Write_ASCII(lfile, 'Computing eigenvalues of Sij...')

            call FI_STRAIN_TENSOR(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            call TENSOR_EIGENVALUES(imax, jmax, kmax, txc(1, 1), txc(1, 7)) ! txc7-txc9

            ifield = ifield + 1; vars(ifield)%field => txc(:, 7); vars(ifield)%tag = 'Lambda1'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'Lambda2'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 9); vars(ifield)%tag = 'Lambda3'; ibc(ifield) = 2

            ! ###################################################################
            ! eigenframe of rate-of-strain tensor
            ! ###################################################################
        case (10)
            call TLab_Write_ASCII(lfile, 'Computing eigenframe of Sij...')

            call FI_STRAIN_TENSOR(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            call TENSOR_EIGENVALUES(imax, jmax, kmax, txc(1, 1), txc(1, 7)) ! txc7-txc9
            call TENSOR_EIGENFRAME(imax, jmax, kmax, txc(1, 1), txc(1, 7)) ! txc1-txc6
            call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 7), txc(1, 8), txc(1, 9), txc(1, 10))

            do ij = 1, isize_field ! local direction cosines of vorticity vector
                dummy = sqrt(txc(ij, 7)*txc(ij, 7) + txc(ij, 8)*txc(ij, 8) + txc(ij, 9)*txc(ij, 9))
                q(ij, 1) = (txc(ij, 7)*txc(ij, 1) + txc(ij, 8)*txc(ij, 2) + txc(ij, 9)*txc(ij, 3))/dummy
                q(ij, 2) = (txc(ij, 7)*txc(ij, 4) + txc(ij, 8)*txc(ij, 5) + txc(ij, 9)*txc(ij, 6))/dummy
                eloc1 = txc(ij, 2)*txc(ij, 6) - txc(ij, 5)*txc(ij, 3)
                eloc2 = txc(ij, 3)*txc(ij, 4) - txc(ij, 6)*txc(ij, 1)
                eloc3 = txc(ij, 1)*txc(ij, 5) - txc(ij, 4)*txc(ij, 2)
                q(ij, 3) = (txc(ij, 7)*eloc1 + txc(ij, 8)*eloc2 + txc(ij, 9)*eloc3)/dummy
            end do

            ifield = ifield + 1; vars(ifield)%field => q(:, 1); vars(ifield)%tag = 'cos(w,lambda1)'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => q(:, 2); vars(ifield)%tag = 'cos(w,lambda2)'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => q(:, 3); vars(ifield)%tag = 'cos(w,lambda3)'; ibc(ifield) = 2

            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, s, txc(1, 7))
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, s, txc(1, 8))
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, s, txc(1, 9))

            do ij = 1, isize_field ! local direction cosines of scalar gradient vector
                dummy = sqrt(txc(ij, 7)*txc(ij, 7) + txc(ij, 8)*txc(ij, 8) + txc(ij, 9)*txc(ij, 9))
                cos1 = (txc(ij, 7)*txc(ij, 1) + txc(ij, 8)*txc(ij, 2) + txc(ij, 9)*txc(ij, 3))/dummy
                cos2 = (txc(ij, 7)*txc(ij, 4) + txc(ij, 8)*txc(ij, 5) + txc(ij, 9)*txc(ij, 6))/dummy
                eloc1 = txc(ij, 2)*txc(ij, 6) - txc(ij, 5)*txc(ij, 3)
                eloc2 = txc(ij, 3)*txc(ij, 4) - txc(ij, 6)*txc(ij, 1)
                eloc3 = txc(ij, 1)*txc(ij, 5) - txc(ij, 4)*txc(ij, 2)
                cos3 = (txc(ij, 7)*eloc1 + txc(ij, 8)*eloc2 + txc(ij, 9)*eloc3)/dummy
                txc(ij, 7) = cos1; txc(ij, 8) = cos2; txc(ij, 9) = cos3
            end do

            ifield = ifield + 1; vars(ifield)%field => txc(:, 7); vars(ifield)%tag = 'cos(G,lambda1)'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'cos(G,lambda2)'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 9); vars(ifield)%tag = 'cos(G,lambda3)'; ibc(ifield) = 2

            ! ###################################################################
            ! Longitudinal velocity derivatives
            ! ###################################################################
        case (11)
            call TLab_Write_ASCII(lfile, 'Computing longitudinal velocity derivatives...')

            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, q(1, 1), txc(1, 1))
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, q(1, 2), txc(1, 2))
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, q(1, 3), txc(1, 3))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'Sxx'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'Syy'; ibc(ifield) = 2
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'Szz'; ibc(ifield) = 2

            ! ###################################################################
            ! Potential vorticity
            ! ###################################################################
        case (12)
            call TLab_Write_ASCII(lfile, 'Computing potential vorticity...')

            call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
            txc(1:isize_field, 6) = txc(1:isize_field, 1)*txc(1:isize_field, 1) &
                                    + txc(1:isize_field, 2)*txc(1:isize_field, 2) &
                                    + txc(1:isize_field, 3)*txc(1:isize_field, 3) ! Enstrophy
            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, s(1, 1), txc(1, 4))
            txc(1:isize_field, 1) = txc(1:isize_field, 1)*txc(1:isize_field, 4)
            txc(1:isize_field, 5) = txc(1:isize_field, 4)*txc(1:isize_field, 4) ! norm grad b
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, s(1, 1), txc(1, 4))
            txc(1:isize_field, 1) = txc(1:isize_field, 1) + txc(1:isize_field, 2)*txc(1:isize_field, 4)
            txc(1:isize_field, 5) = txc(1:isize_field, 5) + txc(1:isize_field, 4)*txc(1:isize_field, 4) ! norm grad b
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, s(1, 1), txc(1, 4))
            txc(1:isize_field, 1) = txc(1:isize_field, 1) + txc(1:isize_field, 3)*txc(1:isize_field, 4)
            txc(1:isize_field, 5) = txc(1:isize_field, 5) + txc(1:isize_field, 4)*txc(1:isize_field, 4) ! norm grad b

            txc(1:isize_field, 5) = sqrt(txc(1:isize_field, 5) + small_wp)
            txc(1:isize_field, 6) = sqrt(txc(1:isize_field, 6) + small_wp)
            txc(1:isize_field, 2) = txc(1:isize_field, 1)/(txc(1:isize_field, 5)*txc(1:isize_field, 6)) ! Cosine of angle between 2 vectors

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'PV'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'Cos'

        end select

        ! ###################################################################
        if (ifield > 0) then
            if (nfield < ifield) then
                call TLab_Write_ASCII(efile, __FILE__//'. Array space nfield incorrect.')
                call TLab_Stop(DNS_ERROR_WRKSIZE)
            end if

            if (kmax_aux*opt_block /= z%size .and. reduce_data) then
                do is = 1, ifield
                    call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, 1, 1, 1, imax, jmax*opt_block, kmax_aux, vars(is)%field)
                end do
            end if

            write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))
            call PDF1V_N(fname, rtime, imax, jmax*opt_block, kmax_aux, &
                         ifield, opt_bins(1), ibc, vmin, vmax, vars, gate_level, gate, z_aux, pdf)

        end if

    end do
    call TLab_Stop(0)

contains
    ! #######################################################################
    ! #######################################################################
    subroutine Pdfs_Initialize()

        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        integer idummy

        ! #######################################################################
        ! Read from tlab.ini
        bakfile = trim(adjustl(ifile))//'.bak'
        block = 'PostProcessing'
        eStr = __FILE__//'. '//trim(adjustl(block))//'. '

        ! -------------------------------------------------------------------
        call ScanFile_Char(bakfile, ifile, block, 'Files', '-1', sRes)
        if (sRes == '-1') then
#ifdef USE_MPI
#else
            write (*, *) 'Iteration numbers ?'
            read (*, '(A512)') sRes
#endif
        end if
        itime_size = itime_size_max
        call LIST_INTEGER(sRes, itime_size, itime_vec)

        if (itime_vec(1) < 0) then ! Check
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Missing Files input.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
        end if

        ! -------------------------------------------------------------------
        opt_main = -1 ! default values
        opt_block = 1
        gate_level = 0
        opt_bins = 16

        call ScanFile_Char(bakfile, ifile, 'PostProcessing', 'ParamPdfs', '-1', sRes)
        iopt_size = iopt_size_max
        call LIST_INTEGER(sRes, iopt_size, opt_vec)
        if (sRes == '-1') then
#ifdef USE_MPI
#else
            write (*, *) 'Option ?'
            write (*, *) ' 1. Main variables'
            write (*, *) ' 2. Scalar gradient G_iG_i/2 equation'
            write (*, *) ' 3. Enstrophy W_iW_i/2 equation'
            write (*, *) ' 4. Strain 2S_ijS_ij/2 equation'
            write (*, *) ' 5. Velocity gradient invariants'
            write (*, *) ' 6. Joint enstrophy and strain'
            write (*, *) ' 7. Joint scalar and scalar gradient'
            write (*, *) ' 8. Scalar gradient components'
            write (*, *) ' 9. Eigenvalues of rate-of-strain tensor'
            write (*, *) '10. Eigenframe of rate-of-strain tensor'
            write (*, *) '11. Longitudinal velocity derivatives'
            write (*, *) '12. Potential vorticity'
            read (*, *) opt_main

            write (*, *) 'Planes block size ?'
            read (*, *) opt_block

            ! write (*, *) 'Gate level to be used ?'
            ! read (*, *) gate_level

            write (*, *) 'Number of PDF bins ?'
            read (*, '(A)') sRes
            idummy = 2
            call LIST_INTEGER(sRes, idummy, opt_bins)

#endif
        else
            opt_main = opt_vec(1)
            if (iopt_size >= 2) opt_block = opt_vec(2)
            if (iopt_size >= 3) opt_bins = opt_vec(4:5)
            ! if (iopt_size >= 3) gate_level = int(opt_vec(3), KIND=1)

        end if

        if (opt_main < 0) then ! Check
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Missing input [ParamPdfs] in tlab.ini.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
        end if

        if (opt_block < 1) then
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Invalid value of opt_block.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
        end if

!         ! -------------------------------------------------------------------
!         ! Defining gate levels for conditioning
!         ! -------------------------------------------------------------------
!         opt_cond = 0 ! default values
!         opt_cond_relative = 0
!         igate_size = 0

!         if (gate_level /= 0) then
! #include "dns_read_partition.h"
!             if (opt_cond > 1) inb_txc = max(inb_txc, 5)
!         end if

        ! ###################################################################
        ! Initialization of array sizes
        ! ###################################################################
        iread_flow = .false.
        iread_scal = .false.
        inb_txc = 0
        nfield = 2
        inb_wrk2d = max(inb_wrk2d, 4)

        if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
            inb_txc = max(inb_txc, 6)
        else
            inb_txc = max(inb_txc, 1)
        end if

        select case (opt_main)
        case (1)
            iread_scal = .true.; iread_flow = .true.; inb_txc = max(inb_txc, 6)
            nfield = 4 + inb_scal_array; isize_pdf = opt_bins(1) + 2
            if (nse_eqns == DNS_EQNS_COMPRESSIBLE) nfield = nfield + 2
        case (2) ! Scalar gradient equation
            iread_scal = .true.; iread_flow = .true.; inb_txc = max(inb_txc, 6)
            nfield = 5; isize_pdf = opt_bins(1) + 2
        case (3) ! Enstrophy equation
            iread_scal = .true.; iread_flow = .true.; inb_txc = max(inb_txc, 8)
            nfield = 7; isize_pdf = opt_bins(1) + 2
        case (4) ! Strain equation
            iread_scal = .true.; iread_flow = .true.; inb_txc = max(inb_txc, 8)
            nfield = 5; isize_pdf = opt_bins(1) + 2
        case (5) ! Invariants
            iread_flow = .true.; inb_txc = max(inb_txc, 6)
            nfield = 3; isize_pdf = opt_bins(1)*opt_bins(2) + 2 + 2*opt_bins(1)
        case (6)
            iread_flow = .true.; inb_txc = max(inb_txc, 4)
            nfield = 2; isize_pdf = opt_bins(1)*opt_bins(2) + 2 + 2*opt_bins(1)
        case (7)
            iread_scal = .true.; iread_flow = .true.; inb_txc = max(inb_txc, 3)
            nfield = 2; isize_pdf = opt_bins(1)*opt_bins(2) + 2 + 2*opt_bins(1)
        case (8)
            iread_scal = .true.; inb_txc = max(inb_txc, 4)
            nfield = 5; isize_pdf = opt_bins(1)*opt_bins(2) + 2 + 2*opt_bins(1)
        case (9) ! eigenvalues
            iread_flow = .true.; inb_txc = max(inb_txc, 9)
            nfield = 3; isize_pdf = opt_bins(1) + 2
        case (10) ! eigenframe
            iread_scal = .true.; iread_flow = .true.; inb_txc = max(inb_txc, 9)
            nfield = 6; isize_pdf = opt_bins(1) + 2
        case (11) ! longitudinal velocity derivatives
            iread_flow = .true.; inb_txc = max(inb_txc, 3)
            nfield = 3; isize_pdf = opt_bins(1) + 2
        case (12) ! potential vorticity
            iread_scal = .true.; iread_flow = .true.; inb_txc = max(inb_txc, 6)
            nfield = 2; isize_pdf = opt_bins(1) + 2
        case (13) ! joint s and v
            iread_scal = .true.; iread_flow = .true.; inb_txc = max(inb_txc, 8)
            nfield = 2; isize_pdf = opt_bins(1)*opt_bins(2) + 2 + 2*opt_bins(1)
        end select

        return
    end subroutine
end program PDFS
