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
    ! use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
    use TLabMPI_Transpose_X, only: TLabMPI_Trp_Initialize_X
#endif
    use IO_Fields
    use TLab_Grid
    use FDM, only: FDM_Initialize
    use NavierStokes
    use Thermodynamics!, only: Thermo_Initialize
    use TLab_Background, only: TLab_Initialize_Background
    use Gravity, only: Gravity_Initialize
    use SpecialForcing, only: SpecialForcing_Initialize
    use Rotation, only: Rotation_Initialize
    use Microphysics, only: Microphysics_Initialize
    use Radiation, only: Radiation_Initialize
    use LargeScaleForcing, only: LargeScaleForcing_Initialize
    use OPR_Partial
    use OPR_Fourier, only: OPR_Fourier_Initialize
    use OPR_Elliptic, only: OPR_Elliptic_Initialize
    use NSE_Burgers, only: NSE_Burgers_Initialize
    use NSE_Pressure
    use FI_VECTORCALCULUS
    use FI_STRAIN_EQN
    use FI_VORTICITY_EQN
    use Tensor
    use Diagnostics
    use Reductions, only: Reduce_Block_InPlace
    use StatsPDFs

    implicit none

    integer(wi), parameter :: itime_size_max = 3000 ! iterations to be processed
    integer(wi) itime_size, it
    integer(wi) itime_vec(itime_size_max)

    integer(wi), parameter :: iopt_size_max = 30    ! options to be processed
    integer(wi) iopt_size, iv
    integer(wi) opt_vec(iopt_size_max)
    character(len=64) opt_name(iopt_size_max)

    ! -------------------------------------------------------------------
    ! Additional local arrays
    real(wp), allocatable :: pdf(:), z_aux(:)

    integer(1), allocatable :: gate(:)

    character*32 fname
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
    ! integer(1) gate_level

    type(pointers_dt), allocatable :: vars(:)
    ! integer :: ifield

    !########################################################################
    !########################################################################
    call TLab_Start()

    call TLab_Initialize_Parameters(ifile)
    call IO_Initialize()

    call TLab_Grid_Initialize()

#ifdef USE_MPI
    call TLabMPI_Initialize(ifile)
    ! call TLabMPI_Trp_Initialize(ifile)
    call TLabMPI_Trp_Initialize_X(ifile)
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

    ! ###################################################################
    ! Postprocess given list of files
    ! ###################################################################
    ibc(:) = 2                  ! default is to consider local interval, with analysis; see pdf1v_n
    allocate (vars(nfield))

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

        iv = 1      ! We only process 1 type per run
        select case (trim(adjustl(opt_name(opt_vec(iv)))))
        case ('Main variables')
            write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))

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

        case ('Scalar gradient G_iG_i/2 equation')
            write (fname, *) itime; fname = 'pdfG2'//trim(adjustl(fname))
            call Diagnose_ScalarGradientEquation(is=inb_scal, vars=vars)

        case ('Enstrophy W_iW_i/2 equation')
            write (fname, *) itime; fname = 'pdfW2'//trim(adjustl(fname))
            call Diagnose_EnstrophyEquation(vars=vars)

        case ('Strain 2S_ijS_ij/2 equation')
            write (fname, *) itime; fname = 'pdfS2'//trim(adjustl(fname))
            call Diagnose_StrainEquation(vars=vars)

        case ('Velocity gradient invariants')
            write (fname, *) itime; fname = 'pdfInv'//trim(adjustl(fname))

            call FI_INVARIANT_R(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            call FI_INVARIANT_Q(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
            call FI_INVARIANT_P(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 3), txc(1, 4))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'InvP'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'InvQ'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'InvR'

            if (kmax_aux*opt_block /= z%size .and. reduce_data) then ! I already need it here
                do is = 1, ifield
                    call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, 1, 1, 1, imax, jmax*opt_block, kmax_aux, vars(is)%field)
                end do
                reduce_data = .false.
            end if

            write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))//'.RQ'
            call PDF2V(fname, imax, jmax*opt_block, kmax_aux, opt_bins, z_aux, txc(1, 1), txc(1, 2), pdf)

            ! ###################################################################
            ! Joint PDF W^2 and 2S^2
            ! ###################################################################
        case ('Joint enstrophy and strain')
            write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))
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
            call PDF2V(fname, imax, jmax*opt_block, kmax_aux, opt_bins, z_aux, txc(1, 1), txc(1, 2), pdf)

        case ('Joint scalar and scalar gradient')
            write (fname, *) itime; fname = 'pdfSGi'//trim(adjustl(fname))

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
            call PDF2V(fname, imax, jmax*opt_block, kmax_aux, opt_bins, s(1, 1), txc(1, 2), z_aux, pdf)

            ! write (fname, *) itime; fname = 'cavgGiGi'//trim(adjustl(fname))
            ! call CAVG1V_N(fname, rtime, imax*opt_block, kmax_aux, kmax, &
            !               1, opt_bins(1), ibc, vmin, vmax, vars, gate_level, gate, txc(1, 1), z_aux, pdf)

            ! write (fname, *) itime; fname = 'cavgLnGiGi'//trim(adjustl(fname))
            ! call CAVG1V_N(fname, rtime, imax*opt_block, kmax_aux, kmax, &
            !               1, opt_bins(1), ibc, vmin, vmax, vars, gate_level, gate, txc(1, 2), z_aux, pdf)

        case ('Scalar gradient components')
            write (fname, *) itime; fname = 'pdfGi'//trim(adjustl(fname))

            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, s, txc(1, 1))
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, s, txc(1, 2))
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, s, txc(1, 3))
            ! Angles; s array is overwritten to save space
            do ij = 1, isize_field
                dummy = txc(ij, 3)/sqrt(txc(ij, 1)*txc(ij, 1) + txc(ij, 2)*txc(ij, 2) + txc(ij, 3)*txc(ij, 3))
                txc(ij, 4) = asin(dummy)                    ! with Oz
                s(ij, 1) = atan2(txc(ij, 2), txc(ij, 1))    ! with Ox in plane xOy
            end do

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'Gx'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'Gy'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'Gz'
            ifield = ifield + 1; vars(ifield)%field => s(:, 1); vars(ifield)%tag = 'Gtheta'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 4); vars(ifield)%tag = 'Gphi'
            ibc(1:ifield) = 2

            write (fname, *) itime; fname = 'pdf'//trim(adjustl(fname))//'.GphiS'
            call PDF2V(fname, imax, jmax*opt_block, kmax_aux, opt_bins, s(1, 1), txc(1, 4), z_aux, pdf)

        case ('Eigenvalues of rate-of-strain tensor')
            write (fname, *) itime; fname = 'pdfEig'//trim(adjustl(fname))

            call FI_STRAIN_TENSOR(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
            call TENSOR_EIGENVALUES(imax, jmax, kmax, txc(1, 1), txc(1, 7)) ! txc7-txc9

            ifield = ifield + 1; vars(ifield)%field => txc(:, 7); vars(ifield)%tag = 'Lambda1'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'Lambda2'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 9); vars(ifield)%tag = 'Lambda3'

        case ('Eigenframe of rate-of-strain tensor')
            write (fname, *) itime; fname = 'pdfCos'//trim(adjustl(fname))

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

            ifield = ifield + 1; vars(ifield)%field => q(:, 1); vars(ifield)%tag = 'cos(w,lambda1)'
            ifield = ifield + 1; vars(ifield)%field => q(:, 2); vars(ifield)%tag = 'cos(w,lambda2)'
            ifield = ifield + 1; vars(ifield)%field => q(:, 3); vars(ifield)%tag = 'cos(w,lambda3)'

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

            ifield = ifield + 1; vars(ifield)%field => txc(:, 7); vars(ifield)%tag = 'cos(G,lambda1)'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 8); vars(ifield)%tag = 'cos(G,lambda2)'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 9); vars(ifield)%tag = 'cos(G,lambda3)'

        case ('Longitudinal velocity derivatives')
            write (fname, *) itime; fname = 'pdfUDer'//trim(adjustl(fname))

            call OPR_Partial_X(OPR_P1, imax, jmax, kmax, q(1, 1), txc(1, 1))
            call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, q(1, 2), txc(1, 2))
            call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, q(1, 3), txc(1, 3))

            ifield = ifield + 1; vars(ifield)%field => txc(:, 1); vars(ifield)%tag = 'Sxx'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 2); vars(ifield)%tag = 'Syy'
            ifield = ifield + 1; vars(ifield)%field => txc(:, 3); vars(ifield)%tag = 'Szz'

        case ('Potential vorticity')
            write (fname, *) itime; fname = 'pdfPV'//trim(adjustl(fname))
            call Diagnose_PotentialEnstrophy(vars=vars)

        case ('Thermodynamics')
            write (fname, *) itime; fname = 'pdfThermo'//trim(adjustl(fname))
            call Diagnose_Thermodynamics(vars=vars)

        case ('Atmospheric Thermodynamics')
            select case (imode_thermo)
            case (THERMO_TYPE_ANELASTIC)
                write (fname, *) itime; fname = 'pdfThermoEnergies'//trim(adjustl(fname))
                call Diagnose_Energies_Anelastic(vars)
                if (kmax_aux*opt_block /= z%size .and. reduce_data) then
                    do is = 1, size(vars)
                        call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, 1, 1, 1, imax, jmax*opt_block, kmax_aux, vars(is)%field)
                    end do
                end if
                call PDF1V_N(fname, imax, jmax*opt_block, kmax_aux, &
                             size(vars), opt_bins(1), ibc, vmin, vmax, vars, z_aux, pdf)

                write (fname, *) itime; fname = 'pdfThermoThetas'//trim(adjustl(fname))
                call Diagnose_Thetas_Anelastic(vars)
                if (kmax_aux*opt_block /= z%size .and. reduce_data) then
                    do is = 1, size(vars)
                        call REDUCE_BLOCK_INPLACE(imax, jmax, kmax, 1, 1, 1, imax, jmax*opt_block, kmax_aux, vars(is)%field)
                    end do
                end if
                call PDF1V_N(fname, imax, jmax*opt_block, kmax_aux, &
                             size(vars), opt_bins(1), ibc, vmin, vmax, vars, z_aux, pdf)

                write (fname, *) itime; fname = 'pdfThermoMoist'//trim(adjustl(fname))
                call Diagnose_Moisture_Anelastic(vars)

            case (THERMO_TYPE_COMPRESSIBLE)
            end select

        end select

        ! ###################################################################
        if (ifield == 0) ifield = size(vars)

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

            call PDF1V_N(fname, imax, jmax*opt_block, kmax_aux, &
                         ifield, opt_bins(1), ibc, vmin, vmax, vars, z_aux, pdf)

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
        iv = 0
        iv = iv + 1; opt_name(iv) = 'Main variables'
        iv = iv + 1; opt_name(iv) = 'Scalar gradient G_iG_i/2 equation'
        iv = iv + 1; opt_name(iv) = 'Enstrophy W_iW_i/2 equation'
        iv = iv + 1; opt_name(iv) = 'Strain 2S_ijS_ij/2 equation'
        iv = iv + 1; opt_name(iv) = 'Velocity gradient invariants'
        iv = iv + 1; opt_name(iv) = 'Joint enstrophy and strain'
        iv = iv + 1; opt_name(iv) = 'Joint scalar and scalar gradient'
        iv = iv + 1; opt_name(iv) = 'Scalar gradient components'
        iv = iv + 1; opt_name(iv) = 'Eigenvalues of rate-of-strain tensor'
        iv = iv + 1; opt_name(iv) = 'Eigenframe of rate-of-strain tensor'
        iv = iv + 1; opt_name(iv) = 'Longitudinal velocity derivatives'
        iv = iv + 1; opt_name(iv) = 'Potential vorticity'
        iv = iv + 1; opt_name(iv) = 'Thermodynamics'
        iv = iv + 1; opt_name(iv) = 'Atmospheric Thermodynamics'
        if (iv > iopt_size_max) then ! Check
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Increase number of options.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
        end if

        call ScanFile_Char(bakfile, ifile, 'PostProcessing', 'ParamPdfs', '-1', sRes)
        if (sRes == '-1') then
#ifdef USE_MPI
#else
            write (*, '(A)') 'Option?'
            do is = 1, iv
                write (*, '(I2,A)') is, '. '//trim(adjustl(opt_name(is)))
            end do
            read (*, '(A512)') sRes
#endif
        end if
        opt_vec(:) = -1
        iopt_size = iopt_size_max
        call LIST_INTEGER(sRes, iopt_size, opt_vec)

        if (opt_vec(1) < 0) then ! Check
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Missing ParamPdfs input.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
        end if

        ! -------------------------------------------------------------------
        opt_main = opt_vec(1) ! default values
        opt_block = 1
        opt_bins = 16
        ! gate_level = 0

        call ScanFile_Char(bakfile, ifile, 'PostProcessing', 'ParamPdfs', '-1', sRes)
        if (sRes == '-1') then
#ifdef USE_MPI
#else
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
            if (iopt_size >= 2) opt_block = opt_vec(2)
            if (iopt_size >= 3) opt_bins = opt_vec(3:4)
            ! if (iopt_size >= 3) gate_level = int(opt_vec(3), KIND=1)

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
        isize_pdf = opt_bins(1) + 2

        if (any([DNS_EQNS_BOUSSINESQ, DNS_EQNS_ANELASTIC] == nse_eqns)) then
            inb_txc = max(inb_txc, 6)
        else
            inb_txc = max(inb_txc, 1)
        end if

        iv = 1
        select case (trim(adjustl(opt_name(opt_vec(iv)))))
        case ('Main variables')
            iread_scal = .true.; iread_flow = .true.; inb_txc = max(inb_txc, 6)
            nfield = 4 + inb_scal_array; isize_pdf = opt_bins(1) + 2
            if (nse_eqns == DNS_EQNS_COMPRESSIBLE) nfield = nfield + 2
        case ('Scalar gradient G_iG_i/2 equation')
            iread_scal = .true.; iread_flow = .true.; inb_txc = max(inb_txc, 7)
            nfield = 5; isize_pdf = opt_bins(1) + 2
        case ('Enstrophy W_iW_i/2 equation')
            iread_scal = .true.; iread_flow = .true.; inb_txc = max(inb_txc, 8)
            nfield = 7; isize_pdf = opt_bins(1) + 2
        case ('Strain 2S_ijS_ij/2 equation')
            iread_scal = .true.; iread_flow = .true.; inb_txc = max(inb_txc, 8)
            nfield = 5; isize_pdf = opt_bins(1) + 2
        case ('Velocity gradient invariants')
            iread_flow = .true.; inb_txc = max(inb_txc, 6)
            nfield = 3; isize_pdf = opt_bins(1)*opt_bins(2) + 2 + 2*opt_bins(1)
        case ('Joint enstrophy and strain')
            iread_flow = .true.; inb_txc = max(inb_txc, 4)
            nfield = 2; isize_pdf = opt_bins(1)*opt_bins(2) + 2 + 2*opt_bins(1)
        case ('Joint scalar and scalar gradient')
            iread_scal = .true.; iread_flow = .true.; inb_txc = max(inb_txc, 3)
            nfield = 2; isize_pdf = opt_bins(1)*opt_bins(2) + 2 + 2*opt_bins(1)
        case ('Scalar gradient components')
            iread_scal = .true.; inb_txc = max(inb_txc, 4)
            nfield = 5; isize_pdf = opt_bins(1)*opt_bins(2) + 2 + 2*opt_bins(1)
        case ('Eigenvalues of rate-of-strain tensor')
            iread_flow = .true.; inb_txc = max(inb_txc, 9)
            nfield = 3; isize_pdf = opt_bins(1) + 2
        case ('Eigenframe of rate-of-strain tensor')
            iread_scal = .true.; iread_flow = .true.; inb_txc = max(inb_txc, 9)
            nfield = 6; isize_pdf = opt_bins(1) + 2
        case ('Longitudinal velocity derivatives')
            iread_flow = .true.; inb_txc = max(inb_txc, 3)
            nfield = 3; isize_pdf = opt_bins(1) + 2
        case ('Potential vorticity')
            iread_scal = .true.; iread_flow = .true.; inb_txc = max(inb_txc, 6)
            nfield = 2
        case ('Thermodynamics')
            iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 4)
            nfield = 2
        case ('Atmospheric Thermodynamics')
            iread_scal = .true.; inb_txc = max(inb_txc, 3)
            nfield = 3
        end select

        return
    end subroutine
end program PDFS
