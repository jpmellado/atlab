#include "tlab_error.h"

program VISUALS
    use TLab_Constants, only: wp, wi, small_wp, MAX_PARS, fmt_r
    use TLab_Constants, only: ifile, gfile, lfile, efile, wfile, tag_flow, tag_scal
    use TLab_Time, only: itime, rtime
    use TLab_Memory, only: imax, jmax, kmax, inb_scal_array, inb_txc, inb_flow, inb_scal, isize_field
    use TLab_Arrays
    use TLab_WorkFlow
    use TLab_Memory, only: TLab_Initialize_Memory
    use TLab_Pointers, only: pointers_dt
#ifdef USE_MPI
    use TLabMPI_PROCS, only: TLabMPI_Initialize
    ! use TLabMPI_Transpose, only: TLabMPI_Trp_Initialize
    use TLabMPI_Transpose_X, only: TLabMPI_Trp_Initialize_X
#endif
    use IO_Fields
    use TLab_Grid
    use FDM, only: FDM_Initialize
    use NavierStokes!, only: NavierStokes_Initialize_Parameters
    use Thermodynamics!, only: Thermo_Initialize
    use TLab_Background, only: TLab_Initialize_Background
    use Gravity, only: Gravity_Initialize
    use Rotation, only: Rotation_Initialize
    use Thermo_Anelastic, only: ribackground, Thermo_Anelastic_Weight_InPlace
    use Radiation !, only: Radiation_Initialize, infraredProps
    use Microphysics !, only: Microphysics_Initialize, sedimentationProps
    use LargeScaleForcing, only: LargeScaleForcing_Initialize
    use OPR_Partial
    use OPR_Fourier, only: OPR_Fourier_Initialize
    use OPR_Elliptic, only: OPR_Elliptic_Initialize
    use NSE_Burgers, only: NSE_Burgers_Initialize
    use FI_VECTORCALCULUS
    use FI_STRAIN_EQN
    use Tensor
    use Diagnostics

    implicit none

    integer(wi), parameter :: itime_size_max = 3000 ! iterations to be processed
    integer(wi) itime_size, it
    integer(wi) itime_vec(itime_size_max)

    integer(wi), parameter :: iopt_size_max = 30    ! options to be processed
    integer(wi) iopt_size, iv
    integer(wi) opt_vec(iopt_size_max)
    character(len=64) opt_name(iopt_size_max)

    integer :: opt_format                           ! File format
    integer, parameter :: FORMAT_SINGLE = 1         ! Single precision, no headers
    integer, parameter :: FORMAT_GENERAL = 2        ! General IO format

    integer(wi) subdomain(6)                        ! Subdomain to be saved

    character(len=32) time_str                      ! Time stamp
    integer, parameter :: MaskSize = 6

    character*32 flow_file, scal_file, plot_file
    character*64 str

    integer, parameter :: IO_SUBARRAY_VISUALS_XOY = 1
    integer, parameter :: IO_SUBARRAY_VISUALS_YOZ = 2
    integer, parameter :: IO_SUBARRAY_VISUALS_XOZ = 3
    type(io_subarray_dt) :: io_subarrays(3)

    integer(wi) ij, is
    logical iread_flow, iread_scal
    real(wp) params(MAX_PARS)

    ! ! Gates for the definition of the intermittency function (partition of the fields)
    ! integer(wi) opt_cond, opt_cond_scal, opt_cond_relative
    ! integer(wi), parameter :: igate_size_max = 8
    ! integer(wi) igate_size, ig
    ! real(wp) gate_threshold(igate_size_max)
    ! integer(1), dimension(:), allocatable, save :: gate

    type(pointers_dt), allocatable :: vars(:)
    integer :: ifield

    interface Write_Visuals
        procedure Write_Visuals_Rank1, Write_Visuals_Rank2
    end interface

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
    call Radiation_Initialize(ifile)
    call Microphysics_Initialize(ifile)
    call LargeScaleForcing_Initialize(ifile)

    call TLab_Consistency_Check()

    call Visuals_Initialize()

    ! #######################################################################
    call TLab_Initialize_Memory(__FILE__)

    call FDM_Initialize(ifile)

    call OPR_Partial_Initialize(ifile)
    call OPR_Fourier_Initialize()
    call OPR_Elliptic_Initialize(ifile)
    call OPR_Check()

    call TLab_Initialize_Background(ifile)
    call NSE_Burgers_Initialize(ifile)

    ! allocate (gate(isize_field))

    ! ###################################################################
    ! Postprocess given list of files
    ! ###################################################################
    do it = 1, itime_size
        itime = itime_vec(it)

        write (str, *) itime; str = 'Processing iteration It'//trim(adjustl(str))//'.'
        call TLab_Write_ASCII(lfile, str)

        if (scal_on .and. iread_scal) then ! Scalar variables
            write (scal_file, *) itime; scal_file = trim(adjustl(tag_scal))//trim(adjustl(scal_file))
            call IO_Read_Fields(scal_file, imax, jmax, kmax, itime, inb_scal, 0, s, params(1:1))
            rtime = params(1)
        elseif (.not. scal_on) then
            s = 0.0_wp
        end if

        if (iread_flow) then ! Flow variables
            write (flow_file, *) itime; flow_file = trim(adjustl(tag_flow))//trim(adjustl(flow_file))
            call IO_Read_Fields(flow_file, imax, jmax, kmax, itime, inb_flow, 0, q, params(1:1))
            rtime = params(1)
        end if

        call TLab_Diagnostic(imax, jmax, kmax, s)

        write (str, fmt_r) rtime; str = 'Physical time '//trim(adjustl(str))
        call TLab_Write_ASCII(lfile, str)

        ! ! -------------------------------------------------------------------
        ! ! Calculate intermittency
        ! ! -------------------------------------------------------------------
        ! if (opt_cond == 1) then ! Read external file
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
        ! end if

        ! -------------------------------------------------------------------
        ! define time string
        ! -------------------------------------------------------------------
        do ij = MaskSize, 1, -1
            time_str(ij:ij) = '0'
        end do
        write (plot_file, '(I10)') itime
        time_str(MaskSize - len_trim(adjustl(plot_file)) + 1:Masksize) = trim(adjustl(plot_file))

        ! -------------------------------------------------------------------
        ! Loop over options
        ! -------------------------------------------------------------------
        do iv = 1, iopt_size
            plot_file = trim(adjustl(opt_name(opt_vec(iv))))//time_str(1:MaskSize)

            select case (trim(adjustl(opt_name(opt_vec(iv)))))
            case ('VelocityX')
                txc(1:isize_field, 1) = q(1:isize_field, 1)
                call Write_Visuals(plot_file, txc(:, 1:1))

            case ('VelocityY')
                txc(1:isize_field, 1) = q(1:isize_field, 2)
                call Write_Visuals(plot_file, txc(:, 1:1))

            case ('VelocityZ')
                txc(1:isize_field, 1) = q(1:isize_field, 3)
                call Write_Visuals(plot_file, txc(:, 1:1))

            case ('VelocityVector')
                txc(1:isize_field, 1:3) = q(1:isize_field, 1:3)
                call Write_Visuals(plot_file, txc(:, 1:3))

            case ('VelocityMagnitude')
                txc(1:isize_field, 1) = sqrt(Tensor_Dot(q(1:isize_field, 1:3), q(1:isize_field, 1:3)))
                call Write_Visuals(plot_file, txc(:, 1:1))

            case ('Thermodynamics')
                call Diagnose_Thermodynamics(vars=vars)
                do ifield = 1, size(vars)
                    call Write_Visuals(trim(adjustl(vars(ifield)%tag))//time_str(1:MaskSize), &
                                       vars(ifield)%field(:))
                end do

            case ('PressureAnalysis')
                call Diagnose_PressureForce(vars=vars)
                do ifield = 1, size(vars)
                    call Write_Visuals(trim(adjustl(vars(ifield)%tag))//time_str(1:MaskSize), &
                                       vars(ifield)%field(:))
                end do

                call Diagnose_PressurePartition(vars=vars)
                do ifield = 1, size(vars)
                    call Write_Visuals(trim(adjustl(vars(ifield)%tag))//time_str(1:MaskSize), &
                                       vars(ifield)%field(:))
                end do

            case ('Buoyancy')
                call Diagnose_Buoyancy(vars)
                do ifield = 1, size(vars)
                    call Write_Visuals(trim(adjustl(vars(ifield)%tag))//time_str(1:MaskSize), &
                                       vars(ifield)%field(:))
                end do

            case ('Scalars')
                do is = 1, inb_scal_array
                    write (str, *) is; plot_file = 'Scalar'//trim(adjustl(str))//time_str(1:MaskSize)

                    txc(1:isize_field, 1) = s(1:isize_field, is)
                    call Write_Visuals(plot_file, txc(:, 1:1))

                end do

            case ('ScalarGradientVector')
                do is = 1, inb_scal_array
                    write (str, *) is; str = 'Scalar'//trim(adjustl(str))

                    plot_file = trim(adjustl(str))//'GradientVector'//time_str(1:MaskSize)
                    call OPR_Partial_X(OPR_P1, imax, jmax, kmax, s(1, is), txc(1, 1))
                    call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, s(1, is), txc(1, 2))
                    call OPR_Partial_Z(OPR_P1, imax, jmax, kmax, s(1, is), txc(1, 3))
                    call Write_Visuals(plot_file, txc(:, 1:3))
                end do

            case ('ScalarGradient G_iG_i (Log)')
                do is = 1, inb_scal_array
                    write (str, *) is; str = 'Scalar'//trim(adjustl(str))

                    plot_file = 'Log'//trim(adjustl(str))//'Gradient'//time_str(1:MaskSize)
                    call FI_GRADIENT(imax, jmax, kmax, s(1, is), txc(1, 1), txc(1, 2))
                    txc(1:isize_field, 1) = log10(txc(1:isize_field, 1) + small_wp)
                    call Write_Visuals(plot_file, txc(:, 1:1))

                end do

            case ('ScalarGradientEquation')
                do is = 1, inb_scal_array
                    write (str, *) is; str = 'ScalarGradEqn'//trim(adjustl(str))

                    call Diagnose_ScalarGradientEquation(is, vars)
                    do ifield = 1, size(vars)
                        plot_file = trim(adjustl(str))//trim(adjustl(vars(ifield)%tag))//time_str(1:MaskSize)
                        call Write_Visuals(plot_file, vars(ifield)%field(:))
                    end do

                end do

            case ('VorticityVector')
                call FI_CURL(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
                call Write_Visuals(plot_file, txc(:, 1:3))

            case ('Enstrophy W_iW_i (Log)')
                call Diagnose_PotentialEnstrophy(vars=vars)
                do ifield = 1, size(vars)
                    call Write_Visuals(trim(adjustl(vars(ifield)%tag))//time_str(1:MaskSize), &
                                       vars(ifield)%field(:))
                end do

            case ('EnstrophyEquation')
                call Diagnose_EnstrophyEquation(vars=vars)
                do ifield = 1, size(vars)
                    plot_file = 'EnstrophyEqn'//trim(adjustl(vars(ifield)%tag))//time_str(1:MaskSize)
                    call Write_Visuals(plot_file, vars(ifield)%field(:))
                end do

            case ('StrainTensor')
                call FI_STRAIN_TENSOR(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                call Write_Visuals(plot_file, txc(:, 1:6))

            case ('Strain 2S_ijS_ij (Log)')
                plot_file = 'LogStrain'//time_str(1:MaskSize)
                call FI_STRAIN(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3))
                txc(1:isize_field, 1) = 2.0_wp*txc(1:isize_field, 1)
                txc(1:isize_field, 1) = log10(txc(1:isize_field, 1) + small_wp)
                call Write_Visuals(plot_file, txc(:, 1:1))

            case ('StrainEquation')
                call Diagnose_StrainEquation(vars=vars)
                do ifield = 1, size(vars)
                    plot_file = 'StrainEqn'//trim(adjustl(vars(ifield)%tag))//time_str(1:MaskSize)
                    call Write_Visuals(plot_file, vars(ifield)%field(:))
                end do

            case ('VelocityGradientInvariants')
                plot_file = 'InvariantP'//time_str(1:MaskSize)
                call FI_INVARIANT_P(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2))
                call Write_Visuals(plot_file, txc(:, 1:1))

                plot_file = 'InvariantQ'//time_str(1:MaskSize)
                call FI_INVARIANT_Q(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), txc(1, 2), txc(1, 3), txc(1, 4))
                call Write_Visuals(plot_file, txc(:, 1:1))

                plot_file = 'InvariantR'//time_str(1:MaskSize)
                call FI_INVARIANT_R(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                                    txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5), txc(1, 6))
                call Write_Visuals(plot_file, txc(:, 1:1))

            case ('HorizontalDivergence')
                call OPR_Partial_X(OPR_P1, imax, jmax, kmax, q(1, 1), txc(1, 1))
                call OPR_Partial_Y(OPR_P1, imax, jmax, kmax, q(1, 2), txc(1, 2))
                txc(1:isize_field, 1) = txc(1:isize_field, 1) + txc(1:isize_field, 2)
                call Write_Visuals(plot_file, txc(:, 1:1))

            case ('Turbulent quantities')
                plot_file = 'Tke'//time_str(1:MaskSize)
                txc(1:isize_field, 1) = q(1:isize_field, 1); call FI_Fluctuation_InPlace(imax, jmax, kmax, txc(:, 1))
                txc(1:isize_field, 2) = q(1:isize_field, 2); call FI_Fluctuation_InPlace(imax, jmax, kmax, txc(:, 2))
                txc(1:isize_field, 3) = q(1:isize_field, 3); call FI_Fluctuation_InPlace(imax, jmax, kmax, txc(:, 3))
                txc(:, 4) = 0.5_wp*Tensor_Dot(txc(:, 1:3), txc(:, 1:3))
                call Write_Visuals(plot_file, txc(:, 4:4))

                plot_file = 'ReynoldsTensor'//time_str(1:MaskSize)
                txc(1:isize_field, 4) = txc(1:isize_field, 1)*txc(1:isize_field, 2)
                txc(1:isize_field, 5) = txc(1:isize_field, 1)*txc(1:isize_field, 3)
                txc(1:isize_field, 6) = txc(1:isize_field, 2)*txc(1:isize_field, 3)
                txc(1:isize_field, 1) = txc(1:isize_field, 1)*txc(1:isize_field, 1)
                txc(1:isize_field, 2) = txc(1:isize_field, 2)*txc(1:isize_field, 2)
                txc(1:isize_field, 3) = txc(1:isize_field, 3)*txc(1:isize_field, 3)
                call Write_Visuals(plot_file, txc(:, 1:6))

                ! plot_file = 'LogDissipation'//time_str(1:MaskSize)
                ! call FI_DISSIPATION(imax, jmax, kmax, q(1, 1), q(1, 2), q(1, 3), txc(1, 1), &
                !                     txc(1, 2), txc(1, 3), txc(1, 4), txc(1, 5))
                ! txc(1:isize_field, 1) = txc(1:isize_field, 1)*visc
                ! txc(1:isize_field, 1) = log10(txc(1:isize_field, 1) + small_wp)
                ! call Write_Visuals(plot_file, txc(:,1:1))

            case ('Atmospheric Thermodynamics')
                select case (imode_thermo)
                case (THERMO_TYPE_ANELASTIC)
                    call Diagnose_Energies_Anelastic(vars)
                    do ifield = 1, size(vars)
                        call Write_Visuals(trim(adjustl(vars(ifield)%tag))//time_str(1:MaskSize), &
                                           vars(ifield)%field(:))
                    end do

                    call Diagnose_Thetas_Anelastic(vars)
                    do ifield = 1, size(vars)
                        call Write_Visuals(trim(adjustl(vars(ifield)%tag))//time_str(1:MaskSize), &
                                           vars(ifield)%field(:))
                    end do

                    call Diagnose_Moisture_Anelastic(vars)
                    do ifield = 1, size(vars)
                        call Write_Visuals(trim(adjustl(vars(ifield)%tag))//time_str(1:MaskSize), &
                                           vars(ifield)%field(:))
                    end do

                case (THERMO_TYPE_COMPRESSIBLE)
                end select

                ! ###################################################################
                ! Radiation
                ! ###################################################################
            case ('Infrared')
                do is = 1, inb_scal
                    if (infraredProps%active(is)) then
                        write (str, *) is; plot_file = 'Infrared'//trim(adjustl(str))//time_str(1:MaskSize)
                        call Radiation_Infrared_Z(infraredProps, imax, jmax, kmax, s, &
                                                  txc(:, 1), txc(:, 2), txc(:, 3), txc(:, 4), txc(:, 5), txc(:, 6))
                        if (nse_eqns == DNS_EQNS_ANELASTIC) then
                            call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, ribackground, txc(:, 1))
                        end if
                        call Write_Visuals(plot_file, txc(:, 1:1))

                        write (str, *) is; plot_file = 'InfraredFlux'//trim(adjustl(str))//time_str(1:MaskSize)
                        txc(1:isize_field, 5) = txc(1:isize_field, 6) - txc(1:isize_field, 5)
                        call Write_Visuals(plot_file, txc(:, 5:5))

                        write (str, *) is; plot_file = 'InfraredFluxUp'//trim(adjustl(str))//time_str(1:MaskSize)
                        call Write_Visuals(plot_file, txc(:, 6:6))
                    end if

                end do

                ! ###################################################################
                ! Microphysics
                ! ###################################################################
            case ('Sedimentation')
                do is = 1, inb_scal
                    if (infraredProps%active(is)) then
                        write (str, *) is; plot_file = 'Sedimentation'//trim(adjustl(str))//time_str(1:MaskSize)
                        call Microphysics_Sedimentation_Z(sedimentationProps, imax, jmax, kmax, is, s, &
                                                          txc(:, 1), txc(:, 2), txc(:, 3))
                        if (nse_eqns == DNS_EQNS_ANELASTIC) then
                            call Thermo_Anelastic_Weight_InPlace(imax, jmax, kmax, ribackground, txc(:, 1))
                        end if
                        call Write_Visuals(plot_file, txc(:, 1:1))

                        write (str, *) is; plot_file = 'SedimentationFlux'//trim(adjustl(str))//time_str(1:MaskSize)
                        call Write_Visuals(plot_file, txc(:, 1:3))

                    end if

                end do

            end select

        end do

    end do

    call TLab_Stop(0)

contains
    ! ###################################################################
    ! ###################################################################
    subroutine Visuals_Initialize()
#ifdef USE_MPI
        use mpi_f08, only: MPI_COMM_WORLD, MPI_REAL4
        use TLabMPI_VARS
#endif
        character(len=32) bakfile, block
        character(len=128) eStr
        character(len=512) sRes

        ! -----------------------------------------------------------------------
#ifdef USE_MPI
        integer(wi) nz
#endif

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
        iv = iv + 1; opt_name(iv) = 'VelocityX'
        iv = iv + 1; opt_name(iv) = 'VelocityY'
        iv = iv + 1; opt_name(iv) = 'VelocityZ'
        iv = iv + 1; opt_name(iv) = 'VelocityVector'
        iv = iv + 1; opt_name(iv) = 'VelocityMagnitude'
        ! iv = iv + 1; opt_name(iv) = 'Density'   ! to be removed into thermodynamics
        ! iv = iv + 1; opt_name(iv) = 'Temperature'  ! to be removed into thermodynamics
        ! iv = iv + 1; opt_name(iv) = 'Pressure'
        iv = iv + 1; opt_name(iv) = 'Thermodynamics'
        iv = iv + 1; opt_name(iv) = 'Void'
        iv = iv + 1; opt_name(iv) = 'PressureAnalysis'
        iv = iv + 1; opt_name(iv) = 'Scalars'
        iv = iv + 1; opt_name(iv) = 'ScalarGradientVector'
        iv = iv + 1; opt_name(iv) = 'ScalarGradient G_iG_i (Log)'
        iv = iv + 1; opt_name(iv) = 'ScalarGradientEquation'
        iv = iv + 1; opt_name(iv) = 'VorticityVector'
        iv = iv + 1; opt_name(iv) = 'Enstrophy W_iW_i (Log)'
        iv = iv + 1; opt_name(iv) = 'EnstrophyEquation'
        iv = iv + 1; opt_name(iv) = 'StrainTensor'
        iv = iv + 1; opt_name(iv) = 'Strain 2S_ijS_ij (Log)'
        iv = iv + 1; opt_name(iv) = 'StrainEquation'
        iv = iv + 1; opt_name(iv) = 'VelocityGradientInvariants'
        iv = iv + 1; opt_name(iv) = 'HorizontalDivergence'
        iv = iv + 1; opt_name(iv) = 'Turbulent quantities'
        iv = iv + 1; opt_name(iv) = 'Buoyancy'
        iv = iv + 1; opt_name(iv) = 'Sedimentation'
        iv = iv + 1; opt_name(iv) = 'Infrared'
        iv = iv + 1; opt_name(iv) = 'Atmospheric Thermodynamics'
        iv = iv + 1; opt_name(iv) = 'Evaporation'
        if (iv > iopt_size_max) then ! Check
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Increase number of options.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
        end if

        ! 17, '. Relative humidity'
        ! 19, '. Analysis of B and V'
        ! 20, '. Total Stress Tensor'

        call ScanFile_Char(bakfile, ifile, block, 'ParamVisuals', '-1', sRes)
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
        iopt_size = iopt_size_max
        call LIST_INTEGER(sRes, iopt_size, opt_vec)

        if (opt_vec(1) < 0) then ! Check
            call TLab_Write_ASCII(efile, trim(adjustl(eStr))//'Missing ParamVisuals input.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
        end if

        ! -------------------------------------------------------------------
        call ScanFile_Char(bakfile, ifile, block, 'Subdomain', '-1', sRes)
        if (sRes == '-1') then
#ifdef USE_MPI
#else
            write (*, *) 'Subdomain limits (full domain, default)?'
            read (*, '(A64)') sRes
#endif
        end if
        is = 6
        call LIST_INTEGER(sRes, is, subdomain)

        if (is < 6) then ! default
            subdomain(1) = 1; subdomain(2) = x%size
            subdomain(3) = 1; subdomain(4) = y%size
            subdomain(5) = 1; subdomain(6) = z%size
        end if

        ! -------------------------------------------------------------------
        call ScanFile_Int(bakfile, ifile, block, 'Format', '-1', opt_format)
        if (opt_format == -1) then
#ifdef USE_MPI
#else
            write (*, *) 'File Format ?'
            write (*, *) ' 1. Raw, single precision, no header (default)'
            write (*, *) ' 2. General restart format'
            read (*, '(A64)') sRes
#endif
        end if
        if (len_trim(adjustl(sRes)) > 0) then
            if (trim(adjustl(sRes)) == 'single') then; opt_format = FORMAT_SINGLE
            else if (trim(adjustl(sRes)) == 'general') then; opt_format = FORMAT_GENERAL
            else
                read (sRes, *) opt_format
            end if
        end if

        if (opt_format < 0) opt_format = FORMAT_SINGLE ! default is single precission, no header

!     ! -------------------------------------------------------------------
!     ! Defining gate levels for conditioning
!     ! -------------------------------------------------------------------
!     opt_cond = 0 ! default values
!     opt_cond_relative = 0
!     igate_size = 0

!     do iv = 1, iopt_size
!         if (opt_vec(iv) == iscal_offset + 11) then
! #include "dns_read_partition.h"
!             if (opt_cond > 1) inb_txc = max(inb_txc, 5)
!             exit
!         end if
!     end do

        ! ###################################################################
        ! Initialization of array sizes
        ! ###################################################################
        iread_flow = .false.
        iread_scal = .false.
        inb_txc = 0

        do iv = 1, iopt_size
            select case (trim(adjustl(opt_name(opt_vec(iv)))))
            case ('VelocityX')
                iread_flow = .true.; inb_txc = max(inb_txc, 1)
            case ('VelocityY')
                iread_flow = .true.; inb_txc = max(inb_txc, 1)
            case ('VelocityZ')
                iread_flow = .true.; inb_txc = max(inb_txc, 1)
            case ('VelocityVector')
                iread_flow = .true.; inb_txc = max(inb_txc, 3)
            case ('VelocityMagnitude')
                iread_flow = .true.; inb_txc = max(inb_txc, 1)
                ! case ('Density')
                !     iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 1)
                ! case ('Temperature')
                !     iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 1)
            case ('Thermodynamics')
                iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 4)
            case ('PressureAnalysis')
                iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 9)
            case ('Scalars')
                iread_scal = .true.; inb_txc = max(inb_txc, 1)
            case ('ScalarGradientVector')
                iread_scal = .true.; inb_txc = max(inb_txc, 3)
            case ('ScalarGradient G_iG_i (Log)')
                iread_scal = .true.; inb_txc = max(inb_txc, 2)
            case ('ScalarGradientEquation')
                iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 7)
            case ('VorticityVector')
                iread_flow = .true.; inb_txc = max(inb_txc, 4)
            case ('Enstrophy W_iW_i (Log)')
                iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 7)
            case ('EnstrophyEquation')
                iread_flow = .true.; inb_txc = max(inb_txc, 8)
            case ('StrainTensor')
                iread_flow = .true.; inb_txc = max(inb_txc, 6)
            case ('Strain 2S_ijS_ij (Log)')
                iread_flow = .true.; inb_txc = max(inb_txc, 3)
            case ('StrainEquation')
                iread_flow = .true.; inb_txc = max(inb_txc, 8)
            case ('VelocityGradientInvariants')
                iread_flow = .true.; inb_txc = max(inb_txc, 6)
            case ('HorizontalDivergence')
                iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 2)
            case ('Turbulent quantities')
                iread_flow = .true.; inb_txc = max(inb_txc, 6)
            case ('Buoyancy')
                iread_flow = .true.; iread_scal = .true.; inb_txc = max(inb_txc, 5)
            case ('Atmospheric Thermodynamics')
                iread_scal = .true.; inb_txc = max(inb_txc, 3)
            case ('Infrared')
                iread_scal = .true.; inb_txc = max(inb_txc, 6)
            case ('Sedimentation')
                iread_scal = .true.; inb_txc = max(inb_txc, 3)
            end select
        end do

        ! #######################################################################
#ifdef USE_MPI

        io_subarrays(:)%active = .false.
        io_subarrays(:)%offset = 0
        io_subarrays(:)%precision = IO_TYPE_SINGLE

        nz = subdomain(6) - subdomain(5) + 1

        ! ###################################################################
        ! Saving full vertical xOz planes; using subdomain(3) to define the plane
        if (yMpi%rank == ((subdomain(3) - 1)/jmax)) io_subarrays(IO_SUBARRAY_VISUALS_XOZ)%active = .true.
        io_subarrays(IO_SUBARRAY_VISUALS_XOZ)%communicator = xMpi%comm
        io_subarrays(IO_SUBARRAY_VISUALS_XOZ)%subarray = IO_Create_Subarray_XOZ(imax, nz, MPI_REAL4)

        ! Saving full vertical yOz planes; using subiddomain(1) to define the plane
        if (xMpi%rank == ((subdomain(1) - 1)/imax)) io_subarrays(IO_SUBARRAY_VISUALS_YOZ)%active = .true.
        io_subarrays(IO_SUBARRAY_VISUALS_YOZ)%communicator = yMpi%comm
        io_subarrays(IO_SUBARRAY_VISUALS_YOZ)%subarray = IO_Create_Subarray_YOZ(jmax, nz, MPI_REAL4)

        ! Saving full blocks xOy planes
        io_subarrays(IO_SUBARRAY_VISUALS_XOY)%active = .true.
        io_subarrays(IO_SUBARRAY_VISUALS_XOY)%communicator = MPI_COMM_WORLD
        io_subarrays(IO_SUBARRAY_VISUALS_XOY)%subarray = IO_Create_Subarray_XOY(imax, jmax, nz, MPI_REAL4)

#else
        io_subarrays(:)%offset = 0
        io_subarrays(:)%precision = IO_TYPE_SINGLE

#endif

        return
    end subroutine Visuals_Initialize

!########################################################################
!########################################################################
    subroutine Write_Visuals_Rank2(fname, field)
        use Reductions, only: Reduce_Block_InPlace

        character(len=*) fname
        real(wp), intent(inout) :: field(:, :)

        ! -------------------------------------------------------------------
        integer(wi) nx, ny, nz, ifield, i
        integer(wi) sizes(5)
        character*32 varname(16)
        integer subarray_plan
        integer nfield

        ! ###################################################################
        ! I think this first block should go into Visuals_Initialize
        nfield = size(field, 2)

        nx = subdomain(2) - subdomain(1) + 1
        ny = subdomain(4) - subdomain(3) + 1
        nz = subdomain(6) - subdomain(5) + 1

        sizes(5) = nfield
        sizes(1) = size(field, 1)           ! array size
        sizes(2) = 1                        ! lower bound
        if (subdomain(2) - subdomain(1) + 1 == x%size .and. &
            subdomain(4) - subdomain(3) + 1 == 1) then              ! xOz plane
            subarray_plan = IO_SUBARRAY_VISUALS_XOZ
            sizes(3) = imax*nz              ! upper bound
            sizes(4) = 1                    ! stride
            nx = imax

        else if (subdomain(4) - subdomain(3) + 1 == y%size .and. &
                 subdomain(2) - subdomain(1) + 1 == 1) then         ! yOz plane
            subarray_plan = IO_SUBARRAY_VISUALS_YOZ
            sizes(3) = jmax*nz              ! upper bound
            sizes(4) = 1                    ! stride
            ny = jmax

        else if (subdomain(2) - subdomain(1) + 1 == x%size .and. &
                 subdomain(4) - subdomain(3) + 1 == y%size) then    ! xOy blocks
            subarray_plan = IO_SUBARRAY_VISUALS_XOY
            sizes(3) = imax*jmax*nz         ! upper bound
            sizes(4) = 1                    ! stride
            nx = imax
            ny = jmax

        else
#ifdef USE_MPI
            call TLab_Write_ASCII(efile, __FILE__//'. Invalid subdomain in parallel mode.')
            call TLab_Stop(DNS_ERROR_INVALOPT)
#else
            sizes(3) = nx*ny*nz             ! upper bound
            sizes(4) = 1                    ! stride
#endif

        end if

        ! ###################################################################
        select case (opt_format)
        case (FORMAT_GENERAL)
            if (nfield > 1) then ! IO_Write_Fields expects field to be aligned by isize_field (instead of isize_txc_field)
                do ifield = 2, nfield
                    do i = 1, isize_field
                        field((ifield - 1)*isize_field + i, 1) = field(i, ifield)
                    end do
                end do
            end if
            call IO_Write_Fields(fname, imax, jmax, kmax, itime, nfield, field)

        case (FORMAT_SINGLE)
            do ifield = 1, nfield
                call Reduce_Block_InPlace(imax, jmax, kmax, &
                                          mod(subdomain(1) - 1, imax) + 1, &    ! starting node
                                          mod(subdomain(3) - 1, jmax) + 1, &    ! starting node
                                          subdomain(5), &                       ! starting node
                                          nx, ny, nz, field(:, ifield))
            end do

            varname = ''
            if (nfield > 1) then
                do ifield = 1, nfield; write (varname(ifield), *) ifield; varname(ifield) = trim(adjustl(varname(ifield)))
                end do
            end if
            call IO_Write_Subarray(io_subarrays(subarray_plan), fname, varname, field, sizes)

        end select

        return
    end subroutine Write_Visuals_Rank2

    subroutine Write_Visuals_Rank1(fname, field)
        character(len=*) fname
        real(wp), intent(inout) :: field(:)

        target :: field
        real(wp), pointer :: field_2r(:, :) => null()

        field_2r(1:size(field), 1:1) => field(1:size(field))
        call Write_Visuals_Rank2(fname, field_2r)
        nullify (field_2r)

        return
    end subroutine

end program VISUALS
