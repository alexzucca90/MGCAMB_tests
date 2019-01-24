module MGCAMB_cache
    use precision

    Type :: MGCAMBFlags
        ! new model selection flags
        integer :: MG_flag
        integer :: pure_MG_flag
        integer :: alt_MG_flag
        integer :: QSA_flag
        integer :: mugamma_par
        integer :: muSigma_par
        integer :: QR_par

        ! DE model flag
        integer :: DE_model

        real(dl) :: GRtrans !< scale factor at which MG is switched on

    end Type MGCAMBFlags

    Type :: MGCAMBModelParams

        ! BZ parametrization (and QS f(R))
        real(dl) :: B1
        real(dl) :: B2
        real(dl) :: lambda1_2
        real(dl) :: lambda2_2
        real(dl) :: ss

        ! Planck Parametrization
        real(dl) :: E11
        real(dl) :: E22

        ! Q-R parametrization 1
        real(dl) :: MGQfix
        real(dl) :: MGRfix

        ! Q-R parametrization 2
        real(dl) :: Qnot
        real(dl) :: Rnot
        real(dl) :: sss

        ! Growth rate gamma
        real(dl) :: Linder_gamma

        ! Symmetron
        real(dl) :: beta_star
        real(dl) :: a_star
        real(dl) :: xi_star

        ! Dilaton
        real(dl) :: beta0
        real(dl) :: xi0
        real(dl) :: DilR
        real(dl) :: DilS

        ! Hu-Sawicki f(R) gravity
        real(dl) :: F_R0
        real(dl) :: FRn

        ! DES parametrization
        real(dl) :: mu0
        real(dl) :: sigma0

        ! DE model parameters
        real(dl) :: w0DE              !< w0 parameters for DE
        real(dl) :: waDE              !< waDE parameters for DE

    end Type MGCAMBModelParams

    character(len=(10)) :: MGCAMB_version = 'v 3.1'


    ! define the type MGCAMB_par_cache
    type :: MGCAMB_parameter_cache
        real(dl) :: omegab
        real(dl) :: omegac
        real(dl) :: omegav
        real(dl) :: h0
        real(dl) :: h0_Mpc
        character(len=30) :: output_root
    end type MGCAMB_parameter_cache

    type(MGCAMB_parameter_cache) :: mgcamb_par_cache

    ! define the tyoe MGCAMB_timestep_cache
    type :: MGCAMB_timestep_cache

        ! 1. Background quantities
        real(dl) :: adotoa
        real(dl) :: Hdot
        real(dl) :: grho
        real(dl) :: gpres
        real(dl) :: grhob_t
        real(dl) :: grhoc_t
        real(dl) :: grhog_t
        real(dl) :: grhor_t
        real(dl) :: grhov_t
        real(dl) :: gpresv_t
        real(dl) :: grhonu_t
        real(dl) :: gpresnu_t

        ! 2. Perturbation quantities
        real(dl) :: k
        real(dl) :: k2
        real(dl) :: dgrho
        real(dl) :: dgq
        real(dl) :: pidot_sum
        real(dl) :: dgpi_w_sum
        real(dl) :: dgpi
        real(dl) :: dgpi_diff
        real(dl) :: dgpidot
        real(dl) :: rhoDelta
        real(dl) :: rhoDeltadot

        ! 3. MG functions
        real(dl) :: mu
        real(dl) :: mudot
        real(dl) :: gamma
        real(dl) :: gammadot
        real(dl) :: q
        real(dl) :: qdot
        real(dl) :: r
        real(dl) :: rdot

        !> 4. Perturbations evolution variables
        real(dl) :: z
        real(dl) :: sigma
        real(dl) :: sigmadot
        real(dl) :: etak
        real(dl) :: etadot

        !> 5. ISW and lensing realted quantities
        real(dl) :: MG_alpha
        real(dl) :: MG_alphadot
        real(dl) :: MG_phi
        real(dl) :: MG_phidot
        real(dl) :: MG_psi
        real(dl) :: MG_psidot
        real(dl) :: MG_ISW
        real(dl) :: MG_lensing
        real(dl) :: source1
        real(dl) :: source3


    end type MGCAMB_timestep_cache

#ifdef DEBUG
    logical , parameter :: DebugMGCAMB = .true.              !< MGCAMB debug flag.This will turn on printing of many things to aid debugging the code.
#else
    logical , parameter :: DebugMGCAMB = .false.             !< MGCAMB debug flag.This will turn on printing of many things to aid debugging the code.
#endif


contains

    !---------------------------------------------------------------------------------------------
    !> Subroutine that reads the MGCAMB model parameters
    subroutine MGCAMB_read_model_params( mg_par_cache, mg_flags, mg_model_pars )
        use IniFile
        Type(MGCAMB_parameter_cache), intent(in)    :: mg_par_cache  !< cache containing the parameters
        Type(MGCAMBFlags), intent(out)              :: mg_flags         !< containing the MGCAMB flags
        Type(MGCAMBModelParams), intent(out)        :: mg_model_pars    !< containing the model parameters

        ! 1. MG_flag
        mg_flags%MG_flag = Ini_Read_Int('MG_flag', 0)

        if ( mg_flags%MG_flag /= 0 ) then
            call print_MGCAMB_header
            write(*,*)
            write(*,*) 'MG_flag:', mg_flags%MG_flag

            write(*,*) 'Debug:', DebugMGCAMB


            ! read GRtrans
            mg_flags%GRtrans = Ini_Read_Double('GRtrans',0.01_dl)
            write(*,*) '    GRtrans:', mg_flags%GRtrans

            ! 1. pure MG models
            if ( mg_flags%MG_flag == 1 ) then

                mg_flags%pure_MG_flag = Ini_Read_Int('pure_MG_flag', 1)

                if ( mg_flags%pure_MG_flag == 1 ) then ! mu-gamma
                    write(*,*) '    MGCAMB: mu-gamma parametrization'
                    mg_flags%mugamma_par = Ini_Read_Int('mugamma_par' , 1)
                    if ( mg_flags%mugamma_par == 1 ) then
                        write(*,*) '        BZ parametrization'
                        mg_model_pars%B1= Ini_Read_Double('B1',0._dl)
                        mg_model_pars%B2= Ini_Read_Double('B2',0._dl)
                        mg_model_pars%lambda1_2= Ini_Read_Double('lambda1_2',0._dl)
                        mg_model_pars%lambda2_2= Ini_Read_Double('lambda2_2',0._dl)
                        mg_model_pars%ss= Ini_Read_Double('ss',0._dl)
                    else if ( mg_flags%mugamma_par == 2 ) then
                        write(*,*) '        Planck parametrization'
                        mg_model_pars%E11     = Ini_Read_Double('E11', 0._dl)
                        mg_model_pars%E22     = Ini_Read_Double('E22', 0._dl)
                        write(*,*) 'E11, E22', mg_model_pars%E11, mg_model_pars%E22
                    else
                        write(*,*) ' write your own mu-gamma parametrization in mgcamb.f90'
                        stop
                    end if


                else if ( mg_flags%pure_MG_flag == 2 ) then ! mu-Sigma
                    write(*,*) '    MGCAMB: mu-Sigma parametrization'
                    mg_flags%muSigma_par = Ini_Read_Int('musigma_par', 1)
                    if ( mg_flags%muSigma_par == 1 ) then
                        write(*,*) '        DES parametrization'
                        mg_model_pars%mu0     = Ini_Read_Double('mu0', 0._dl)
                        mg_model_pars%sigma0  = Ini_Read_Double('sigma0', 0._dl)
                        write(*,*) 'mu0, sigma0:', mg_model_pars%mu0, mg_model_pars%sigma0
                    else if ( mg_flags%muSigma_par == 2 ) then
                        write(*,*) 'write you own mu-sigma parametrization in mgcamb.f90'
                        stop
                    else
                        write(*,*) 'Please choose a model in params_MG.ini'
                        stop
                    end if

                else if ( mg_flags%pure_MG_flag == 3 ) then ! Q-R
                    write(*,*) '    MGCAMB: Q-R parametrization'
                    mg_flags%QR_par = Ini_Read_Int('QR_par', 1)
                    if ( mg_flags%QR_par == 1 ) then
                        mg_model_pars%MGQfix=Ini_Read_Double('MGQfix', 0._dl)
                        mg_model_pars%MGRfix=Ini_Read_Double('MGRfix', 0._dl)
                    else if ( mg_flags%QR_par == 2 ) then
                        mg_model_pars%Qnot=Ini_Read_Double('Qnot', 0._dl)
                        mg_model_pars%Rnot=Ini_Read_Double('Rnot', 0._dl)
                        mg_model_pars%sss=Ini_Read_Double('sss', 0._dl)
                    else if ( mg_flags%QR_par == 3 ) then
                        write(*,*) 'write your own QR parametrization in mgcamb.f90'
                        stop
                    else
                        write(*,*) 'Please choose a model in params_MG.ini'
                        stop
                    end if

                end if

                ! Checking DE Model
                mg_flags%DE_model = Ini_Read_Int('DE_model', 0)

                write(*,*) 'DE_model:', mg_flags%DE_model

                if ( mg_flags%DE_model == 1 ) then
                    mg_model_pars%w0DE = Ini_Read_Double('w0DE', -1._dl)
                else if ( mg_flags%DE_model == 2 ) then
                    mg_model_pars%w0DE = Ini_Read_Double('w0DE', -1._dl)
                    mg_model_pars%waDE = Ini_Read_Double('waDE', 0._dl)
                else if ( mg_flags%DE_model == 3 ) then
                    write(*,*) 'This will contain the reconstruction of w_DE(a)'
                    write(*,*) 'Not implemented yet'
                    stop
                else if ( mg_flags%DE_model == 4 ) then
                    write(*,*) 'This will contain the reconstruction of rho_DE(a)'
                    write(*,*) 'Not implemented yet'
                    stop
                else if ( mg_flags%DE_model /= 0 ) then
                    write(*,*) 'Please choose a DE model'
                    stop
                end if



            else if ( mg_flags%MG_flag == 2 ) then
                mg_flags%alt_MG_flag = Ini_Read_Int('alt_MG_flag', 1)
                if ( mg_flags%alt_MG_flag == 1 ) then
                    write(*,*) '    MGCAMB: Linder Gamma'
                    mg_model_pars%Linder_gamma = Ini_Read_Double('Linder_gamma', 0._dl)
                else if ( mg_flags%alt_MG_flag == 2 ) then
                    write(*,*) 'Please write your alternative MG model in mgcamb.f90'
                    stop
                else
                    write(*,*) 'Please choose a model in params_MG.ini'
                    stop
                end if

                ! Checking DE Model
                mg_flags%DE_model = Ini_Read_Int('DE_model', 0)

                if ( mg_flags%DE_model /= 0 ) then
                    write(*,*) 'alternative MG models supported only with cosmological constant!'
                end if


            else if ( mg_flags%MG_flag == 3 ) then
                write(*,*) '    MGCAMB: quasi-static models'
                mg_flags%QSA_flag = Ini_Read_Int('QSA_flag', 1)
                if ( mg_flags%QSA_flag ==  1 ) then
                    write(*,*) '        QSA f(R)'
                    mg_model_pars%B1 = 4._dl/3._dl
                    mg_model_pars%lambda1_2= Ini_Read_Double('B0',0._dl) ! it is considered as the B0 parameter here
                    mg_model_pars%lambda1_2 = (lambda1_2*(299792458.d-3)**2)/(2._dl*mg_par_cache%H0**2)
                    mg_model_pars%B2 = 0.5d0
                    mg_model_pars%lambda2_2 = B1* lambda1_2
                    mg_model_pars%ss = 4._dl

                else if ( mg_flags%QSA_flag ==  2 ) then
                    write(*,*) '        QSA Symmetron'
                    mg_model_pars%beta_star = Ini_Read_Double('beta_star', 0._dl)
                    mg_model_pars%xi_star = Ini_Read_Double ('xi_star', 0._dl)
                    mg_model_pars%a_star = Ini_Read_Double('a_star', 0._dl)
                    mg_flags%GRtrans = mg_model_pars%a_star

                else if ( mg_flags%QSA_flag ==  3 ) then
                    write(*,*) '        QSA Dilaton'
                    ! GENERALIZED DILATON
                    mg_model_pars%beta0 = Ini_Read_Double('beta0', 0._dl)
                    mg_model_pars%xi0 = Ini_Read_Double('xi0', 0._dl)
                    mg_model_pars%DilR = Ini_Read_Double('DilR', 0._dl)
                    mg_model_pars%DilS = Ini_Read_Double('DilS', 0._dl)

                else if ( mg_flags%QSA_flag ==  4 ) then
                    write(*,*) '        QSA Hu-Sawicki f(R)'
                    mg_model_pars%F_R0 = Ini_Read_Double('F_R0', 0._dl)
                    mg_model_pars%FRn = Ini_Read_Double('FRn', 0._dl)
                    mg_model_pars%beta0 = 1._dl/sqrt(6._dl)
                else if ( mg_flags%QSA_flag ==  5 ) then
                    write(*,*) 'Please write your QSA model in mgcamb.f90'
                    stop

                end if

                ! Checking DE Model
                mg_flags%DE_model = Ini_Read_Int('DE_model', 0)

                if ( mg_flags%DE_model /= 0 ) then
                    write(*,*) 'QSA models supported only with cosmological constant!'
                end if

            else
                write(*,*) ' Please choose a model'
                stop
            end if


        end if

    end subroutine MGCAMB_read_model_params

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the MGCAMB header.
    subroutine print_MGCAMB_header

        implicit none

        ! print the header:
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') "     __  _________  ________   __  ______  "
        write(*,'(a)') "    /  \/  / ____/ / ___/ _ | /  |/  / _ ) "
        write(*,'(a)') "   / /\_/ / /_,-, / /__/ __ |/ /|_/ / _  | "
        write(*,'(a)') "  /_/  /_/_____/  \___/_/ |_/_/  /_/____/  "//" "//MGCAMB_version
        write(*,'(a)') "  "
        write(*,'(a)') "        Modified Growth with CAMB "
        write(*,'(a)') "  "
        write(*,'(a)') "***************************************************************"

    end subroutine print_MGCAMB_header

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the mgcamb_cache to zero
    subroutine MGCAMB_timestep_cache_nullify( mg_cache )
        use precision
        implicit none

        type(MGCAMB_timestep_cache),  intent(inout) :: mg_cache      !< cache containing the time-dependent quantities

        ! 1. Background quantities
        mg_cache%adotoa     = 0._dl
        mg_cache%Hdot       = 0._dl
        mg_cache%grho       = 0._dl
        mg_cache%gpres      = 0._dl
        mg_cache%grhob_t    = 0._dl
        mg_cache%grhoc_t    = 0._dl
        mg_cache%grhog_t    = 0._dl
        mg_cache%grhor_t    = 0._dl
        mg_cache%grhov_t    = 0._dl
        mg_cache%gpresv_t   = 0._dl
        mg_cache%grhonu_t   = 0._dl
        mg_cache%gpresnu_t  = 0._dl

        ! 2. Perturbation quantities
        mg_cache%k          = 0._dl
        mg_cache%k2         = 0._dl
        mg_cache%dgrho      = 0._dl
        mg_cache%dgq        = 0._dl
        mg_cache%pidot_sum  = 0._dl
        mg_cache%dgpi_w_sum = 0._dl
        mg_cache%dgpi       = 0._dl
        mg_cache%dgpi_diff  = 0._dl
        mg_cache%dgpidot    = 0._dl
        mg_cache%rhoDelta   = 0._dl
        mg_cache%rhoDeltadot= 0._dl

        ! 3. MG functions
        mg_cache%mu         = 0._dl
        mg_cache%mudot      = 0._dl
        mg_cache%gamma      = 0._dl
        mg_cache%gammadot   = 0._dl
        mg_cache%q          = 0._dl
        mg_cache%qdot       = 0._dl
        mg_cache%r          = 0._dl
        mg_cache%rdot       = 0._dl

        !> 4. Perturbations evolution variables
        mg_cache%z          = 0._dl
        mg_cache%sigma      = 0._dl
        mg_cache%sigmadot   = 0._dl
        mg_cache%etak       = 0._dl
        mg_cache%etadot     = 0._dl

        !> 5. ISW and lensing realted quantities
        mg_cache%MG_alpha   = 0._dl
        mg_cache%MG_alphadot= 0._dl
        mg_cache%MG_phi     = 0._dl
        mg_cache%MG_phidot  = 0._dl
        mg_cache%MG_psi     = 0._dl
        mg_cache%MG_psidot  = 0._dl
        mg_cache%MG_ISW     = 0._dl
        mg_cache%MG_lensing = 0._dl
        mg_cache%source1    = 0._dl
        mg_cache%source3    = 0._dl

    end subroutine

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that opens the MGCAMB cache files (for Debug)
    subroutine MGCAMB_open_cache_files
        use precision
        implicit none

        ! 1. Open sources file
        open(unit=111, file=TRIM(mgcamb_par_cache%output_root) // 'MGCAMB_debug_sources.dat', status="new", &
            & action="write")
        write(111,*)  'k  ', 'a  ', 'MG_ISW  ', 'MG_Lensing  ', 'S_T  ', 'S_lensing'

        ! 2 Open MG functions file
        open(unit=222, file=TRIM(mgcamb_par_cache%output_root) // 'MGCAMB_debug_MG_fncs.dat', status="new",&
            & action="write")
        write(222,*)  'k  ', 'a  ', 'mu  ', 'gamma ', 'Q ', 'R ', 'Phi ', 'Psi ', 'dPhi ', 'dPsi '

        ! 3. Open Einstein solutions file
        open(unit=333, file=TRIM(mgcamb_par_cache%output_root) // 'MGCAMB_debug_EinsteinSol.dat', status="new",&
            & action="write")
        write(333,*) 'k  ', 'a  ', 'etak  ', 'z  ', 'sigma  ', 'etadot  ', 'sigmadot  '

        ! 4. Open Perturbation solution file
        open(unit=444, file=TRIM(mgcamb_par_cache%output_root) // 'MGCAMB_debug_PerturbSol.dat', status="new",&
        & action="write")
        write(444,*)  'k  ', 'a  ', 'dgrho  ', 'dgq  ', 'rhoDelta  ', 'dgpi  ', 'pidot_sum  ', 'dgpi_w_sum  '

        ! 5. Open Background file
        open(unit=555, file=TRIM(mgcamb_par_cache%output_root) // 'MGCAMB_debug_Background.dat', status="new",&
            & action="write")
        write(555,*)  'k  ', 'a  ', 'H  ', 'Hdot  ', 'grhov_t  ', 'gpresv_t  '

    end subroutine MGCAMB_open_cache_files


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that closes the MGCAMB cache files (for Debug)
    subroutine MGCAMB_close_cache_files
        use precision
        implicit none

        close(111);close(222); close(333);close(444);close(555)

    end subroutine MGCAMB_close_cache_files


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints the MGCAMB cache on a file
    subroutine MGCAMB_dump_cache( a, mg_cache )
        use precision
        implicit none

        real(dl), intent(in) :: a   !< scale factor
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache      !< cache containing the time-dependent quantities
        character(*), parameter :: cache_output_format = 'e18.8'


        ! 1. Write the sources
        write(111,'(14'//cache_output_format//')') mg_cache%k, a, mg_cache%MG_ISW, mg_cache%MG_Lensing,&
                                                    & mg_cache%source1, mg_cache%source3

        ! 2. Write the MG functions and the potentials
        write(222,'(14'//cache_output_format//')') mg_cache%k, a, mg_cache%mu, mg_cache%gamma, mg_cache%q, mg_cache%r, &
                                                & mg_cache%MG_phi, mg_cache%MG_psi, mg_cache%MG_phidot, mg_cache%MG_psidot

        ! 3. Write the Einstein equations solutions
        write(333,'(14'//cache_output_format//')') mg_cache%k, a, mg_cache%etak, mg_cache%z, mg_cache%sigma,&
                                                & mg_cache%etadot,mg_cache%sigmadot

        ! 4. Write the Perturbations Solutions
        write(444,'(14'//cache_output_format//')') mg_cache%k, a, mg_cache%dgrho, mg_cache%dgq, mg_cache%rhoDelta,&
                                                    & mg_cache%dgpi, mg_cache%pidot_sum, mg_cache%dgpi_w_sum

        !5. Write the background
        write(555,'(14'//cache_output_format//')') mg_cache%k, a, mg_cache%adotoa, mg_cache%Hdot, mg_cache%grhov_t,&
                                                    & mg_cache%gpresv_t



    end subroutine MGCAMB_dump_cache



end module MGCAMB_cache

