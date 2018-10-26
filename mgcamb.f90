module MGCAMB
    use precision

    integer :: model
    integer :: DE_model

    real(dl) :: GRtrans                     !< scale factor at which MG is switched on
    real(dl) :: B1, B2, lambda1_2, lambda2_2, ss
    real(dl) :: MGQfix, MGRfix, Qnot, Rnot, sss
    real(dl) :: Linder_gamma                !< for model 6 (Linder Gamma)
    real(dl) :: beta_star, a_star, xi_star  !< for model 7 (symmetron)
    real(dl) :: beta0, xi0, DilR, DilS, A_2 !< for model 8 and 10 (dilaton)
    real(dl) :: F_R0, FRn                   !< for model 9 (large curvature f(R))

    real(dl) :: wDE             !< constant wDE
    real(dl) :: w0, wa          !< w0,wa parameters for DE

    character(len=(10)) :: MGCAMB_version = 'v 3.0'


    ! define the type MGCAMB_par_cache
    type :: MGCAMB_parameter_cache
        real(dl) :: omegab
        real(dl) :: omegac
        real(dl) :: omegav
        real(dl) :: h0
        real(dl) :: h0_Mpc
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

    end type MGCAMB_timestep_cache

contains

    !---------------------------------------------------------------------------
    !> this subroutine computes the MG functions at a time-step
    subroutine MGCAMB_compute_MG_functions( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache), intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        if ( model == 1 .or. &
            model == 4 .or. &
            model == 5 .or. &
            model == 6 .or. &
            model == 7 .or. &
            model == 8 .or. &
            model == 9 .or. &
            model == 10) then

            mg_cache%mu         = MGCAMB_Mu( a, mg_par_cache, mg_cache )
            mg_cache%mudot      = MGCAMB_MuDot( a, mg_par_cache, mg_cache )
            mg_cache%gamma      = MGCAMB_Gamma( a, mg_par_cache, mg_cache )
            mg_cache%gammadot   = MGCAMB_GammaDot( a, mg_par_cache, mg_cache )

            ! other EFT functions are zero
            mg_cache%q      = 0.d0
            mg_cache%qdot   = 0.d0
            mg_cache%r      = 0.d0
            mg_cache%rdot   = 0.d0

        else if ( model == 2 .or. &
                  model == 3 ) then

            mg_cache%q      = MGCAMB_Q( a, mg_par_cache, mg_cache )
            mg_cache%qdot   = MGCAMB_Qdot( a, mg_par_cache, mg_cache )
            mg_cache%r      = MGCAMB_R( a, mg_par_cache, mg_cache )
            mg_cache%rdot   = MGCAMB_Rdot( a, mg_par_cache, mg_cache )

            ! other MG functions are zero
            mg_cache%mu         = 0.d0
            mg_cache%mudot      = 0.d0
            mg_cache%gamma      = 0.d0
            mg_cache%gammadot   = 0.d0

        end if

    end subroutine MGCAMB_compute_MG_functions


    !---------------------------------------------------------------------------
    !> this subroutine computes the shear sigma in MG
    subroutine MGCAMB_compute_sigma( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache), intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        if ( model == 1 .or. &
            model == 4 .or. &
            model == 5 .or. &
            model == 6 .or. &
            model == 7 .or. &
            model == 8 .or. &
            model == 9 .or. &
            model == 10) then

            ! first calculate MG_alpha
            mg_cache%MG_alpha = ( mg_cache%etak/mg_cache%k + mg_cache%mu * ( mg_cache%gamma*mg_cache%rhoDelta+ &
                                ( mg_cache%gamma -1.d0 )*2.d0* mg_cache%dgpi)/(2.d0*mg_cache%k2)) / mg_cache%adotoa

            ! then calculate sigma
            mg_cache%sigma = mg_cache%k * mg_cache%MG_alpha

        else if ( model == 2 .or. &
                  model == 3 ) then

            mg_cache%MG_phi      = - mg_cache%rhoDelta * mg_cache%q/(2.d0*mg_cache%k2)
            mg_cache%sigma       = (mg_cache%etak - mg_cache%k * mg_cache%MG_phi)/mg_cache%adotoa
            mg_cache%MG_alpha    = mg_cache%sigma/mg_cache%k

        end if

    end subroutine MGCAMB_compute_sigma

    !---------------------------------------------------------------------------
    !> this subroutine computes the shear sigma in MG
    subroutine MGCAMB_compute_z( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache), intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        !> other parameters
        real(dl) :: fmu
        real(dl) :: f1
        real(dl) :: fQ
        real(dl) :: term1
        real(dl) :: term2
        real(dl) :: term3
        real(dl) :: term4
        real(dl) :: term5
        real(dl) :: term6
        real(dl) :: k2alpha

        if ( model == 1 .or. &
            model == 4 .or. &
            model == 5 .or. &
            model == 6 .or. &
            model == 7 .or. &
            model == 8 .or. &
            model == 9 .or. &
            model == 10) then

            !> adding the massive neutrinos contibutions
            fmu = mg_cache%k2+0.5d0*mg_cache%gamma*mg_cache%mu*(3.d0*(mg_cache%grhoc_t+mg_cache%grhob_t) &
                & + 4.d0*(mg_cache%grhog_t+mg_cache%grhor_t) +3.d0 * (mg_cache%grhonu_t + mg_cache%gpresnu_t ))

            !> adding massive neutrinos contributions

            !f1 = mg_cache%k2+3.d0*( mg_cache%adotoa**2 - mg_cache%Hdot )
            f1 = mg_cache%k2+0.5d0*(3.d0*(mg_cache%grhoc_t+mg_cache%grhob_t) &
                & + 4.d0*(mg_cache%grhog_t+mg_cache%grhor_t) + 3.d0*(mg_cache%grhonu_t+mg_cache%gpresnu_t) &
                & + 3.d0*(mg_cache%grhov_t+mg_cache%gpresv_t))

            term1 = mg_cache%gamma*mg_cache%mu* f1 * mg_cache%dgq/mg_cache%k

            !> adding massive neutrinos contribution, if w_DE /= -1 this has to be changed
            term2 = mg_cache%k2*mg_cache%MG_alpha* ((mg_cache%mu* mg_cache%gamma- 1.d0)*(mg_cache%grhoc_t+mg_cache%grhob_t&
                    & +(4.d0/3.d0)*(mg_cache%grhog_t+mg_cache%grhor_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t)) &
                    & - (mg_cache%grhov_t+ mg_cache%gpresv_t))

            !term2 = mg_cache%k2*mg_cache%MG_alpha* (mg_cache%mu* mg_cache%gamma*( mg_cache%grhoc_t+mg_cache%grhob_t&
            !        & +(4.d0/3.d0)*(mg_cache%grhog_t+mg_cache%grhor_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t) ) &
            !        & - 2.d0*(mg_cache%adotoa**2 - mg_cache%Hdot))

            term3= (mg_cache%mu * ( mg_cache%gamma -1.d0)* mg_cache%adotoa - mg_cache%gamma*mg_cache%mudot &
                    & - mg_cache%gammadot*mg_cache%mu )*mg_cache%rhoDelta

            ! typo corrected here
            term4 = 2.d0*mg_cache%mu*(mg_cache%gamma - 1.d0)*mg_cache%adotoa*mg_cache%dgpi_w_sum

            ! separated fromt the previous term
            term5 = -2.d0*((mg_cache%gamma-1.d0)*mg_cache%mudot -mg_cache%gammadot*mg_cache%mu)*mg_cache%dgpi

            !> adding massive neutrinos contribution
            term6= 2.d0 * mg_cache%mu*(1.d0 - mg_cache%gamma)* mg_cache%pidot_sum

            !> calculate etadot
            mg_cache%etadot = (term1 + term2 + term3 + term4 + term5 + term6)/( 2.d0 * fmu)

            !> finally calculate Z
            mg_cache%z = mg_cache%sigma - 3.d0 * mg_cache%etadot/mg_cache%k

            !> Calculate the Newtonian potential
            mg_cache%MG_psi = - mg_cache%mu * ( mg_cache%rhoDelta + 2.d0* mg_cache%dgpi)/(2.d0*mg_cache%k2)

            !> calculate the curvature perturbation potential
            mg_cache%MG_phi = mg_cache%gamma * mg_cache%MG_psi + mg_cache%mu* 1.d0*mg_cache%dgpi/mg_cache%k2

            mg_cache%MG_phidot = mg_cache%etadot - mg_cache%adotoa * (mg_cache%MG_psi - mg_cache%adotoa * mg_cache%MG_alpha) &
                                & - mg_cache%Hdot * mg_cache%MG_alpha

        else if ( model == 2 .or. &
                  model == 3 ) then

            ! adding massive neutrinos contributions
            fQ = mg_cache%k2 + 0.5d0*mg_cache%q * (3.d0*(mg_cache%grhob_t+mg_cache%grhoc_t)+&
                & 4.d0*(mg_cache%grhor_t+mg_cache%grhog_t)+3.d0*(mg_cache%grhonu_t + mg_cache%gpresnu_t))

            ! fixed for w_DE /= -1
            !f1=mg_cache%k2+3.d0*( mg_cache%adotoa**2 - mg_cache%Hdot )
            f1 = mg_cache%k2+0.5d0*(3.d0*(mg_cache%grhoc_t+mg_cache%grhob_t) &
                & + 4.d0*(mg_cache%grhog_t+mg_cache%grhor_t) + 3.d0*(mg_cache%grhonu_t+mg_cache%gpresnu_t) &
                & + 3.d0*(mg_cache%grhov_t+mg_cache%gpresv_t))

            k2alpha= mg_cache%k * mg_cache%sigma

            term1 = mg_cache%q * f1 * mg_cache%dgq/mg_cache%k

            term2 = k2alpha * ((mg_cache%q - 1.d0) * ( mg_cache%grhob_t+mg_cache%grhoc_t+(4.d0/3.d0) &
                    & *(mg_cache%grhor_t+mg_cache%grhog_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t) &
                    & ) -mg_cache%grhov_t - mg_cache%gpresv_t)

            !term2 = k2alpha * ((mg_cache%q) * ( mg_cache%grhob_t+mg_cache%grhoc_t+(4.d0/3.d0) &
            !    & *(mg_cache%grhor_t+mg_cache%grhog_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t)) &
            !    & - 2.d0 *(mg_cache%adotoa**2 - mg_cache%Hdot))

            term3 = -( mg_cache%qdot + (mg_cache%r-1.d0) * mg_cache%q * mg_cache%adotoa ) * mg_cache%rhoDelta

            mg_cache%etadot = (term1 + term2 + term3)/( 2.d0 * fQ )

            mg_cache%z = mg_cache%sigma - 3.d0 * mg_cache%etadot/mg_cache%k

            !calculating also ISW related quantities
            mg_cache%MG_psi     = mg_cache%r * mg_cache%MG_phi - mg_cache%q * 1.d0 * mg_cache%dgpi/mg_cache%k2
            mg_cache%MG_phidot  = mg_cache%etadot - mg_cache%adotoa * (mg_cache%MG_psi - mg_cache%adotoa * mg_cache%MG_alpha) &
                                & - mg_cache%Hdot * mg_cache%MG_alpha

        end if

        ! calculate sigmadot
        mg_cache%sigmadot = mg_cache%k * (mg_cache%MG_psi - mg_cache%adotoa * mg_cache%MG_alpha)

    end subroutine MGCAMB_compute_z

    !---------------------------------------------------------------------------
    !> this subroutine computes the ISW term in MG
    subroutine MGCAMB_compute_ISW( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache), intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        !local variables
        real(dl) :: term0

        term0 = mg_cache%k2 + 3.d0* (mg_cache%adotoa**2.d0 - mg_cache%Hdot)

        !adding MG_rhoDeltadot
        mg_cache%rhoDeltadot = -term0 * mg_cache%dgq/mg_cache%k - (mg_cache%grho + mg_cache%gpres)* mg_cache%k*mg_cache%z &
                            & - mg_cache%adotoa * mg_cache%rhoDelta - 2.d0 * mg_cache%adotoa * mg_cache%dgpi

        !adding dgpidot
        mg_cache%dgpidot = mg_cache%pidot_sum - (2.d0*mg_cache%dgpi+ mg_cache%dgpi_diff )*mg_cache%adotoa

        if (model==1 .or. &
            model==4 .or. &
            model==5 .or. &
            model==6 .or. &
            model==7 .or. &
            model==8 .or. &
            model==9 .or. &
            model==10) then

            mg_cache%MG_psidot = - 0.5d0*mg_cache%mu/mg_cache%k2*(mg_cache%rhoDeltadot+2.d0*mg_cache%dgpidot) &
                                & - 0.5d0*mg_cache%mudot/mg_cache%k2*(mg_cache%rhoDelta+2.d0*mg_cache%dgpi)

        else if (model==2 .or. &
                 model==3) then

            mg_cache%MG_psidot = mg_cache%R * mg_cache%MG_phidot + mg_cache%Rdot * mg_cache%MG_phi - &
                            & mg_cache%Qdot*mg_cache%dgpi/mg_cache%k2 - mg_cache%Q * mg_cache%dgpidot /mg_cache%k2

        end if

        mg_cache%MG_ISW = mg_cache%MG_phidot+mg_cache%MG_psidot

        mg_cache%MG_alphadot = mg_cache%MG_psi - mg_cache%adotoa * mg_cache%MG_alpha

    end subroutine MGCAMB_compute_ISW

    !---------------------------------------------------------------------------
    !> this subroutine computes the lensing term in MG
    subroutine MGCAMB_compute_lensing( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache), intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        mg_cache%MG_lensing = mg_cache%MG_phi + mg_cache%MG_psi

    end subroutine MGCAMB_compute_lensing


    !-----------------------------------------------
    !> mu(a,k) function
    function MGCAMB_Mu( a, mg_par_cache, mg_cache )
        !use ModelParams
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Mu                                       !< MG mu function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        ! local variables
        real(dl) :: LKA1 ! \lambda_1^2 k^2 a^s
        real(dl) :: LKA2 ! \lambda_1^2 k^2 a^s
        real(dl) :: t1, t2, t1dot, t2dot
        real(dl) :: omm, ommdot

        ! beta, m parametrization
        real(dl) :: beta, m

        if( model==1 .or. &     !< usual mu-gamma parametrization
            model==4 .or. &
            model==5 ) then

            LKA1 = lambda1_2 * mg_cache%k2 * a**ss
            LKA2 = lambda2_2 * mg_cache%k2 * a**ss

            MGCAMB_Mu = (1.d0 + B1 * LKA1)/(1.d0 + LKA1)

            if (model ==4) then ! correction for f(R) mu function.
                MGCAMB_Mu = MGCAMB_Mu/(1.d0 - 1.4d-8 * lambda1_2 * a**3)
            end if


        else if (model ==6) then ! Linder Gamma. This has to be adapted to arbitrary background
                omm=(mg_par_cache%omegab+mg_par_cache%omegac)/((mg_par_cache%omegab+mg_par_cache%omegac) &
                    & + (1-mg_par_cache%omegab-mg_par_cache%omegac)*a**3)
                ommdot=-3.d0*omm**2*a**3*mg_cache%adotoa*(1-mg_par_cache%omegab-mg_par_cache%omegac) &
                        & /(mg_par_cache%omegab+mg_par_cache%omegac)

                MGCAMB_Mu=2.d0/3.d0*omm**(Linder_gamma-1.d0)*&
                (omm**Linder_gamma+2-3.d0*Linder_gamma+3.d0*(Linder_gamma-0.5d0)*omm)


        else if (model == 7 .or. &      !< mapping beta,mu into mu,gamma
                model==8 .or. &
                model==9 .or. &
                model==10) then

            beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
            m       = MGCAMB_M( a, mg_par_cache, mg_cache )
            t1      = (2.d0*beta**2.d0)*mg_cache%k2
            t2      = (m**2.d0)*a**2.d0

            MGCAMB_Mu = (mg_cache%k2 + t1 + t2)/(mg_cache%k2 + t2)

        end if
    end function MGCAMB_Mu

    !-----------------------------------------------
    !> \dot{mu}(a,k) function
    function MGCAMB_Mudot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Mudot                                    !< MG mudot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        ! local variables
        real(dl) :: LKA1 ! \lambda_1^2 k^2 a^s
        real(dl) :: LKA2 ! \lambda_1^2 k^2 a^s
        real(dl) :: k2, t1,t2,t1dot,t2dot
        real(dl) :: omm, ommdot

        ! mapping beta,m into mu,gamma
        real(dl) :: beta, betadot, m, mdot
        real(dl) :: mu

        if( model==1 .or. &
            model==4 .or. &
            model==5) then

            LKA1 = lambda1_2 * mg_cache%k2 * a**ss
            LKA2 = lambda2_2 * mg_cache%k2 * a**ss

            MGCAMB_Mudot = ((B1 - 1.d0) * mg_cache%adotoa * ss * LKA1) / ((1.d0+LKA1)**2.d0)

            if ( model==4 ) then ! correction for f(R) mu function.

                mu = MGCAMB_Mu( a, mg_par_cache, mg_cache )

                MGCAMB_Mudot = MGCAMB_Mudot/(1.d0 - 1.4d-8 * lambda1_2 * a**3) + 3.d0 * &
                            mu* mg_cache%adotoa *a**3 *(1.4d-8 * lambda1_2 ) &
                            /(1.d0 - 1.4d-8 * lambda1_2 * a**3)
            end if

        else if ( model ==6) then

            mu = MGCAMB_Mu( a, mg_par_cache, mg_cache )

            omm=(mg_par_cache%omegab+mg_par_cache%omegac)/((mg_par_cache%omegab+mg_par_cache%omegac) &
                & +(1-mg_par_cache%omegab-mg_par_cache%omegac)*a**3)
            ommdot=-3.d0*omm**2*a**3*mg_cache%adotoa*(1-mg_par_cache%omegab-mg_par_cache%omegac) &
                & /(mg_par_cache%omegab+mg_par_cache%omegac)

            MGCAMB_Mudot = mu/omm*(Linder_gamma-1.d0)*ommdot+&
                    2.d0/3.d0*omm**(Linder_gamma-1.d0)*ommdot*&
                    (Linder_gamma*omm**(Linder_gamma-1.d0)+3.d0*(Linder_gamma-0.5d0))

        else if (model==7 .or. &
                model==8 .or. &
                model==9 .or. &
                model==10) then

            beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
            m       = MGCAMB_M( a, mg_par_cache, mg_cache )
            betadot = MGCAMB_Betadot( a, mg_par_cache, mg_cache )
            mdot    = MGCAMB_Mdot( a, mg_par_cache, mg_cache )

            t1 = (2.d0*beta**2.d0)*k2
            t2 = (m**2.d0)*a**2.d0
            t1dot = 4.d0*beta*betadot*mg_cache%k2
            t2dot = (2.d0*a**2.d0)*(m*mdot+ (m**2.d0) *mg_cache%adotoa)

            MGCAMB_Mudot = (t1dot*(mg_cache%k2 + t2) - t1*t2dot)/((mg_cache%k2 + t2)**2.d0)

        end if

    end function MGCAMB_Mudot

    !-----------------------------------------------
    ! gamma(a,k) function
    function MGCAMB_Gamma( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Gamma                                    !< MG gamma function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        real(dl) :: LKA1 ! \lambda_1^2 k^2 a^s
        real(dl) :: LKA2 ! \lambda_1^2 k^2 a^s
        real(dl) :: t1,t2, t1dot, t2dot

        real(dl) :: beta, m

        if( model==1 .or. model==4 .or. model==5 ) then
            LKA1 = lambda1_2 * mg_cache%k2 * a**ss
            LKA2 = lambda2_2 * mg_cache%k2 * a**ss

            MGCAMB_Gamma = (1.d0 + B2 * LKA2)/(1.d0 +LKA2)

        else if ( model==6 ) then
                MGCAMB_Gamma = 1.d0

        else if ( model==7 .or. &
                model==8 .or. &
                model==9 .or. &
                model==10 ) then

            beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
            m       = MGCAMB_M( a, mg_par_cache, mg_cache )

            t1 = (2.d0*beta**2.d0)*mg_cache%k2
            t2 = (m**2.d0)*a**2.d0

            MGCAMB_Gamma = (mg_cache%k2 - t1 + t2)/(mg_cache%k2 + t1 + t2)

        end if

    end function MGCAMB_Gamma


    !-----------------------------------------------
    ! \dot{gamma}(a,k) function
    function MGCAMB_Gammadot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Gammadot                                 !< MG gammadot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters


        real(dl) :: LKA1 ! \lambda_1^2 k^2 a^s
        real(dl) :: LKA2 ! \lambda_1^2 k^2 a^s
        real(dl) :: k2
        real(dl) :: t1,t2,t1dot,t2dot

        real(dl) :: beta, betadot, m, mdot

        if(model==1 .or.model==4 .or.model==5) then
            LKA1 = lambda1_2 * mg_cache%k2 * a**ss
            LKA2 = lambda2_2 * mg_cache%k2 * a**ss

            MGCAMB_Gammadot = ((B2 -1.d0)*mg_cache%adotoa * ss* LKA2)/((1.d0+LKA2)**2.d0)

        else if ( model ==6) then
                MGCAMB_Gammadot = 0.d0

        else if (model == 7 .or. model==8 .or. model==9 .or. model==10) then

            beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
            m       = MGCAMB_M( a, mg_par_cache, mg_cache )
            betadot = MGCAMB_Betadot( a, mg_par_cache, mg_cache )
            mdot    = MGCAMB_Mdot( a, mg_par_cache, mg_cache )

            t1      = (2.d0*beta**2.d0)*mg_cache%k2
            t2      = (m**2.d0)*a**2.d0
            t1dot   = 4.d0*beta*betadot*mg_cache%k2
            t2dot   = (2.d0*a**2.d0)*(m*mdot + (m**2.d0) *mg_cache%adotoa)

            MGCAMB_Gammadot = 2.d0*(t1*t2dot-t1dot*(mg_cache%k2 + t2))/((mg_cache%k2 + t1 + t2)**2.d0)

        end if


    end function MGCAMB_Gammadot


!----------------------------------------------------------------------------------------------
!> MGCAMB (beta, m) parametrization, QSA for scalar-tensor models

    !-----------------------------------------------
    !> m(a) function
    function MGCAMB_M( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_M                                        !< MG m function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        real(dl) :: FRm0

        ! SYMMETRON
        if(model == 7) then
            MGCAMB_M = (mg_par_cache%H0/3.0D05) / (xi_star) * sqrt(1.d0-(a_star/a)**3.d0)

        ! DILATON: based on 1206.3568
        else if (model==8) then
            MGCAMB_M = (mg_par_cache%H0/3.0D05) /(xi0) * a**(- DilR)

        ! Hu-Sawicki f(R) model: m, beta parametrization as in 1305.5647
        else if (model == 9)then
            FRm0 = (mg_par_cache%h0/3.0D05)*sqrt((4.d0*mg_par_cache%omegav + mg_par_cache%omegab + mg_par_cache%omegac) &
                    & /((FRn+1.d0)*F_R0))!note factor of c here
            MGCAMB_M = FRm0 * ((4.d0 * mg_par_cache%omegav + (mg_par_cache%omegab + mg_par_cache%omegac)*a**(-3.d0)) &
                    & /(4.d0 * mg_par_cache%omegav + mg_par_cache%omegab + mg_par_cache%omegac))**(FRn/2.d0+1.d0)

        ! Simpler DILATON model
        else if ( model == 10 )then
            MGCAMB_M = sqrt(3.d0*A_2)*(mg_cache%adotoa/a)  ! H(a) = da/dtau/a**2 = adotoa/a

        end if

    end function MGCAMB_M


    !-----------------------------------------------
    !> \dot{m}(a) function
    function MGCAMB_Mdot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Mdot                                     !< MG mdot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        real(dl) :: FRm0
        real(dl) :: m

        m = MGCAMB_M( a, mg_par_cache, mg_cache )

        ! SYMMETRON
        if(model == 7) then
            MGCAMB_Mdot = 1.5d0*(mg_par_cache%H0/3.0D05)/(xi_star)*((a_star/a)**3.d0*mg_cache%adotoa)/(sqrt(1.d0-(a_star/a)**3.d0))

        ! DILATON
        else if (model==8) then
            MGCAMB_Mdot = - DilR * m * mg_cache%adotoa

        ! Hu-Sawicki f(R) model
        else if (model == 9)then

            FRm0 = (mg_par_cache%h0/3.0D05)*sqrt((4.d0*mg_par_cache%omegav + mg_par_cache%omegab + mg_par_cache%omegac)/ &
                    & ((FRn+1.d0)*F_R0))
            MGCAMB_Mdot = m / (4.d0 * mg_par_cache%omegav + (mg_par_cache%omegab + mg_par_cache%omegac)*a**(-3.d0)) &
                    & * (-3.d0*FRn/2.d0-3.d0)*((mg_par_cache%omegab + mg_par_cache%omegac)* a**(-3.d0)*mg_cache%adotoa)!/(4.d0 * mg_par_cache%omegav + mg_par_cache%omegab + mg_par_cache%omegac)) ! complete this

        ! Simple DILATON model
        else if (model ==10)then
            MGCAMB_Mdot = sqrt(3.d0*A_2)*(mg_cache%Hdot- mg_cache%adotoa**2.d0)/a !/3.0D05

        end if

    end function MGCAMB_Mdot

    !-----------------------------------------------
    !> beta(a) function
    function MGCAMB_Beta( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Beta                                     !< MG beta function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        ! SYMMETRON
        if(model == 7) then
            MGCAMB_Beta =  beta_star * sqrt(1.d0-(a_star/a)**3.d0)

        ! DILATON
        else if (model==8) then
            MGCAMB_Beta = beta0 * exp((DilS)/(2.d0* DilR - 3.d0)*(a**(2.d0* DilR - 3.d0)-1.d0))

        ! Hu-Sawicki f(R) model
        else if (model == 9)then
            MGCAMB_Beta = beta0

        ! Simple DILATON model
        else if (model ==10)then
            MGCAMB_Beta = beta0*(a**3.d0)
        end if

    end function MGCAMB_Beta

    !-----------------------------------------------
    !> \dot{beta}(a) function
    function MGCAMB_Betadot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Betadot                                  !< MG betadot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        real(dl) :: beta

        beta = MGCAMB_Beta( a, mg_par_cache, mg_cache )

        ! SYMMETRON
        if(model == 7) then
            MGCAMB_Betadot = 1.5d0 * (beta_star * (a_star/a)**3.d0 * mg_cache%adotoa) /( sqrt(1.d0-(a_star/a)**3.d0))

        ! DILATON
        else if (model==8) then
            MGCAMB_Betadot = beta * (DilS * a**(2.d0* DilR - 3.d0) *  mg_cache%adotoa)

        ! Hu-Sawicki f(R) model
        else if (model == 9)then
            MGCAMB_Betadot = 0.d0

        ! Simple DILATON model
        else if (model ==10)then
            MGCAMB_Betadot = 3.d0 *beta* mg_cache%adotoa

        end if

    end function MGCAMB_Betadot

    !----------------------------------------------------------------------------------------------
    !> MGCAMB (Q,R) parametrization, QSA for scalar-tensor models

    !-----------------------------------------------
    !> Q(a,k) function
    function MGCAMB_Q( a, mg_par_cache, mg_cache )
    implicit none
    real(dl) :: a                                               !< scale factor
    real(dl) :: MGCAMB_Q                                        !< MG Q function
    type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
    type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        if (model ==2) then
            MGCAMB_Q = MGQfix

        else if (model ==3) then
            MGCAMB_Q = 1.d0 + (Qnot - 1.d0)* a**sss

        end if

    end function MGCAMB_Q

    !-----------------------------------------------
    !> \dot{Q}(a,k) function
    function MGCAMB_Qdot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Qdot                                     !< MG Qdot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        if (model ==2) then
            MGCAMB_Qdot = 0.d0

        else if (model ==3) then
            MGCAMB_Qdot = (Qnot - 1.d0)*mg_cache%adotoa* sss* a**(sss)
        end if

    end function MGCAMB_Qdot

    !-----------------------------------------------
    ! R(a,k) function
    function MGCAMB_R( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_R                                        !< MG R function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters


        if (model ==2) then
            MGCAMB_R=MGRfix

        else if (model ==3) then
            MGCAMB_R = 1.d0 + (Rnot - 1.d0)* a**sss

        end if

    end function MGCAMB_R

    !-----------------------------------------------
    ! \dot{R}(a,k) function
    function MGCAMB_Rdot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Rdot                                     !< MG Rdot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        if (model ==2) then
            MGCAMB_Rdot = 0.d0

        else if (model ==3) then
            MGCAMB_Rdot = (Rnot - 1.d0)*mg_cache%adotoa* sss* a**(sss)

        end if

    end function MGCAMB_Rdot

    ! ---------------------------------------------------------------------------------------------
    !> Modifying the background
    subroutine MGCAMB_DarkEnergy( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache),  intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        real(dl) :: wnow

        if ( DE_model == 0 ) then
            mg_cache%grhov_t = 3.d0*mg_par_cache%h0_Mpc**2 * mg_par_cache%omegav *a**2
            mg_cache%gpresv_t = - mg_cache%grhov_t
        else if ( DE_model == 1 ) then
            mg_cache%grhov_t = 3.d0*mg_par_cache%h0_Mpc**2*mg_par_cache%omegav*a**(-1.d0-3.d0*wDE)
            mg_cache%gpresv_t = mg_cache%grhov_t * wDE
        else if (DE_model == 2 ) then
            wnow = w0+(1.d0-a)*wa
            mg_cache%grhov_t = 3.d0*mg_par_cache%h0_Mpc**2*mg_par_cache%omegav*a**(-1.d0-3.d0*wnow)
            mg_cache%gpresv_t = mg_cache%grhov_t * wnow
        else
            write(*,*) 'choose a DE model'
            stop
        end if

    end subroutine MGCAMB_DarkEnergy

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


end module MGCAMB

