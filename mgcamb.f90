module MGCAMB
    use precision
    use MGCAMB_cache
    use ModelParams


contains

    !---------------------------------------------------------------------------
    !> this subroutine computes the MG functions at a time-step
    subroutine MGCAMB_compute_MG_functions( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache), intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        ! Divide the cases here
        if (( CP%mgcamb_flags%MG_flag == 1 .and. CP%mgcamb_flags%pure_MG_flag /= 3 ) & ! all the mu, gamma parametrizations
            .or. CP%mgcamb_flags%MG_flag == 2 &
            .or. CP%mgcamb_flags%MG_flag == 3 ) then

            mg_cache%mu         = MGCAMB_Mu( a, mg_par_cache, mg_cache )
            mg_cache%mudot      = MGCAMB_MuDot( a, mg_par_cache, mg_cache )
            mg_cache%gamma      = MGCAMB_Gamma( a, mg_par_cache, mg_cache )
            mg_cache%gammadot   = MGCAMB_GammaDot( a, mg_par_cache, mg_cache )

            !write(*,*) 'a, k, mu, mudot, gamma, gammadot', a, mg_cache%k, mg_cache%mu,&
            !            mg_cache%mudot, mg_cache%gamma, mg_cache%gammadot

            ! other EFT functions are zero
            mg_cache%q      = 0._dl
            mg_cache%qdot   = 0._dl
            mg_cache%r      = 0._dl
            mg_cache%rdot   = 0._dl

        else if (  CP%mgcamb_flags%MG_flag == 1 .and. CP%mgcamb_flags%pure_MG_flag == 3  ) then ! the Q,R parametrization

            mg_cache%q      = MGCAMB_Q( a, mg_par_cache, mg_cache )
            mg_cache%qdot   = MGCAMB_Qdot( a, mg_par_cache, mg_cache )
            mg_cache%r      = MGCAMB_R( a, mg_par_cache, mg_cache )
            mg_cache%rdot   = MGCAMB_Rdot( a, mg_par_cache, mg_cache )

            ! other MG functions are zero
            mg_cache%mu         = 0._dl
            mg_cache%mudot      = 0._dl
            mg_cache%gamma      = 0._dl
            mg_cache%gammadot   = 0._dl

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

        if (( CP%mgcamb_flags%MG_flag == 1 .and. CP%mgcamb_flags%pure_MG_flag /= 3 ) & ! all the mu, gamma parametrizations
            .or. CP%mgcamb_flags%MG_flag == 2 &
            .or. CP%mgcamb_flags%MG_flag == 3 ) then

            ! first calculate MG_alpha
            mg_cache%MG_alpha = ( mg_cache%etak/mg_cache%k + mg_cache%mu * ( mg_cache%gamma*mg_cache%rhoDelta+ &
                                ( mg_cache%gamma -1._dl )*2._dl* mg_cache%dgpi)/(2._dl*mg_cache%k2)) / mg_cache%adotoa

            ! then calculate sigma
            mg_cache%sigma = mg_cache%k * mg_cache%MG_alpha

        else if (  CP%mgcamb_flags%MG_flag == 1 .and. CP%mgcamb_flags%pure_MG_flag == 3  ) then

            mg_cache%MG_phi      = - mg_cache%rhoDelta * mg_cache%q/(2._dl*mg_cache%k2)
            mg_cache%sigma       = (mg_cache%etak - mg_cache%k * mg_cache%MG_phi)/mg_cache%adotoa
            mg_cache%MG_alpha    = mg_cache%sigma/mg_cache%k

        end if

    end subroutine MGCAMB_compute_sigma

    !---------------------------------------------------------------------------
    !> this subroutine computes the perturbation Z in MG
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

        if (( CP%mgcamb_flags%MG_flag == 1 .and. CP%mgcamb_flags%pure_MG_flag /= 3 ) & ! all the mu, gamma parametrizations
            .or. CP%mgcamb_flags%MG_flag == 2 &
            .or. CP%mgcamb_flags%MG_flag == 3 ) then

            !> adding the massive neutrinos contibutions
            fmu = mg_cache%k2+0.5d0*mg_cache%gamma*mg_cache%mu*(3._dl*(mg_cache%grhoc_t+mg_cache%grhob_t) &
                & + 4._dl*(mg_cache%grhog_t+mg_cache%grhor_t) +3._dl * (mg_cache%grhonu_t + mg_cache%gpresnu_t ))

            !> adding massive neutrinos contributions

            f1 = mg_cache%k2+3._dl*( mg_cache%adotoa**2 - mg_cache%Hdot )
            !f1 = mg_cache%k2+0.5d0*(3._dl*(mg_cache%grhoc_t+mg_cache%grhob_t) &
            !    & + 4._dl*(mg_cache%grhog_t+mg_cache%grhor_t) + 3._dl*(mg_cache%grhonu_t+mg_cache%gpresnu_t) &
            !    & + 3._dl*(mg_cache%grhov_t+mg_cache%gpresv_t))

            term1 = mg_cache%gamma*mg_cache%mu* f1 * mg_cache%dgq/mg_cache%k

            !> adding massive neutrinos contribution, if w_DE /= -1 this has to be changed
            !term2 = mg_cache%k2*mg_cache%MG_alpha* ((mg_cache%mu* mg_cache%gamma- 1._dl)*(mg_cache%grhoc_t+mg_cache%grhob_t&
            !        & +(4._dl/3._dl)*(mg_cache%grhog_t+mg_cache%grhor_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t)) &
            !        & - (mg_cache%grhov_t+ mg_cache%gpresv_t))

            term2 = mg_cache%k2*mg_cache%MG_alpha* (mg_cache%mu* mg_cache%gamma*( mg_cache%grhoc_t+mg_cache%grhob_t   &
                    & +(4._dl/3._dl)*(mg_cache%grhog_t+mg_cache%grhor_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t) ) &
                    & - 2._dl*(mg_cache%adotoa**2 - mg_cache%Hdot))

            term3= (mg_cache%mu * ( mg_cache%gamma -1._dl)* mg_cache%adotoa - mg_cache%gamma*mg_cache%mudot &
                    & - mg_cache%gammadot*mg_cache%mu )*mg_cache%rhoDelta

            ! typo corrected here
            term4 = 2._dl*mg_cache%mu*(mg_cache%gamma - 1._dl)*mg_cache%adotoa*mg_cache%dgpi_w_sum

            ! separated fromt the previous term
            term5 = -2._dl*((mg_cache%gamma-1._dl)*mg_cache%mudot -mg_cache%gammadot*mg_cache%mu)*mg_cache%dgpi

            !> adding massive neutrinos contribution
            term6= 2._dl * mg_cache%mu*(1._dl - mg_cache%gamma)* mg_cache%pidot_sum

            !> calculate etadot
            mg_cache%etadot = (term1 + term2 + term3 + term4 + term5 + term6)/( 2._dl * fmu)

            !> finally calculate Z
            mg_cache%z = mg_cache%sigma - 3._dl * mg_cache%etadot/mg_cache%k

            !> Calculate the Newtonian potential
            mg_cache%MG_psi = - mg_cache%mu * ( mg_cache%rhoDelta + 2._dl* mg_cache%dgpi)/(2._dl*mg_cache%k2)

            !> calculate the curvature perturbation potential
            mg_cache%MG_phi = mg_cache%gamma * mg_cache%MG_psi + mg_cache%mu* 1._dl*mg_cache%dgpi/mg_cache%k2

            mg_cache%MG_phidot = mg_cache%etadot - mg_cache%adotoa * (mg_cache%MG_psi - mg_cache%adotoa * mg_cache%MG_alpha) &
                                & - mg_cache%Hdot * mg_cache%MG_alpha

        else if (  CP%mgcamb_flags%MG_flag == 1 .and. CP%mgcamb_flags%pure_MG_flag == 3  ) then

            ! adding massive neutrinos contributions
            fQ = mg_cache%k2 + 0.5d0*mg_cache%q * (3._dl*(mg_cache%grhob_t+mg_cache%grhoc_t)+&
                & 4._dl*(mg_cache%grhor_t+mg_cache%grhog_t)+3._dl*(mg_cache%grhonu_t + mg_cache%gpresnu_t))

            ! fixed for w_DE /= -1
            !f1=mg_cache%k2+3._dl*( mg_cache%adotoa**2 - mg_cache%Hdot )
            f1 = mg_cache%k2+0.5d0*(3._dl*(mg_cache%grhoc_t+mg_cache%grhob_t) &
                & + 4._dl*(mg_cache%grhog_t+mg_cache%grhor_t) + 3._dl*(mg_cache%grhonu_t+mg_cache%gpresnu_t) &
                & + 3._dl*(mg_cache%grhov_t+mg_cache%gpresv_t))

            k2alpha= mg_cache%k * mg_cache%sigma

            term1 = mg_cache%q * f1 * mg_cache%dgq/mg_cache%k

            term2 = k2alpha * ((mg_cache%q - 1._dl) * ( mg_cache%grhob_t+mg_cache%grhoc_t+(4._dl/3._dl) &
                    & *(mg_cache%grhor_t+mg_cache%grhog_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t) &
                    & ) -mg_cache%grhov_t - mg_cache%gpresv_t)

            !term2 = k2alpha * ((mg_cache%q) * ( mg_cache%grhob_t+mg_cache%grhoc_t+(4._dl/3._dl) &
            !    & *(mg_cache%grhor_t+mg_cache%grhog_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t)) &
            !    & - 2._dl *(mg_cache%adotoa**2 - mg_cache%Hdot))

            term3 = -( mg_cache%qdot + (mg_cache%r-1._dl) * mg_cache%q * mg_cache%adotoa ) * mg_cache%rhoDelta

            mg_cache%etadot = (term1 + term2 + term3)/( 2._dl * fQ )

            mg_cache%z = mg_cache%sigma - 3._dl * mg_cache%etadot/mg_cache%k

            !calculating also ISW related quantities
            mg_cache%MG_psi     = mg_cache%r * mg_cache%MG_phi - mg_cache%q * 1._dl * mg_cache%dgpi/mg_cache%k2
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

        term0 = mg_cache%k2 + 3._dl* (mg_cache%adotoa**2._dl - mg_cache%Hdot)

        !adding MG_rhoDeltadot
        mg_cache%rhoDeltadot = -term0 * mg_cache%dgq/mg_cache%k-(mg_cache%grho + mg_cache%gpres)* mg_cache%k*mg_cache%z&
                            &-mg_cache%adotoa*mg_cache%rhoDelta-2._dl * mg_cache%adotoa * mg_cache%dgpi

        !adding dgpidot
        mg_cache%dgpidot = mg_cache%pidot_sum - (2._dl*mg_cache%dgpi+ mg_cache%dgpi_diff )*mg_cache%adotoa

        if (( CP%mgcamb_flags%MG_flag == 1 .and. CP%mgcamb_flags%pure_MG_flag /= 3 ) & ! all the mu, gamma parametrizations
            .or. CP%mgcamb_flags%MG_flag == 2 &
            .or. CP%mgcamb_flags%MG_flag == 3 ) then

            mg_cache%MG_psidot = - 0.5d0*mg_cache%mu/mg_cache%k2*(mg_cache%rhoDeltadot+2._dl*mg_cache%dgpidot) &
                                & - 0.5d0*mg_cache%mudot/mg_cache%k2*(mg_cache%rhoDelta+2._dl*mg_cache%dgpi)

        else if (  CP%mgcamb_flags%MG_flag == 1 .and. CP%mgcamb_flags%pure_MG_flag == 3  ) then

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

        real(dl) :: omegaDE_t

        ! beta, m parametrization
        real(dl) :: beta, m

        !> pure MG models
        if ( CP%mgcamb_flags%MG_flag == 1 .and. CP%mgcamb_flags%pure_MG_flag /= 3 ) then

            if ( CP%mgcamb_flags%pure_MG_flag == 1 ) then ! mu-gamma
                if ( CP%mgcamb_flags%mugamma_par == 1 ) then ! BZ parametrization
                    LKA1 = CP%mgcamb_model_pars%lambda1_2 * mg_cache%k2 * a**CP%mgcamb_model_pars%ss
                    LKA2 = CP%mgcamb_model_pars%lambda2_2 * mg_cache%k2 * a**CP%mgcamb_model_pars%ss

                    MGCAMB_Mu = (1._dl + CP%mgcamb_model_pars%B1 * LKA1)/(1._dl + LKA1)

                else if ( CP%mgcamb_flags%mugamma_par == 2 ) then ! Planck parametrization

                    ! changing the following
                    !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2

                    omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                    MGCAMB_Mu = 1._dl + CP%mgcamb_model_pars%E11*omegaDE_t

                else if ( CP%mgcamb_flags%mugamma_par == 3 ) then
                    MGCAMB_Mu = 1._dl

                end if

            else if ( CP%mgcamb_flags%pure_MG_flag == 2 ) then ! mu-Sigma

                if ( CP%mgcamb_flags%muSigma_par == 1 ) then ! DES parametrization

                    !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                    !MGCAMB_Mu = 1._dl + mu0 * omegaDE_t/mg_par_cache%omegav

                    ! this is being changed
                    omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                    MGCAMB_Mu = 1._dl + CP%mgcamb_model_pars%mu0 * omegaDE_t/mg_par_cache%omegav


                else if ( CP%mgcamb_flags%muSigma_par == 2 ) then
                    MGCAMB_Mu = 1._dl

                end if

            end if

        !> alternative MG
        else if ( CP%mgcamb_flags%MG_flag == 2 ) then

            if (CP%mgcamb_flags%alt_MG_flag == 1) then !(Linder Gamma)
                omm=(mg_par_cache%omegab+mg_par_cache%omegac)/((mg_par_cache%omegab+mg_par_cache%omegac) &
                    & + (1-mg_par_cache%omegab-mg_par_cache%omegac)*a**3)
                ommdot=-3._dl*omm**2*a**3*mg_cache%adotoa*(1-mg_par_cache%omegab-mg_par_cache%omegac) &
                    & /(mg_par_cache%omegab+mg_par_cache%omegac)

                MGCAMB_Mu=2._dl/3._dl*omm**(CP%mgcamb_model_pars%Linder_gamma-1._dl)*&
                    & (omm**CP%mgcamb_model_pars%Linder_gamma+2-3._dl*CP%mgcamb_model_pars%Linder_gamma&
                    & +3._dl*(CP%mgcamb_model_pars%Linder_gamma-0.5d0)*omm)

            else if ( CP%mgcamb_flags%alt_MG_flag == 2 ) then
                MGCAMB_Mu = 1._dl
            end if


        !> QSA models
        else if ( CP%mgcamb_flags%MG_flag == 3 ) then

            if ( CP%mgcamb_flags%QSA_flag == 1 ) then ! f(R)
                LKA1 = CP%mgcamb_model_pars%lambda1_2 * mg_cache%k2 * a**CP%mgcamb_model_pars%ss
                LKA2 = CP%mgcamb_model_pars%lambda2_2 * mg_cache%k2 * a**CP%mgcamb_model_pars%ss
                MGCAMB_Mu = (1._dl + CP%mgcamb_model_pars%B1 * LKA1)/(1._dl + LKA1)
                MGCAMB_Mu = MGCAMB_Mu/(1._dl - 1.4d-8 * CP%mgcamb_model_pars%lambda1_2 * a**3)

            else if ( CP%mgcamb_flags%QSA_flag == 2 .or. &  ! beta, m parametrization
                      CP%mgcamb_flags%QSA_flag == 3 .or. &
                      CP%mgcamb_flags%QSA_flag == 4 ) then
                beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                m       = MGCAMB_M( a, mg_par_cache, mg_cache )
                t1      = (2._dl*beta**2._dl)*mg_cache%k2
                t2      = (m**2._dl)*a**2._dl

                MGCAMB_Mu = (mg_cache%k2 + t1 + t2)/(mg_cache%k2 + t2)


            else if ( CP%mgcamb_flags%QSA_flag == 5 )  then
                MGCAMB_Mu = 1._dl

            end if

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

        real(dl) :: omegaDEdot

        !> pure MG models
        if ( CP%mgcamb_flags%MG_flag == 1 .and. CP%mgcamb_flags%pure_MG_flag /= 3 ) then

            if ( CP%mgcamb_flags%pure_MG_flag == 1 ) then ! mu-gamma
                if ( CP%mgcamb_flags%mugamma_par == 1 ) then ! BZ parametrization
                    LKA1 = CP%mgcamb_model_pars%lambda1_2 * mg_cache%k2 * a**CP%mgcamb_model_pars%ss
                    LKA2 = CP%mgcamb_model_pars%lambda2_2 * mg_cache%k2 * a**CP%mgcamb_model_pars%ss

                    MGCAMB_Mudot = ((CP%mgcamb_model_pars%B1 - 1._dl) * mg_cache%adotoa * CP%mgcamb_model_pars%ss * LKA1)&
                        & / ((1._dl+LKA1)**2._dl)

                else if ( CP%mgcamb_flags%mugamma_par == 2 ) then ! Planck parametrization

                    ! changingh the following quantity
                    !omegaDEdot = - 3._dl * mg_cache%adotoa * (mg_cache%grhov_t + mg_cache%gpresv_t) &
                    !            & / a**2 / 3._dl / mg_par_cache%h0_Mpc**2

                    omegaDEdot=-(mg_cache%grhov_t+3._dl*mg_cache%gpresv_t)/3._dl/mg_cache%adotoa &
                            & - 2._dl*mg_cache%Hdot/3._dl/mg_cache%adotoa**3*mg_cache%grhov_t

                    MGCAMB_Mudot = CP%mgcamb_model_pars%E11*omegaDEdot

                else if ( CP%mgcamb_flags%mugamma_par == 3 ) then
                    MGCAMB_Mudot = 0._dl

                end if

            else if ( CP%mgcamb_flags%pure_MG_flag == 2 ) then ! mu-Sigma

                if ( CP%mgcamb_flags%muSigma_par == 1 ) then ! DES parametrization
                    ! changing the following
                    !omegaDEdot = - 3._dl * mg_cache%adotoa * (mg_cache%grhov_t + mg_cache%gpresv_t) &
                    !            & / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                    omegaDEdot=-(mg_cache%grhov_t+3._dl*mg_cache%gpresv_t)/3._dl/mg_cache%adotoa &
                                & - 2._dl*mg_cache%Hdot/3._dl/mg_cache%adotoa**3*mg_cache%grhov_t

                    MGCAMB_Mudot =  CP%mgcamb_model_pars%mu0 * omegaDEdot/mg_par_cache%omegav

                else if ( CP%mgcamb_flags%muSigma_par == 2 ) then
                    MGCAMB_Mudot = 0._dl

                end if

            end if

        !> alternative MG
        else if ( CP%mgcamb_flags%MG_flag == 2 ) then

            if ( CP%mgcamb_flags%alt_MG_flag == 1 ) then !(Linder Gamma)
                mu = MGCAMB_Mu( a, mg_par_cache, mg_cache )

                omm=(mg_par_cache%omegab+mg_par_cache%omegac)/((mg_par_cache%omegab+mg_par_cache%omegac) &
                    & +(1-mg_par_cache%omegab-mg_par_cache%omegac)*a**3)
                ommdot=-3._dl*omm**2*a**3*mg_cache%adotoa*(1-mg_par_cache%omegab-mg_par_cache%omegac) &
                    & /(mg_par_cache%omegab+mg_par_cache%omegac)

                MGCAMB_Mudot = mu/omm*(CP%mgcamb_model_pars%Linder_gamma-1._dl)*ommdot+&
                    &2._dl/3._dl*omm**(CP%mgcamb_model_pars%Linder_gamma-1._dl)*ommdot*&
                    &(CP%mgcamb_model_pars%Linder_gamma*omm**(CP%mgcamb_model_pars%Linder_gamma-1._dl)&
                    &+3._dl*(CP%mgcamb_model_pars%Linder_gamma-0.5d0))

            else if ( CP%mgcamb_flags%alt_MG_flag == 2 ) then
                MGCAMB_Mudot = 0._dl
            end if


        !> QSA models
        else if ( CP%mgcamb_flags%MG_flag == 3 ) then

            if ( CP%mgcamb_flags%QSA_flag == 1 ) then ! f(R)
                LKA1 = CP%mgcamb_model_pars%lambda1_2 * mg_cache%k2 * a**CP%mgcamb_model_pars%ss
                LKA2 = CP%mgcamb_model_pars%lambda2_2 * mg_cache%k2 * a**CP%mgcamb_model_pars%ss

                MGCAMB_Mudot = ((CP%mgcamb_model_pars%B1 - 1._dl) * mg_cache%adotoa * CP%mgcamb_model_pars%ss&
                    & * LKA1) / ((1._dl+LKA1)**2._dl)
                mu = MGCAMB_Mu( a, mg_par_cache, mg_cache )
                MGCAMB_Mudot = MGCAMB_Mudot/(1._dl - 1.4d-8 * CP%mgcamb_model_pars%lambda1_2 * a**3) + 3._dl * &
                                mu* mg_cache%adotoa *a**3 *(1.4d-8 * CP%mgcamb_model_pars%lambda1_2 ) &
                                /(1._dl - 1.4d-8 * CP%mgcamb_model_pars%lambda1_2 * a**3)

            else if ( CP%mgcamb_flags%QSA_flag == 2 .or. &  ! beta, m parametrization
                    CP%mgcamb_flags%QSA_flag == 3 .or. &
                    CP%mgcamb_flags%QSA_flag == 4 ) then

                beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                m       = MGCAMB_M( a, mg_par_cache, mg_cache )
                betadot = MGCAMB_Betadot( a, mg_par_cache, mg_cache )
                mdot    = MGCAMB_Mdot( a, mg_par_cache, mg_cache )

                t1 = (2._dl*beta**2._dl)*k2
                t2 = (m**2._dl)*a**2._dl
                t1dot = 4._dl*beta*betadot*mg_cache%k2
                t2dot = (2._dl*a**2._dl)*(m*mdot+ (m**2._dl) *mg_cache%adotoa)

                MGCAMB_Mudot = (t1dot*(mg_cache%k2 + t2) - t1*t2dot)/((mg_cache%k2 + t2)**2._dl)


            else if ( CP%mgcamb_flags%QSA_flag == 5 )  then
                MGCAMB_Mudot = 0._dl

            end if

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
        real(dl) :: omegaDE_t

        real(dl) :: sigma_t
        real(dl) :: mu_t

        !> pure MG models
        if ( CP%mgcamb_flags%MG_flag == 1 .and. CP%mgcamb_flags%pure_MG_flag /= 3 ) then

            if ( CP%mgcamb_flags%pure_MG_flag == 1 ) then ! mu-gamma
                if ( CP%mgcamb_flags%mugamma_par == 1 ) then ! BZ parametrization
                    LKA1 = CP%mgcamb_model_pars%lambda1_2 * mg_cache%k2 * a**CP%mgcamb_model_pars%ss
                    LKA2 = CP%mgcamb_model_pars%lambda2_2 * mg_cache%k2 * a**CP%mgcamb_model_pars%ss

                    MGCAMB_Gamma = (1._dl + CP%mgcamb_model_pars%B2 * LKA2)/(1._dl +LKA2)

                else if ( CP%mgcamb_flags%mugamma_par == 2 ) then ! Planck parametrization
                    ! changing the following
                    !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                    omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                    MGCAMB_Gamma = 1._dl+CP%mgcamb_model_pars%E22*omegaDE_t

                else if ( CP%mgcamb_flags%mugamma_par == 3 ) then
                    MGCAMB_Gamma = 1._dl

                end if

            else if ( CP%mgcamb_flags%pure_MG_flag == 2 ) then ! mu-Sigma

                if ( CP%mgcamb_flags%muSigma_par == 1 ) then ! DES parametrization
                    ! changing the following
                    !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                    omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                    sigma_t = 1._dl + CP%mgcamb_model_pars%sigma0 * omegaDE_t / mg_par_cache%omegav
                    mu_t    = 1._dl + CP%mgcamb_model_pars%mu0 * omegaDE_t / mg_par_cache%omegav
                    MGCAMB_Gamma = 2._dl * sigma_t / mu_t - 1._dl

                else if ( CP%mgcamb_flags%muSigma_par == 2 ) then
                    MGCAMB_Gamma = 1._dl

                end if

            end if

        !> alternative MG
        else if ( CP%mgcamb_flags%MG_flag == 2 ) then

            if ( CP%mgcamb_flags%alt_MG_flag == 1 ) then !(Linder Gamma)
                MGCAMB_Gamma = 1._dl

            else if ( CP%mgcamb_flags%alt_MG_flag == 2 ) then
                MGCAMB_Gamma = 1._dl
            end if


        !> QSA models
        else if ( CP%mgcamb_flags%MG_flag == 3 ) then

            if ( CP%mgcamb_flags%QSA_flag == 1 ) then ! f(R)
                LKA1 = CP%mgcamb_model_pars%lambda1_2 * mg_cache%k2 * a**CP%mgcamb_model_pars%ss
                LKA2 = CP%mgcamb_model_pars%lambda2_2 * mg_cache%k2 * a**CP%mgcamb_model_pars%ss

                MGCAMB_Gamma = (1._dl + CP%mgcamb_model_pars%B2 * LKA2)/(1._dl +LKA2)

            else if ( CP%mgcamb_flags%QSA_flag == 2 .or. &  ! beta, m parametrization
                    CP%mgcamb_flags%QSA_flag == 3 .or. &
                    CP%mgcamb_flags%QSA_flag == 4 ) then

                beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                m       = MGCAMB_M( a, mg_par_cache, mg_cache )

                t1 = (2._dl*beta**2._dl)*mg_cache%k2
                t2 = (m**2._dl)*a**2._dl

                MGCAMB_Gamma = (mg_cache%k2 - t1 + t2)/(mg_cache%k2 + t1 + t2)


            else if ( CP%mgcamb_flags%QSA_flag == 5 )  then
                MGCAMB_Gamma = 1._dl

            end if

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
        real(dl) :: omegaDE_t, omegaDEdot
        real(dl) :: sigma_t, sigmadot_t
        real(dl) :: mu_t, mudot_t



        !> pure MG models
        if ( CP%mgcamb_flags%MG_flag == 1 .and. CP%mgcamb_flags%pure_MG_flag /= 3 ) then

            if ( CP%mgcamb_flags%pure_MG_flag == 1 ) then ! mu-gamma
                if ( CP%mgcamb_flags%mugamma_par == 1 ) then ! BZ parametrization
                    LKA1 = CP%mgcamb_model_pars%lambda1_2 * mg_cache%k2 * a**CP%mgcamb_model_pars%ss
                    LKA2 = CP%mgcamb_model_pars%lambda2_2 * mg_cache%k2 * a**CP%mgcamb_model_pars%ss

                    MGCAMB_Gammadot = ((CP%mgcamb_model_pars%B2-1._dl)*mg_cache%adotoa*CP%mgcamb_model_pars%ss*&
                        &LKA2)/((1._dl+LKA2)**2._dl)

                else if ( CP%mgcamb_flags%mugamma_par == 2 ) then ! Planck parametrization
                    ! changing the following
                    !omegaDEdot = - 3._dl * mg_cache%adotoa * (mg_cache%grhov_t + mg_cache%gpresv_t) &
                    !            & / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                    omegaDEdot=-(mg_cache%grhov_t+3._dl*mg_cache%gpresv_t)/3._dl/mg_cache%adotoa &
                                & - 2._dl*mg_cache%Hdot/3._dl/mg_cache%adotoa**3*mg_cache%grhov_t

                    MGCAMB_Gammadot = CP%mgcamb_model_pars%E22*omegaDEdot

                else if ( CP%mgcamb_flags%mugamma_par == 3 ) then
                    MGCAMB_Gammadot = 0._dl

                end if

            else if ( CP%mgcamb_flags%pure_MG_flag == 2 ) then ! mu-Sigma

                if ( CP%mgcamb_flags%muSigma_par == 1 ) then ! DES parametrization

                ! changing the following
                !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                !omegaDEdot = - 3._dl * mg_cache%adotoa * (mg_cache%grhov_t + mg_cache%gpresv_t) &
                !            & / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                omegaDEdot  =-(mg_cache%grhov_t+3._dl*mg_cache%gpresv_t)/3._dl/mg_cache%adotoa &
                                & - 2._dl*mg_cache%Hdot/3._dl/mg_cache%adotoa**3*mg_cache%grhov_t
                sigma_t     = 1._dl + CP%mgcamb_model_pars%sigma0 * omegaDE_t / mg_par_cache%omegav
                sigmadot_t  = CP%mgcamb_model_pars%sigma0 * omegaDEdot / mg_par_cache%omegav
                mu_t        = 1._dl + CP%mgcamb_model_pars%mu0 * omegaDE_t / mg_par_cache%omegav
                mudot_t     = CP%mgcamb_model_pars%mu0 * omegaDEdot / mg_par_cache%omegav
                MGCAMB_Gammadot = 2._dl * sigmadot_t / mu_t - 2._dl *sigma_t*mudot_t/mu_t**2

                else if ( CP%mgcamb_flags%muSigma_par == 2 ) then
                    MGCAMB_Gammadot = 0._dl

            end if

            end if

        !> alternative MG
        else if ( CP%mgcamb_flags%MG_flag == 2 ) then

            if ( CP%mgcamb_flags%alt_MG_flag == 1 ) then !(Linder Gamma)
                MGCAMB_Gammadot = 0._dl

            else if ( CP%mgcamb_flags%alt_MG_flag == 2 ) then
                MGCAMB_Gammadot = 0._dl
            end if


        !> QSA models
        else if ( CP%mgcamb_flags%MG_flag == 3 ) then

            if ( CP%mgcamb_flags%QSA_flag == 1 ) then ! f(R)
                LKA1 = CP%mgcamb_model_pars%lambda1_2 * mg_cache%k2 * a**CP%mgcamb_model_pars%ss
                LKA2 = CP%mgcamb_model_pars%lambda2_2 * mg_cache%k2 * a**CP%mgcamb_model_pars%ss

                MGCAMB_Gammadot =((CP%mgcamb_model_pars%B2-1._dl)*mg_cache%adotoa*CP%mgcamb_model_pars%ss*LKA2)/&
                    &((1._dl+LKA2)**2._dl)

            else if ( CP%mgcamb_flags%QSA_flag == 2 .or. &  ! beta, m parametrization
                    CP%mgcamb_flags%QSA_flag == 3  .or. &
                    CP%mgcamb_flags%QSA_flag == 4 ) then

                beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                m       = MGCAMB_M( a, mg_par_cache, mg_cache )
                betadot = MGCAMB_Betadot( a, mg_par_cache, mg_cache )
                mdot    = MGCAMB_Mdot( a, mg_par_cache, mg_cache )

                t1      = (2._dl*beta**2._dl)*mg_cache%k2
                t2      = (m**2._dl)*a**2._dl
                t1dot   = 4._dl*beta*betadot*mg_cache%k2
                t2dot   = (2._dl*a**2._dl)*(m*mdot + (m**2._dl) *mg_cache%adotoa)

                MGCAMB_Gammadot = 2._dl*(t1*t2dot-t1dot*(mg_cache%k2 + t2))/((mg_cache%k2 + t1 + t2)**2._dl)

            else if ( CP%mgcamb_flags%QSA_flag == 5 )  then
                MGCAMB_Gammadot = 0._dl

            end if

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
        if( CP%mgcamb_flags%QSA_flag ==  2 ) then
            MGCAMB_M = (mg_par_cache%H0/3.0D05) / (CP%mgcamb_model_pars%xi_star)*sqrt(1._dl-(CP%mgcamb_model_pars%a_star/a)**3._dl)

        ! DILATON: based on 1206.3568
        else if ( CP%mgcamb_flags%QSA_flag ==  3 ) then
            MGCAMB_M = (mg_par_cache%H0/3.0D05) /(CP%mgcamb_model_pars%xi0) * a**(- CP%mgcamb_model_pars%DilR)

        ! Hu-Sawicki f(R) model: m, beta parametrization as in 1305.5647
        else if ( CP%mgcamb_flags%QSA_flag ==  4 )then
            FRm0 = (mg_par_cache%h0/3.0D05)*sqrt((4._dl*mg_par_cache%omegav + mg_par_cache%omegab + mg_par_cache%omegac) &
                    & /((CP%mgcamb_model_pars%FRn+1._dl)*CP%mgcamb_model_pars%F_R0))!note factor of c here
            MGCAMB_M = FRm0 * ((4._dl * mg_par_cache%omegav + (mg_par_cache%omegab + mg_par_cache%omegac)*a**(-3._dl)) &
                    & /(4._dl * mg_par_cache%omegav + mg_par_cache%omegab + mg_par_cache%omegac))**&
                    &(CP%mgcamb_model_pars%FRn/2._dl+1._dl)

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
        if( CP%mgcamb_flags%QSA_flag ==  2 ) then
            MGCAMB_Mdot = 1.5d0*(mg_par_cache%H0/3.0D05)/(CP%mgcamb_model_pars%xi_star)*&
                    &((CP%mgcamb_model_pars%a_star/a)**3._dl*mg_cache%adotoa)/&
                    & (sqrt(1._dl-(CP%mgcamb_model_pars%a_star/a)**3._dl))

        ! DILATON
        else if ( CP%mgcamb_flags%QSA_flag ==  3 ) then
            MGCAMB_Mdot = - CP%mgcamb_model_pars%DilR * m * mg_cache%adotoa

        ! Hu-Sawicki f(R) model
        else if ( CP%mgcamb_flags%QSA_flag ==  4 )then

            FRm0 = (mg_par_cache%h0/3.0D05)*sqrt((4._dl*mg_par_cache%omegav+mg_par_cache%omegab+mg_par_cache%omegac)/ &
                    & ((CP%mgcamb_model_pars%FRn+1._dl)*CP%mgcamb_model_pars%F_R0))
            MGCAMB_Mdot = m / (4._dl * mg_par_cache%omegav + (mg_par_cache%omegab + mg_par_cache%omegac)*a**(-3._dl)) &
                    & * (-3._dl*CP%mgcamb_model_pars%FRn/2._dl-3._dl)*((mg_par_cache%omegab + mg_par_cache%omegac)&
                    &* a**(-3._dl)*mg_cache%adotoa)!/(4._dl * mg_par_cache%omegav + mg_par_cache%omegab + mg_par_cache%omegac)) ! complete this

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
        if( CP%mgcamb_flags%QSA_flag == 2 ) then
            MGCAMB_Beta =  CP%mgcamb_model_pars%beta_star * sqrt(1._dl-(CP%mgcamb_model_pars%a_star/a)**3._dl)

        ! DILATON
        else if ( CP%mgcamb_flags%QSA_flag == 3 ) then
            MGCAMB_Beta = CP%mgcamb_model_pars%beta0 * exp((CP%mgcamb_model_pars%DilS)/(2._dl* &
                &CP%mgcamb_model_pars%DilR - 3._dl)*(a**(2._dl* CP%mgcamb_model_pars%DilR - 3._dl)-1._dl))

        ! Hu-Sawicki f(R) model
        else if ( CP%mgcamb_flags%QSA_flag == 4 )then
            MGCAMB_Beta = CP%mgcamb_model_pars%beta0

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
        if( CP%mgcamb_flags%QSA_flag == 2 ) then
            MGCAMB_Betadot = 1.5d0 * (CP%mgcamb_model_pars%beta_star * (CP%mgcamb_model_pars%a_star/a)**3._dl &
                &* mg_cache%adotoa) /( sqrt(1._dl-(CP%mgcamb_model_pars%a_star/a)**3._dl))

        ! DILATON
        else if ( CP%mgcamb_flags%QSA_flag == 3 ) then
            MGCAMB_Betadot = beta * (CP%mgcamb_model_pars%DilS * a**(2._dl* CP%mgcamb_model_pars%DilR - 3._dl) &
                &*  mg_cache%adotoa)

        ! Hu-Sawicki f(R) model
        else if ( CP%mgcamb_flags%QSA_flag == 4 )then
            MGCAMB_Betadot = 0._dl


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

        if ( CP%mgcamb_flags%QR_par == 1) then
            MGCAMB_Q = CP%mgcamb_model_pars%MGQfix

        else if ( CP%mgcamb_flags%QR_par == 2 ) then
            MGCAMB_Q = 1._dl + (CP%mgcamb_model_pars%Qnot - 1._dl)* a**CP%mgcamb_model_pars%sss

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

        if ( CP%mgcamb_flags%QR_par == 1 ) then
            MGCAMB_Qdot = 0._dl

        else if ( CP%mgcamb_flags%QR_par == 2 ) then
            MGCAMB_Qdot = (CP%mgcamb_model_pars%Qnot - 1._dl)*mg_cache%adotoa* CP%mgcamb_model_pars%sss* &
                &a**(CP%mgcamb_model_pars%sss)
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


        if ( CP%mgcamb_flags%QR_par == 1 ) then
            MGCAMB_R=CP%mgcamb_model_pars%MGRfix

        else if ( CP%mgcamb_flags%QR_par == 2 ) then
            MGCAMB_R = 1._dl + (CP%mgcamb_model_pars%Rnot - 1._dl)* a**CP%mgcamb_model_pars%sss

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

        if ( CP%mgcamb_flags%QR_par == 1 ) then
            MGCAMB_Rdot = 0._dl

        else if ( CP%mgcamb_flags%QR_par == 2 ) then
            MGCAMB_Rdot = (CP%mgcamb_model_pars%Rnot - 1._dl)*mg_cache%adotoa* CP%mgcamb_model_pars%sss*&
                & a**(CP%mgcamb_model_pars%sss)

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

        if ( CP%mgcamb_flags%DE_model == 0 ) then
            mg_cache%grhov_t = 3._dl*mg_par_cache%h0_Mpc**2 * mg_par_cache%omegav *a**2
            mg_cache%gpresv_t = - mg_cache%grhov_t
        else if ( CP%mgcamb_flags%DE_model == 1 ) then
            mg_cache%grhov_t = 3._dl*mg_par_cache%h0_Mpc**2*mg_par_cache%omegav*&
                &a**(-1._dl-3._dl*CP%mgcamb_model_pars%w0DE)
            mg_cache%gpresv_t = mg_cache%grhov_t * CP%mgcamb_model_pars%w0DE
        else if ( CP%mgcamb_flags%DE_model == 2 ) then
            wnow = CP%mgcamb_model_pars%w0DE+(1._dl-a)*CP%mgcamb_model_pars%waDE
            mg_cache%grhov_t = 3._dl*mg_par_cache%h0_Mpc**2*mg_par_cache%omegav*a**(-1._dl-3._dl*wnow)
            mg_cache%gpresv_t = mg_cache%grhov_t * wnow
        else if ( CP%mgcamb_flags%DE_model == 3 ) then
            write(*,*) 'This will contain the reconstruction of w_DE(a)'
            write(*,*) 'Not implemented yet'
            stop
        else if ( CP%mgcamb_flags%DE_model == 4 ) then
            write(*,*) 'This will contain the reconstruction of rho_DE(a)'
            write(*,*) 'Not implemented yet'
            stop
        else
            write(*,*) 'choose a DE model'
            stop
        end if

    end subroutine MGCAMB_DarkEnergy


end module MGCAMB

