    module mghandles
    use CAMB
    use Precision
    use ModelParams
    use Transfer
    use iso_c_binding
    use MGCAMB
    use mGCAMB_cache
    implicit none

    contains

        subroutine MGCAMB_setparams(Params, flags, model_params, par_cache)
            type(CAMBparams)  :: Params
            Type(MGCAMBFlags) :: flags
            Type(MGCAMBModelParams) :: model_params
            Type(MGCAMB_parameter_cache) :: par_cache

            Params%mgcamb_flags = flags
            Params%mgcamb_model_pars = model_params
            Params%mgcamb_par_cache = par_cache

        end subroutine MGCAMB_setparams


    end module mghandles
