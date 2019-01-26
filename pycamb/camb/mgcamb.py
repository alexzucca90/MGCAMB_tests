# MGCAMB parameters

from .baseconfig import CAMB_Structure, CAMBError, camblib
from ctypes import c_int, c_double, c_bool, POINTER, byref
import numpy as np
from . import constants


class MGCAMBFlags(CAMB_Structure):
    """
        Object to store the MGCAMB flags.
    """

    _fields_ = [
                
                ("MG_flag", c_int),
                ("pure_MG_flag", c_int),
                ("alt_MG_flag", c_int),
                ("QSA_flag", c_int),
                ("mugamma_par", c_int),
                ("muSigma_par", c_int),
                ("QR_par", c_int),
                ("DE_model", c_int),
                ("GRtrans", c_double)
                ]

    def set_params(self, MG_flag=0, pure_MG_flag=1, alt_MG_flag=1, QSA_flag=1,
                   mugamma_par=2, muSigma_par=1, QR_par=1, DE_model=0, GRtrans=0.001):
        """
           Set the MGCAMB flags
        """

        self.MG_flag = MG_flag
        self.pure_MG_flag = pure_MG_flag
        self.alt_MG_flag = alt_MG_flag
        self.QSA_flag = QSA_flag
        self.mugamma_par = mugamma_par
        self.muSigma_par = muSigma_par
        self.QR_par = QR_par
        self.DE_model = DE_model
        self.GRtrans = GRtrans

        return self


class MGCAMBModelParams(CAMB_Structure):
    """
        Object to store the MGCAMB flags.
    """
    _fields_ = [
                ("B1", c_double),
                ("B2", c_double),
                ("lambda1_2", c_double),
                ("lambda2_2", c_double),
                ("ss", c_double),
                ("E11", c_double),
                ("E22", c_double),
                ("MGQfix", c_double),
                ("MGRfix", c_double),
                ("Qnot", c_double),
                ("Rnot", c_double),
                ("sss", c_double),
                ("Linder_gamma", c_double),
                ("beta_star", c_double),
                ("a_star", c_double),
                ("xi_star", c_double),
                ("xi0", c_double),
                ("DilR", c_double),
                ("DilS", c_double),
                ("F_R0", c_double),
                ("FRn", c_double),
                ("mu0", c_double),
                ("sigma0", c_double),
                ("w0DE", c_double),
                ("waDE", c_double)
                ]

    def set_params(self, B1=0.0, B2=0.0, lambda1_2=0.0, lambda2_2=0.0, ss=1.0, E11=0.0,
                   E22=0.0, MGQfix=0.0, MGRfix=0.0, Qnot=0.0, Rnot=0.0,sss=1.0, Linder_gamma=1.0,
                   beta_star=0.0, a_star=1.0e-2, xi_star=0.0,xi0=0.0, DilR=1.0, DilS=1.0,
                   F_R0=1.e-4, FRn=1.0, mu0=0.0, sigma0=0.0, w0DE=-1.0, waDE=0.0 ):
                   
        """
            Set the MGCAMB model parameters
        """

        self.B1 = B1
        self.B2 = B2
        self.lambda1_2 = lambda1_2
        self.lambda2_2 = lambda2_2
        self.ss = ss
        self.E11 = E11
        self.E22 = E22
        self.MGQfix = MGQfix
        self.MGRfix = MGRfix
        self.Qnot = Qnot
        self.Rnot = Rnot
        self.sss = sss
        self.Linder_gamma = Linder_gamma
        self.beta_star = beta_star
        self.a_star = a_star
        self.xi_star = xi_star
        self.xi0 = xi0
        self.DilR = DilR
        self.DilS = DilS
        self.F_R0 = F_R0
        self.FRn = FRn
        self.mu0 = mu0
        self.sigma0=sigma0
        self.w0DE = w0DE
        self.waDE = waDE

        return self

class MGCAMB_parameter_cache(CAMB_Structure):
    """
        Object to store the MGCAMB parameter cache
    """

    # the field output root is missing still

    _fields_ = [
                ("omegab", c_double),
                ("omegac", c_double),
                ("omegav", c_double),
                ("h0", c_double),
                ("h0_Mpc", c_double)
                ]

    def set_params(self, omegab=4.6e-2, omegac=2.28e-1, omegav=7.24e-1, h0=7.e1):
        """
            Set the MGCAMB parameter cache
        """

        self.omegab = omegab
        self.omegac = omegac
        self.omegav = omegav
        self.h0 = h0
        self.h0_Mpc = h0 * (1.e3/constants.c)

        return self


