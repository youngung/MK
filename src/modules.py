"""
Modules necessary for region B calculations
"""
import numpy as np
from scipy import optimize

def eps_dot_eq(sig,sig_eq,edot):
    """
    Arguments
    ---------
    sig     (stress tensor)
    sig_eq  (Equivalent stress, a scholar value)
    edot    (strain rate tensor
    """
    wdot = np.tensordot(sig,edot)
    return wdot/sig_eq

## Below is the two contributions: one from strain hardening; one from strain rate sensitivity.
from func_sr import *
from func_hard import *

def find_e11_dot(
        sig_b,func_yld,
        deps_b,eps_b_old,dt,func_F,func_G)
    """
    Main function of MK FLD
    Find e11_dot

    Given the current dsig_b

    Arguments
    ---------
    sig_b      (6D cauchy stress)
    func_yld   (characterized yield function)
    deps_b     (6D strain increment) given.
    eps_b_old  (6D strain at the previous step)
    dt         (time increment)
    func_F     (characterized strain rate function)
    func_G     (characterized strain hardening function)
    """
    sigb_eq      = func_yld(sig_b)
    e11dot_guess = 0. # guess

    ## objective function to minimize
    def objf(e11dot):
        deps_b[5] = e11dot_guess
        wrate = np.tensordot(sig_b,deps_b)
        eps_b_dot_eq_tilde = wrate / sigb_eq
        eps_b_tilde = eps_b_old + eps_b_dot_eq * dt
        return sigb_eq - func_F(eps_b_dot_eq_tilde) * func_G(eps_b_tilde)

    return optmize.minimize(objf,x0=0.)
