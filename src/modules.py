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
        deps_b_ref,eps_b_old,dt,func_F,func_G):
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
    sigb_eq = func_yld(sig_b)
    deps_b  = deps_b_ref.copy()

    ## objective function to minimize
    def objf(e11dot):
        deps_b[0] = e11dot
        deps_b[2] = -deps_b[0]-deps_b[1] ## incompressibility
        # wrate = np.dot(sig_b,deps_b) ## equivalent work rate
        wrate = 0.
        for i in xrange(3):
            wrate = wrate+sig_b[i]  *deps_b[i]
            wrate = wrate+sig_b[i+3]*deps_b[i+3]*2.

        eps_b_dot_eq_tilde = wrate / sigb_eq
        eps_b_tilde = eps_b_old + eps_b_dot_eq_tilde * dt
        F = func_F(eps_b_dot_eq_tilde)
        G = func_G(eps_b_tilde)
        # print '-'*20
        # print 'eps_b_dot_eq_tilde:', eps_b_dot_eq_tilde
        # print 'eps_b_tilde:', eps_b_tilde
        # print 'sigb_eq:', sigb_eq
        # print 'e11dot:', e11dot
        # print 'strain rate F:',F
        # print 'strain hard G:',G
        # print 'F*G:', F*G
        # print '|sigb_eq - F*G|:',abs(F*G-sigb_eq)
        # print '-'*20,'\n'
        return abs(sigb_eq - F*G)

    return objf
