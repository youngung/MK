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
from lib import assoc_flow_c

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

    """
    1.
    sig_b is the stress that B region should have
    at the current incremental step in order to meet
    the force equilibrium condition

    2.
    deps_b should be the incremental strain
    that B region should have at the current incremental
    step in order to meet the compatibility condition
    But only deps_b_11 is unknown and should be determined.
    """

    ## objective function to minimize
    def objf(e11dot):
        deps_b  = deps_b_ref[::]
        deps_b[0] = e11dot
        deps_b[2] = -deps_b[0]-deps_b[1] ## incompressibility


        ## Find correspondin "stress direction"
        ## Find sig_b_tilde that satisfies
        ## such that  e_tilde = assoc_flow_c(s_tilde_6,1.,func_yld)
        objf_sd = find_stress(func_yld,deps_b,sig_b)

        s22=sig_b[1]#10.
        irepeat = True; nit=0; nit_mx =5
        err_tol = 1e-1; dx=1e-4
        verbose = False
        while (irepeat):
            f = objf_sd(s22)
            if abs(f)<err_tol:
                irepeat=False
                break
            elif abs(f)>err_tol:
                f1 = objf_sd(s22-dx)
                f2 = objf_sd(s22+dx)
                jac_new = (f2-f1)/(dx*2.)
                if np.isinf(jac_new):
                    pass
                else:
                    jac = jac_new

                s22 = s22 - f / jac
                nit = nit+ 1
                print 'nit, f1, f2, jac, s22, f'
                print '%4i %4.1f %4.1f %8.1e %8.5f %8.5f'%(
                    nit, f1, f2, jac, s22, f)
                if nit>nit_mx:
                    raise IOError, 'Could not find the solution'

        raise IOError, 'debug'

        ## s22 was found.
        sig_b[1]=s22


        sigb_eq = func_yld(sig_b) ##

        # wrate = np.dot(sig_b,deps_b) ## equivalent work rate
        wrate = 0.
        for i in xrange(3):
            wrate = wrate+sig_b[i]  *deps_b[i]
            wrate = wrate+sig_b[i+3]*deps_b[i+3]*2.

        if wrate<0:
            print 'sigb'
            for j in xrange(6): print '%7.1f '%sig_b[j],
            print
            print 'deps_b'
            for j in xrange(6): print '%7.1e '%deps_b[j],
            print
            #raise IOError,'wrate should be positive'
            print 'wrate should be positive'

        if eps_b_old<0:
            raise IOError, 'epsb_old should be positive'

        eps_b_dot_eq_tilde = wrate / sigb_eq
        # if eps_b_dot_eq_tilde<0:
        #     eps_b_dot_eq_tilde=0.
        #     eps_b_dot_eq_tilde=1e-10 ## very low value

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


def find_stress(yfunc,d6,s6_old):
    """
    """
    d6_n = d6/np.sqrt(((d6**2).sum()))
    s6_n = np.sqrt((s6_old**2).sum())
    def objf(s22):
        s6 = s6_old[::]
        s6[1]=s22

        s6 = s6/np.sqrt((s6**2).sum())*s6_n

        sn=s6[::]
        # sn[:3] = sn[:3] - sn[2]

        if np.isnan(s22) or np.isinf(s22):
            raise IOError

        guess = assoc_flow_c(sn,1.,yfunc)
        guess = guess/np.sqrt((guess**2).sum())

        for i in xrange(6): print '%7.2f '%sn[i],
        print '|',
        for i in xrange(6): print '%7.2f '%guess[i],
        print '|',
        for i in xrange(6): print '%7.2f '%d6_n[i],
        print
        # print '\n\n'

        return np.sqrt(((d6_n-guess)**2).sum())

    return objf
