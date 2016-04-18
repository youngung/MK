import time
import lib;reload(lib)
from lib import *
import yf;reload(yf)
from yf import *
import func_hard; reload(func_hard)
from func_hard import *
from func_sr import *
from MP import progress_bar
uet = progress_bar.update_elapsed_time
import multiprocessing
from multiprocessing import Pool
from scipy import optimize
import numpy as np


## Material characteristics
func_hd   = c_G(0,k=6.00e2, eps_0 = 4.23e-4,n=2.55e-1)
func_sr   = c_F(0,ed0=1e-3,m=1.2)
func_yd  = return_c_yld(iopt=1) ## Von Mises

## maximum strain....
eps  = np.linspace(0,1.2,100)
sigma_bar = lambda eps: func_hd(eps)

def test_FLDA_onepath():
    fig=plt.figure(figsize=(7,6))
    ax1=fig.add_subplot(221)
    ax2=fig.add_subplot(222)
    ax3=fig.add_subplot(223)
    # ax4=fig.add_subplot(224)

    ## Material Characteristics
    func_hd   = c_G(0,k=6.00e2, eps_0 = 4.23e-4,n=2.55e-1)
    func_sr   = c_F(0,ed0=1e-3,m=1.2)
    func_yd  = return_c_yld(iopt=0,R0=1.5,R90=0.7) ## QuadHill
    # func_yd  = return_c_yld(iopt=1) ## Von Mises

    alphs = np.linspace(0,1,5)
    for i in xrange(len(alphs)):
        strain_6, stress_6, strain_eq, \
            stress_eq, strain_rate_6, \
            time_stamps,delta_time = FLDA_onepath(
                alpha=alphs[i],beta=0,sr_eq=1e-3,
                debar=1e-4,
                ebar_mx=0.1,
                func_yd=func_yd,
                func_hd=func_hd,
                func_sr=func_sr)


        ax1.plot(strain_6[:,1],strain_6[:,0],'r-')
        ax2.plot(stress_6[:,1],stress_6[:,0],'b-')
        ax3.plot(time_stamps,strain_rate_6[:,2],'g-')

    fn = 'FLDA_onepath_test.pdf'
    fig.savefig(fn)
    print '%s has been saved'%fn

def FLDA_onepath(
        alpha,beta,sr_eq,
        debar,ebar_mx,
        func_yd,func_hd,func_sr):

    """
    Arguments
    ---------
    alpha, beta (alpha, beta for stress/strain path boundary condition)
    sr_eq       (strain rate)
    func_yd    (yield function)
    func_hd     (strain hardening function)
    func_sr     (strain rate function)

    ebar_mx     (maximum ebar)
    debar       (incremental step of equivalent strain)
    """
    eps6      = np.zeros(6)
    ebar      = 0.
    time_flow = 0.

    strain=[]; stress=[]
    strain_eq=[]; stress_eq=[]
    time_stamps=[]; strain_rate=[]
    delta_time=[];

    i=0
    while (ebar<ebar_mx):
        ebar      = ebar + debar
        sig_bar   = func_hd(ebar) ## static flow stress
        sig_bar   = sig_bar * func_sr(sr_eq)
        delt_work = debar * sig_bar

        if i==0:
            ## when the initial anisotropy is assumed invariable
            ## Cauchy stress vector direction
            sig_norm6 = alph2sig6(alpha,beta)
            ## equivalent stress state
            sig_norm6 = sig_norm6 / func_yd(sig_norm6)
            ## Strains vector direction (deviatoric strain).
            #deps6 = alph2eps(alpha,beta,func_yd,**yld_kwargs)
            deps6 = alph2eps_c(alpha,beta,func_yd)
            pass

        ## Cauchy stress
        sig6        = sig_norm6*sig_bar
        ## Correct deps6 magnitude by conjugating the incremental work
        dw = 0.
        for j in xrange(3):
            dw = dw + deps6[j]*sig6[j]
            dw = dw + 2*deps6[j+3]*sig6[j+3]
            pass
        # delta time
        dt          = debar  / sr_eq
        x           = delt_work  / dw
        delta_eps6  = deps6      * x
        reps6       = delta_eps6 / dt

        ## stamps
        strain.append(eps6)
        stress.append(sig6)
        strain_eq.append(ebar)
        stress_eq.append(sig_bar)
        strain_rate.append(reps6)
        time_stamps.append(time_flow)
        delta_time.append(dt)

        eps6        = eps6       + delta_eps6
        time_flow   = time_flow  + dt

        i=i+1
        pass

    return np.array(strain),np.array(stress),\
        np.array(strain_eq),np.array(stress_eq),\
        np.array(strain_rate),np.array(time_stamps),\
        np.array(delta_time)


if __name__=='__main__':
    test_FLDA_onepath()
