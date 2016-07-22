"""
Functions used in forming limit diagrams
"""
import numpy as np
from numba import jit
cos = np.cos
sin = np.sin
log = np.log
tan = np.tan


## @jit
def calcD2FB(psi,f2b):
    """
    """
    import time
    t0_grand = time.time()
    d2fb = np.zeros((6,6))

    cp = cos(psi);  sp = sin(psi)
    c2 = cp*cp; s2 = sp*sp; sc = sp*cp
    t1 = time.time()-t0_grand


    t0 = time.time()
    d2fb[0,0] = f2b[0,0]*c2+f2b[0,1]*s2+f2b[0,5]*sc/2.
    d2fb[1,0] = f2b[1,0]*c2+f2b[1,1]*s2+f2b[1,5]*sc/2.
    d2fb[5,0] = f2b[5,0]*c2+f2b[5,1]*s2+f2b[5,5]*sc/2.
    d2fb[0,1] = f2b[0,0]*s2+f2b[0,1]*c2-f2b[0,5]*sc/2.
    d2fb[1,1] = f2b[1,0]*s2+f2b[1,1]*c2-f2b[1,5]*sc/2.
    d2fb[5,1] = f2b[5,0]*s2+f2b[5,1]*c2-f2b[5,5]*sc/2.
    d2fb[0,5] = 2.*sc*(f2b[0,1]-f2b[0,0])+f2b[0,5]*(c2-s2)/2.
    d2fb[1,5] = 2.*sc*(f2b[1,1]-f2b[1,0])+f2b[1,5]*(c2-s2)/2.
    d2fb[5,5] = 2.*sc*(f2b[5,1]-f2b[5,0])+f2b[5,5]*(c2-s2)/2.
    t2=time.time()-t0
    dt = time.time()-t0_grand

    # print ('%.1f  '*2)%(t1/dt*100, t2/dt*100)
    return d2fb

# @jit
def calcF2XB(psi0,f2b):
    """
    Originally in func_fld1 function
    Its role is to rotate the second derivative of
    yield function (i.e., <f2b>) into the band axes
    """
    cp = cos(psi0)
    sp = sin(psi0)
    c2 = cp*cp
    s2 = sp*sp
    sc = sp*cp
    f2xb = np.zeros((6,6))
    f2xb[0,0] = f2b[0,0]*c2 + f2b[0,1]*s2 + f2b[0,5]*sc
    f2xb[1,0] = f2b[1,0]*c2 + f2b[1,1]*s2 + f2b[1,5]*sc
    f2xb[5,0] =(f2b[5,0]*c2 + f2b[5,1]*s2 + f2b[5,5]*sc)/2
    f2xb[0,1] = f2b[0,0]*s2 + f2b[0,1]*c2 - f2b[0,5]*sc
    f2xb[1,1] = f2b[1,0]*s2 + f2b[1,1]*c2 - f2b[1,5]*sc
    f2xb[5,1] =(f2b[5,0]*s2 + f2b[5,1]*c2 - f2b[5,5]*sc)/2
    f2xb[0,5] = 2*(f2b[0,1]-f2b[0,0])*sc + f2b[0,5]*(c2-s2)
    f2xb[1,5] = 2*(f2b[1,1]-f2b[1,0])*sc + f2b[1,5]*(c2-s2)
    f2xb[5,5] =(2*(f2b[5,1]-f2b[5,0])*sc + f2b[5,5]*(c2-s2))/2
    return f2xb

def func_fld2(
        ndim,
        T,
        # s,
        b,
        x,
        yold,
        matA,
        matB,
        verbose):
    """
    Arguments
    ---------
    ndim
    T      : axis along integration occurs
    s      : stress of region A
    b      : variables?

        b[0] : unknown
        b[1] : f0
        b[2] : unknown
        b[3] : unknown
        b[4] : unknown
        b[5] : unknown
        b[6] : unknown
        b[7] : unknown
        b[8] : delta t:   delta equivalent strain increment: \dot{lambda}d(time)=dleta lambda^B for region B
        b[9] : psi_old
        b[10]: psi_new

    x      :
       x[0]: d\lambda^A, i.e., the increment of equivalent strain pertaining to region A
       x[1]: s11 component of stress
       x[2]: s22 component of stress
       x[3]: s12 component of stress

    yold: yold defined in pasapas [deps,]
        y[1]: accumulative strain RD
        y[2]: accumulative strain TD
    matA   : A material
    matB   : B material
    verbose

    ** yold
    Initially, it is determined:
    [0,0,0,tzero*fb[0],tzero*fb[1]]

    ** bn
    b[1]  = f0
    b[8]  = deltt
    b[9]  = psi_old
    b[10] = delta_psi

    Returns
    -------
    F,J,fa,fb,b,siga,s (region A strss)
    """
    import os, time
    from mk.library.lib import rot_6d

    t0 = time.time()
    t0_grand = time.time()

    f0     = b[1]
    deltat = b[8] # delta t for region B
    ## stress in region B
    xb      = x[1]
    yb      = x[2]
    zb      = x[3]
    ## (guessed) stress referred in band
    sb_dump = np.array([xb,yb,0.,0.,0.,zb])

    ## Parameters in region A
    matA.update_yld(matA.stress)
    s,phia,fa,f2a = matA.o_yld

    t1 = time.time()-t0_grand
    # ------------------------------------------------------------#
    t0 = time.time()

    ## psi0 * ()
    delta_psi= x[0] * (fa[0]-fa[1])*tan(b[9])/(1+tan(b[9])**2)
    b[10]    = delta_psi*1.
    psi_new  = b[9] + delta_psi
    sa_rot   = rot_6d(s,-psi_new) ## A stress referred in the band axes
    xa,ya,za = sa_rot[0], sa_rot[1],sa_rot[5]
    da_rot   = rot_6d(fa,-psi_new) ## A's dphi/dsig referred in the band axes
    matA.update_hrd(yold[0]+x[0])
    # na,ma,siga,dsiga,dma,qqa = matA.o_hrd
    ma, siga, dsiga, dma, qqa = matA.o_hrd[:5]

    ## parameters in region B
    sb      = rot_6d(sb_dump,psi_new) ## use updated (tilde) psi
    matB.update_yld(sb)
    sb,phib,fb,f2b = matB.o_yld
    db_rot         = rot_6d(fb,-psi_new)
    matB.update_hrd(T+deltat) ## delta t (incremental effective strain)
    mb,sigb,dsigb,dmb,qqb = matB.o_hrd[:5]

    E = -yold[3]-yold[4]-deltat*(fb[0]+fb[1])\
        +yold[1]+yold[2]+  x[0]*(fa[0]+fa[1])

    t2 = time.time()-t0
    # ------------------------------------------------------------#
    t0 = time.time()

    cp = cos(psi_new);  sp = sin(psi_new)
    c2 = cp*cp; s2 = sp*sp; sc = sp*cp

    t_bench_0 = time.time()
    d2fb = calcD2FB(psi=psi_new,f2b=f2b)
    dt_bench = time.time() - t_bench_0

    dxp    = np.zeros(6)
    dxp[0] = 2*sc*(yb-xb)-2*(c2-s2)*zb
    dxp[1] = -dxp[0]
    dxp[5] = (c2-s2)*(xb-yb)-4*sc*zb
    dpe = (fa[0]-fa[1])*tan(psi_new)/(1+tan(psi_new)**2)
    Q   = x[0]/deltat # ratio between A / B equivalent strain rate increment?

    ## x[0] = \delta \Lambda^A?
    ## x[1] = xb
    ## x[2] = yb ???
    ## x[3] = zb


    t3 = time.time()-t0
    # print 'dt_bench/t3',(dt_bench/t3)*100
    # ------------------------------------------------------------#

    t0 = time.time()
    F = np.zeros(4)
    ## conditions to satisfy

    eE = np.exp(E)

    F[0] = f0*eE*sigb*xb - xa*siga*(Q**mb)*qqb**(ma-mb)

    F[1] = xb*za - zb*xa
    F[2] = phib - 1.0
    F[3] = deltat*db_rot[1] - x[0]*da_rot[1]

    J=np.zeros((4,4))
    J[0,0] =  (fa[0]+fa[1])*f0*eE*sigb*xb\
              -xa*dsiga*(Q**mb)*qqb**(ma-mb)\
              -xa*siga*(
                  (mb/deltat)*(Q**(mb-1.))*qqb**(ma-mb)\
                  +(dma/x[0])*np.log(qqb)*qqb*(ma-mb)
              )
    J[0,1] =-deltat*(d2fb[0,0]+d2fb[1,0])*f0*eE*sigb*xb\
        +f0*eE*sigb
    J[0,2] =-deltat*(d2fb[0,1]+d2fb[1,1])*f0*eE*sigb*xb
    J[0,3] =-deltat*(d2fb[0,5]+d2fb[1,5])*f0*eE*sigb*xb
    J[1,1] = za
    J[1,3] =-xa

    J[2,0] = (fb[0]*dxp[0]+fb[1]*dxp[1]+2*fb[5]*dxp[5])*dpe
    J[2,1] = db_rot[0]
    J[2,2] = db_rot[1]
    J[2,3] = 2*db_rot[5]

    J[3,0] = (f2b[0,0]*dxp[0]+f2b[0,1]*dxp[1]+f2b[0,5]*dxp[5]/2.)*s2+(f2b[1,0]*dxp[0]+f2b[1,1]*dxp[1]+f2b[1,5]*dxp[5]/2)*c2
    J[3,0] = (J[3,0]-2*(f2b[5,0]*dxp[0]+f2b[5,1]*dxp[1]+f2b[5,5]*dxp[5]/2)*sc)*deltat*dpe
    J[3,0] = J[3,0]-da_rot[1]+x[0]*(2*sc*(fa[0]-fa[1])-2*fa[5]*(c2-s2))*dpe + deltat*(2*sc*(fb[0]-fb[1])-2*fb[5]*(c2-s2))*dpe
    J[3,1] = deltat*(d2fb[0,0]*s2+d2fb[1,0]*c2-2*d2fb[5,0]*sc)
    J[3,2] = deltat*(d2fb[0,1]*s2+d2fb[1,1]*c2-2*d2fb[5,1]*sc)
    J[3,3] = deltat*(d2fb[0,5]*s2+d2fb[1,5]*c2-2*d2fb[5,5]*sc)


    t4 = time.time()-t0
    # ------------------------------------------------------------#
    dt_grand = time.time() - t0_grand

    t1=t1/dt_grand * 100
    t2=t2/dt_grand * 100
    t3=t3/dt_grand * 100
    t4=t4/dt_grand * 100

    # print '---------------------------------------'
    # print ('%.1f '*5)%(t1,t2,t3,t4, (t1+t2+t3+t4))
    # print '---------------------------------------'
    # raise IOError

    return F,J,fa,fb,b# ,s

def func_fld1(
        ndim,b,x,
        # f_hard,f_yld,
        matA, matB,
        verbose):
    """
    Arguments
    =========
    ndim
    ----

    b
    -
    an array determined in the onepath
      that is initially set as:
      [ psi0, f0, matA.sig, matA.m, matA.qq, sx[5]/sx[0] ]
      sx is initially the stress of A that is rotated by psi
    b is *NOT* iteratively changing within the Newton Raphson
      as func_fld1 is to find the initial state of region B
      that corresponding to the initial state of region A
    b[6], b[7] are changing though - they are functions of fb[0], fb[1]


    x
    -
    initially given as [1,1,0,0] but should be
    iteratively determined
    x[0] = s11; x[1] == s22 and x[2] == s12 for the region B
    x[3] is unknown yet, but seems related with some kind of
        incremental step size (?)
        - potentially x[3] is d\lambda, that is the plastic multiplier.
        - Usually, d\lambda determines the incremental
          size of plastic strain rate

    F_{i} = J_{ij}: x_{j}


    matA
    ----
    Material description of region A
    (including the evolved state varaiables)


    matB
    ----
    Material description of region B
    (including the evolved state varaiables)


    verbose


    Returns
    -------
    F   : objective functions
    J   : jacobian
    fb  : first derivative of the yield function
          pertaining to the region B)
    """
    import os
    from mk.library.lib import rot_6d
    psi0, f0, siga, ma, qq, b6 = b[:6]

    ### currently (guessed) stress states of region B in band axes
    s11,s22,s12 = x[:3]
    s = np.array([s11,s22,0.,0.,0.,s12])

    ## rotate it back to 'RD/TD/ND' axes
    ## This is primarily due to the fact that yield functions
    ## are usually written in the frame of rolled sheet metals.
    ## The stress (as the argument) to the yield function should
    ## be referred in the proper material axes.
    sb = rot_6d(s, psi0)
    ## stress of region B in the rd-td-nd axes
    sb, phib, fb, f2b = matB.f_yld(sb)

    """
    For reminding the defintion of b-array:
    b[0] = psi0
    b[1] = f0
    b[2] = matA.sig
    b[3] = matA.m
    b[4] = matA.qq
    ## stress state ratio within the band
    ## from region A stress state
    b[5] = sx[5]/sx[0]

    if x[3] is d\lambda(
       that means, in the context of func_fld1,
       the first incremental strain.
    b[6] = $\dot{e}_{RD}^B$
    b[7] = $\dot{e}_{TD}^B$
    """
    b[6]=x[3]*fb[0]
    ## equivalent accumulative strain x
    ## (or the first incremental equivalent strain)
    b[7]=x[3]*fb[1]

    db = rot_6d(fb,-psi0) ## Region B's strain rate in the band axes
    if verbose: print 'db:',db

    ## Seems to rotate the second derivative of yield function
    ## into the band axes
    ## May not be required when psi0 = 0,
    ## which is often the case in forming limit simulations
    f2xb = calcF2XB(psi0,f2b)

    mb,sigb,dsigb,dmb,qq=matB.f_hrd(x[3])[:5]
    if verbose:print 'x(4):',x[3]
    f=np.zeros(4)
    f[0] = db[1]              ## ettb (transverse strain rate or region B, should it be ettb - etta???)
    ## force equilibrium condition
    f[1] = x[2] - x[0]*b[5]   ## s12b - s11b * s12a/s11a ( all in band axes)
    f[2] = phib - 1.          ## determine if the stress is on the yield locus (consistency condition?)
    f[3] =-log(b[1]*sigb/b[2])+(b[3]-mb)*log(b[4])-x[3]*db[0]


    if not(np.isfinite(f[3])):
        print '-'*20
        print 'f[3] is not finite:',f[3]
        print 'b[:]:',b[:]
        print 'mb:',mb
        print 'sigb:',sigb
        print 'dsigb:',dsigb
        print 'dmb:',dmb
        print 'qq:',qq
        print 'x[:]',x[:]
        print '-'*20


    cp = cos(psi0)
    sp = sin(psi0)
    c2 = cp*cp
    s2 = sp*sp
    sc = sp*cp

    J=np.zeros((4,4))
    J      = np.zeros((4,4))
    J[0,0] = s2  * f2xb[0,0] + c2 * f2xb[1,0] - 2*sc*f2xb[5,0]
    J[0,1] = s2  * f2xb[0,1] + c2 * f2xb[1,1] - 2*sc*f2xb[5,1]
    J[0,2] = s2  * f2xb[0,5] + c2 * f2xb[1,5] - 2*sc*f2xb[5,5]
    J[1,0] = -b[5]   ## sx[5]/sx[0]
    J[1,2] = 1.
    J[2,0] = db[0]   ## e_nn^B
    J[2,1] = db[1]   ## e_tt^B
    J[2,2] = 2*db[5] ## 2e_tn^B
    J[3,3] = -dsigb / sigb - dmb*log(b[4]) - db[0] ## (dH^B/H^B) - dm^B (qq^A) - enn^B

    return f, J, fb
