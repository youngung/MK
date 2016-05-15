import numpy as np
from numba import jit
cos = np.cos
sin = np.sin
log = np.log
tan = np.tan


@jit
def calcD2FB(psi,f2b):
    """
    """
    d2fb = np.zeros((6,6))
    d2fb=np.zeros((6,6))
    cp = cos(psi);  sp = sin(psi)
    c2 = cp*cp; s2 = sp*sp; sc = sp*cp
    d2fb[0,0] = f2b[0,0]*c2+f2b[0,1]*s2+f2b[0,5]*sc/2
    d2fb[1,0] = f2b[1,0]*c2+f2b[1,1]*s2+f2b[1,5]*sc/2
    d2fb[5,0] = f2b[5,0]*c2+f2b[5,1]*s2+f2b[5,5]*sc/2
    d2fb[0,1] = f2b[0,0]*s2+f2b[0,1]*c2-f2b[0,5]*sc/2
    d2fb[1,1] = f2b[1,0]*s2+f2b[1,1]*c2-f2b[1,5]*sc/2
    d2fb[5,1] = f2b[5,0]*s2+f2b[5,1]*c2-f2b[5,5]*sc/2
    d2fb[0,5] = 2*sc*(f2b[0,1]-f2b[0,0])+f2b[0,5]*(c2-s2)/2
    d2fb[1,5] = 2*sc*(f2b[1,1]-f2b[1,0])+f2b[1,5]*(c2-s2)/2
    d2fb[5,5] = 2*sc*(f2b[5,1]-f2b[5,0])+f2b[5,5]*(c2-s2)/2
    return d2fb

# @jit
# def calcD2FB(psi,f2b):
#     """
#     """
#     d2fb = np.zeros((6,6))
#     d2fb=np.zeros((6,6))
#     cp = cos(psi);  sp = sin(psi)
#     c2 = cp*cp; s2 = sp*sp; sc = sp*cp
#     d2fb[0,0] = f2b[0,0]*c2+f2b[0,1]*s2+f2b[0,5]*sc/2
#     d2fb[1,0] = f2b[1,0]*c2+f2b[1,1]*s2+f2b[1,5]*sc/2
#     d2fb[5,0] =(f2b[5,0]*c2+f2b[5,1]*s2+f2b[5,5]*sc)/2
#     d2fb[0,1] = f2b[0,0]*s2+f2b[0,1]*c2-f2b[0,5]*sc/2
#     d2fb[1,1] = f2b[1,0]*s2+f2b[1,1]*c2-f2b[1,5]*sc/2
#     d2fb[5,1] =(f2b[5,0]*s2+f2b[5,1]*c2-f2b[5,5]*sc)/2
#     d2fb[0,5] = 2*sc*(f2b[0,1]-f2b[0,0])+f2b[0,5]*(c2-s2)/2
#     d2fb[1,5] = 2*sc*(f2b[1,1]-f2b[1,0])+f2b[1,5]*(c2-s2)/2
#     d2fb[5,5] =(2*sc*(f2b[5,1]-f2b[5,0])+f2b[5,5]*(c2-s2))/2
#     return d2fb

def func_fld2(ndim,T,s,b,x,yancien,f_hard,f_yld,verbose):
    """
    Arguments
    ---------
    ndim
    T      : axis along integration occurs
    s      : stress of region A
    b      : variables?               [?,f0,?,?,?,?,?,?,deltat]
    x      : the unknowns passed to func_fld [psi,s11,s22,s12]
    yancien: yancien defined in pasapas [deps,]
    f_hard : hardening function
    f_yld  : yield function
    verbose

    Returns
    -------
    F,J,fa,fb,b,siga,s (region A strss)
    """
    import os
    from lib import rot_6d

    f0     = b[1]
    deltat = b[8]

    ## stress in region B
    xb      = x[1]
    yb      = x[2]
    zb      = x[3]
    sb_dump = np.array([xb,yb,0.,0.,0.,zb])

    ## Parameters in region A
    s,phia,fa,f2a = f_yld(s)
    dpsi     = x[0] * (fa[0]-fa[1])*tan(b[9])/(1+tan(b[9])**2)
    b[10]    = dpsi*1.
    psi_new  = b[9]+dpsi
    sa_rot   = rot_6d(s,-psi_new)
    xa,ya,za = sa_rot[0], sa_rot[1],sa_rot[5]
    da_rot   = rot_6d(fa,-psi_new)
    na,ma,siga,dsiga,dma,qqa=f_hard(x[0]+yancien[0])

    ## parameters in region B
    sb             = rot_6d(sb_dump,psi_new)
    sb,phib,fb,f2b = f_yld(sb)
    db_rot         = rot_6d(fb,-psi_new)
    nb,mb,sigb,dsigb,dmb,qqb = f_hard(T+deltat)

    E = -yancien[3] - yancien[4]-deltat* (fb[0]+fb[1])\
        + yancien[1]+yancien[2]\
        + x[0] * (fa[0]+fa[1])

    cp = cos(psi_new);  sp = sin(psi_new)
    c2 = cp*cp; s2 = sp*sp; sc = sp*cp
    d2fb = calcD2FB(psi=psi_new,f2b=f2b)

    dxp    = np.zeros(6)
    dxp[0] = 2*sc*(yb-xb)-2*(c2-s2)*zb
    dxp[1] = -dxp[0]
    dxp[5] = (c2-s2)*(xb-yb)-4*sc*zb
    dpe = (fa[0]-fa[1])*tan(psi_new)/(1+tan(psi_new)**2)
    Q=x[0]/deltat # x[0]: psi

    F=np.zeros(4)
    ## conditions to satisfy
    F[0] = f0*np.exp(E)*sigb*xb - xa*siga*(Q**mb)*qqb**(ma-mb)
    F[1] = xb*za - zb*xa
    F[2] = phib - 1.0
    F[3] = deltat*db_rot[1] - x[0]*da_rot[1]

    J=np.zeros((4,4))
    J[0,0] = (fa[0]+fa[1])*f0*np.exp(E)*sigb*xb-xa*dsiga*(Q**mb)*qqb**(ma-mb)-xa*siga*((mb/deltat)*(Q**(mb-1.))*qqb**(ma-mb)+(dma/x[0])*np.log(qqb)*qqb*(ma-mb))
    J[0,1] =-deltat*(d2fb[0,0]+d2fb[1,0])*f0*np.exp(E)*sigb*xb\
        +f0*np.exp(E)*sigb
    J[0,2] =-deltat*(d2fb[0,1]+d2fb[1,1])*f0*np.exp(E)*sigb*xb
    J[0,3] =-deltat*(d2fb[0,5]+d2fb[1,5])*f0*np.exp(E)*sigb*xb
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
    return F,J,fa,fb,b,siga,s

def func_fld1(ndim,b,x,f_hard,f_yld,verbose):
    """
    Arguments
    ---------
    ndim
    b
    x
    f_hard
    f_yld
    verbose
    """
    import os
    from lib import rot_6d
    psi0, f0, siga, ma, qq, b6 = b[:6]
    # psi0=15*np.pi/180.
    # x[0]=1.
    # x[1]=2.
    # x[3]=0.

    s11,s22,s12 = x[:3]
    s  = np.array([s11,s22,0.,0.,0.,s12])
    sb = rot_6d(s, psi0)

    if verbose:
        print 'psi0:',psi0
        print 'x:'
        print x[:4]
        print 'sb:'
        print sb
        pass

    sb,phib,fb,f2b=f_yld(sb)
    #if verbose:
    if True:
        print 'phib:',phib
        print 'fb:', fb
        print 'f2b:'
        for i in xrange(6):
            for j in xrange(6):
                print '%6.2f'%f2b[i,j],
            print
            pass
        pass

    b[6]=x[3]*fb[0]
    b[7]=x[3]*fb[1]

    db = rot_6d(fb,-psi0)
    if verbose: print 'db:',db

    f2xb = np.zeros((6,6))
    cp = cos(psi0)
    sp = sin(psi0)
    c2 = cp*cp
    s2 = sp*sp
    sc = sp*cp
    f2xb[0,0] = f2b[0,0]*c2 + f2b[0,1]*s2 + f2b[0,5]*sc
    f2xb[1,0] = f2b[1,0]*c2 + f2b[1,1]*s2 + f2b[1,5]*sc
    f2xb[5,0] =(f2b[5,0]*c2 + f2b[5,1]*s2 + f2b[5,5]*sc)/2
    f2xb[0,1] = f2b[0,0]*s2 + f2b[0,1]*c2 - f2b[0,5]*sc
    f2xb[1,1] = f2b[1,0]*s2 + f2b[1,1]*c2 - f2b[1,5]*sc
    f2xb[5,1] =(f2b[5,0]*s2 + f2b[5,1]*c2 - f2b[5,5]*sc)/2
    f2xb[0,5] = 2*(f2b[0,1]-f2b[0,0])*sc + f2b[0,5]*(c2-s2)
    f2xb[1,5] = 2*(f2b[1,1]-f2b[1,0])*sc + f2b[1,5]*(c2-s2)
    f2xb[5,5] =(2*(f2b[5,1]-f2b[5,0])*sc + f2b[5,5]*(c2-s2))/2

    if verbose:
        print 'f2xb:'
        for i in xrange(6):
            for j in xrange(6):
                print '%6.2f'%f2xb[i,j],
                pass
            print
            pass

    nb,mb,sigb,dsigb,dmb,qq=f_hard(x[3])
    if verbose:print 'x(4):',x[3]
    f=np.zeros(4)
    f[0] = db[1]
    f[1] = x[2] - x[0]*b[5]
    f[2] = phib - 1.
    f[3] =-log(b[1]*sigb/b[2])+(b[3]-mb)*log(b[4])-x[3]*db[0]
    if verbose:
        print '--'
        print 'F in func_fld1:'
        print f
        print 'b[:5]'
        print b[:5]
        np.set_printoptions(precision=3)
        print 'db:'
        print db
        # os._exit(1)
        # print '-log(b[1]*sigb/b[2])',-log(b[1]*sigb/b[2])
        # raise IOError
        # print '(b[3]-mb)*log(b[4])-x[3]*db[0]',(b[3]-mb)*log(b[4])-x[3]*db[0]
        pass


    J=np.zeros((4,4))
    J      = np.zeros((4,4))
    J[0,0] = s2  * f2xb[0,0] + c2 * f2xb[1,0] - 2*sc*f2xb[5,0]
    J[0,1] = s2  * f2xb[0,1] + c2 * f2xb[1,1] - 2*sc*f2xb[5,1]
    J[0,2] = s2  * f2xb[0,5] + c2 * f2xb[1,5] - 2*sc*f2xb[5,5]
    J[1,0] = -b[5]
    J[1,2] = 1.
    J[2,0] = db[0]
    J[2,1] = db[1]
    J[2,2] = 2*db[5]
    J[3,3] = -dsigb / sigb - dmb*log(b[4]) - db[0]
    if verbose:
        print 'dmb:',dmb
        print 'dsigb:',dsigb
        print 'sigb:',sigb
        print 'b(5):',b[4]
        print 'ennb:',db[0]
        print "J:"
        for i in xrange(4):
            for j in xrange(4):
                print '%12.2f '%J[i,j],
                pass
            print
            pass
        #print 'in func_fld1'
        #os._exit(1)

    return f, J, fb
