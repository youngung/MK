"""
Adapted from FB's yld2000-2d subroutines for forming limit diagram predictions
"""

from for_lib import vm
import numpy as np
import time
from MP import progress_bar
uet=progress_bar.update_elapsed_time
cos=np.cos
sin=np.sin
tan=np.tan
log=np.log
atan2=np.arctan2
sqrt=np.sqrt

def main(f0=0.996,fGenPath=None,**kwargs):
    """
    Assumed proportional loadings

    Argument
    --------
    f0
    fGenPath (None or function or kwargs)
    """
    import os
    from mk_lib   import findStressOnYS
    from lib      import gen_tempfile
    f_yld = vm

    if type(fGenPath).__name__=='NoneType':
        from mk_paths import PSRD
        angs,npt,pth,stressLeft,stressRight = PSRD()
    elif type(fGenPath).__name__=='function':
        angs,npt,pth,stressLeft,stressRight = fGenPath()
    else:
        print type(fGenPath).__name__
        raise IOError, 'unexpected type of fGenPath given'

    ang1,ang2,anginc = angs

    logFileName = gen_tempfile(prefix='log',affix='mk')
    logFile     = open(logFileName,'w')

    tTime = 0.
    for npth in xrange(npt):
        s,dphi = findStressOnYS(
            f_yld,stressLeft.copy(),stressRight.copy(),
            pth=pth,verbose=True)

        ## integrate for each path.
        ntotAngle = ((ang2-ang1)/anginc)+1
        rad2deg   = 180./np.pi
        deg2rad   =   1./rad2deg
        psi0s     = np.linspace(ang1,ang2,ntotAngle)*deg2rad

        absciss  = 1e3
        absciss0 = 1e3
        print 'Iteration over the given psi angle'
        head = ('%8s'*9)%('epsRD','epsTD','psi0','psif','sigRD',
                          'sigTD','sigA','T','cmpTime\n')
        logFile.write(head)

        for psi0_at_each in psi0s:
            print 'PSI: %5.1f'%(psi0_at_each*rad2deg)
            t0   = time.time()
            ynew, Ahist, Bhist, absciss,xbb,siga,sx = onepath(
                f_yld=f_yld,sa=s,psi0=psi0_at_each,f0=f0,T=absciss)
            dTime = time.time() - t0
            tTime = tTime+dTime

            psif1=xbb[0]
            # print ('%8s'*6)%('s1','s2','psi0','psif','siga','absciss')
            print 'ynew:',ynew
            cnt = ('%8.3f'*8+'%8i')%(
                ynew[1],ynew[2],
                psi0_at_each*rad2deg,
                psif1*rad2deg,
                sx[0],sx[1],
                siga,
                absciss,dTime)
            print cnt
            logFile.write(cnt+'\n')

            if absciss<absciss0: ## when a smaller total strain is found update the smallest.
                absciss0=absciss
                psi0_min=psi0_at_each
                psif_min=psif1
                y2      =ynew[1]
                y3      =ynew[2]
                ss1     =sx[0]
                ss2     =sx[1]
                siga_fin=siga

            print 'sigma:',siga
            print '*'*50,'\n'
            pass ## end of each path

        # ## exit
        # logFile.close(); os._exit(1)

        # uet(dt,'elapsed time:')
        # print

        print 'ynew:'
        print(ynew)

        print 'RD strain', 'TD strain', 'Angle psi0', 'angle psif','RD stress','TD stress'
        print ynew[1],ynew[2]
        pass ## end of each path

    uet(tTime,'total time spent')

    # ## exit
    # os._exit(1)
    logFile.close()
    return logFileName,tTime

def calcAlphaRho(s,e):
    """
    """
    if e[0]>e[1]: ## RD
        rho = e[1]/e[0]
        alpha = s[1]/s[0]
    else: ## TD
        rho = e[0]/e[1]
        alpha = s[0]/s[1]
    return rho, alpha


def main2(f0=0.999,psi0=0,th=0):
    """
    Assumed proportional loadings

    Argument
    --------
    f0
    psi0         in degree
    th  (epsAng) in degree
    """
    import os
    from mk_lib   import findStressOnYS
    from lib      import gen_tempfile
    from mk_paths import constructBC,findCorrectPath
    rad2deg   = 180./np.pi
    deg2rad   =   1./rad2deg
    f_yld = vm
    stressA_off, dum1, dum2 = constructBC(epsAng=th, f_yld=f_yld,verbose=False)
    stressA, phi, dphi, d2phi = f_yld(stressA_off) ## put the stress on the locus
    np.set_printoptions(precision=3)
    alpha,rho = calcAlphaRho(stressA,dphi)
    print 'stressA:', stressA
    print 'strainA:', dphi
    print 'alpha:','%3.1f'%alpha
    print 'rho:  ','%3.1f'%rho

    logFileName = gen_tempfile(
        prefix='mk-f0%3.3i-th%4.4i-psi%2.2i'%(
            int(f0*1e3),int(th),int(psi0)),
        affix='log')
    logFile  = open(logFileName,'w')

    ## integrate for each path.
    absciss  = 1e3
    absciss0 = 1e3
    print 'Iteration over the given psi angle'
    head = ('%8s'*9)%('epsRD','epsTD','psi0','psif','sigRD',
                      'sigTD','sigA','T','cmpT[s]\n')
    logFile.write(head)
    t0   = time.time()
    ynew, Ahist, Bhist, absciss,xbb,siga,sx = onepath(
        f_yld=f_yld,sa=stressA,
        psi0=psi0*deg2rad,f0=f0,T=absciss)
    dTime = time.time() - t0

    psif1=xbb[0]
    cnt = ('%8.3f'*8+'%8i')%(
        ynew[1],ynew[2],
        psi0*rad2deg,
        psif1*rad2deg,
        sx[0],sx[1],
        siga,
        absciss,dTime)
    print cnt
    logFile.write(cnt+'\n')
    uet(dTime,'total time spent');print
    logFile.close()
    print '%s has been saved'%logFileName
    return logFileName,dTime

def onepath(f_yld,sa,psi0,f0,T):
    """
    Arguments
    ---------
    f_yld
    sa
    psi0
    f0
    T
    """
    import os
    from lib import rot_6d
    from func_hard_for import return_swift
    # from for_lib import swift
    sa,phia,fa,f2a=f_yld(sa)

    # ## debug
    # sa=np.array([1,2,0,0,0,0])
    # psi0=15.*np.pi/180.

    sx = rot_6d(sa,-psi0)
    # print 'sa:,',sa
    # print 'xa0,ya0,za0'
    # print (3*'%6.3f ')%(sx[0],sx[1],sx[5])
    # return

    na = 0.28985			# n        5e-1
    ks = 518.968			# K       500
    ma = 5e-2                           # strain rate sensitivity
    e0 = 0.0007648 		        # eps_0   1e-5
    qq = 1000. ##Strain rate ratio (E_a.dot / E_0.dot)

    f_hard = return_swift(na,ma,ks,e0,qq)

    # print ks,ma,e0
    na,ma,siga,dsiga,dma,qq = f_hard(0.)
    print 'siga,dsiga,ma,dma'
    print siga,dsiga,ma,dma
    # os._exit(1)

    ## initial_conditions

    ndim = 4
    b    = np.zeros(20)
    b[0] = psi0
    b[1] = f0
    b[2] = siga
    b[3] = ma
    b[4] = qq
    b[5] = sx[5]/sx[0]

    xzero = np.array([1,1,0,0])
    for i in xrange(6):
        print 'B%i'%(i+1),'%7.2f'%b[i]
        pass

    xfinal,fb=new_raph_fld(
        ndim=ndim,ncase=1,
        xzero=xzero,b=b,f_hard=f_hard,f_yld=f_yld,
        verbose=False)
    #np.set_printoptions(precision=3)
    #print 'xfinal:', xfinal
    fmt='%5s          %12.4e'
    print fmt%('e0_b', xfinal[3])
    print fmt%('x_b', xfinal[0])
    print fmt%('y_b', xfinal[1])
    print fmt%('z_b', xfinal[2])

    ## Initial values
    tzero = xfinal[3]
    print 'tzero:',tzero

    yzero = np.zeros(5)
    yzero[0] = 0.
    yzero[1] = 0.
    yzero[2] = 0.
    yzero[3] = tzero*fb[0] ## fb: first derivative in region B
    yzero[4] = tzero*fb[1]

    print ('%7s'+'%8.3f'*5)%('yzero:',yzero[0],yzero[1],yzero[2],
                             yzero[3],yzero[4])

    ndds    = 5 ## dimension differential system
    dydx    = np.zeros(ndds)
    dydx[0] = 1.

    xbb    = np.zeros(7)
    xbb[0] = psi0
    xbb[1] = sx[0]
    xbb[2] = sx[1]
    xbb[6] = sx[5]

    print 'xbb:'
    np.set_printoptions(precision=3)
    print ('%7.3f'*7)%(xbb[0],xbb[1],xbb[2],
                       xbb[3],xbb[4],xbb[5],xbb[6])

    t0=time.time()

    ## integrate function
    ynew,Ahist,Bhist,absciss,xbb,siga, SA_fin = pasapas(
        f0,sa,tzero,yzero,ndds,dydx,xbb,f_hard,f_yld,verbose=False)
    psif = xbb[0]
    print ('%8.3f'*5)%(ynew[0],ynew[1],ynew[2],ynew[3],ynew[4])
    uet(time.time()-t0,'Elapsed time in step by step integration')
    print

    print 'After pasapas'
    print 'absciss:',absciss
    # hist_plot(f_yld,Ahist,Bhist)

    ## check the hardening curve?
    return ynew,Ahist,Bhist,absciss,xbb,siga, SA_fin

def hist_plot(f_yld,Ahist,Bhist):
    """
    Arguments
    ---------
    f_yld (common for both regions)
    Ahist
    Bhist
    """
    EA=[]; SA=[]; EB=[]; SB=[]
    if (len(Ahist)!=len(Bhist)):
        raise IOError, 'Unexpected unbalanced step size'

    sigma_A=[]; sigma_B=[]
    for i in xrange(len(Ahist)):
        A = Ahist[i]; B = Bhist[i]
        ea = A.H.eps; sa = A.H.sig
        eb = B.H.eps; sb = B.H.sig

        sig_A = A.stress; sig_B = B.stress
        sigma_A.append(sig_A); sigma_B.append(sig_B)

        EA.append(ea); EB.append(eb)
        SA.append(sa); SB.append(sb)

    EA=np.array(EA); EB=np.array(EB)
    SA=np.array(SA); SB=np.array(SB)
    S6A=np.array(sigma_A); S6B=np.array(sigma_B)

    import matplotlib.pyplot as plt
    fig=plt.figure(figsize=(10,3.5));
    ax1=fig.add_subplot(131)
    ax2=fig.add_subplot(132)
    ax3=fig.add_subplot(133)

    ax1.plot(EA,SA,label='A',ls='-',zorder=99)
    ax1.plot(EB,SB,label='B',ls='-',zorder=100,alpha=0.4)

    ## plot yield locus
    pi = np.pi; sin=np.sin; cos=np.cos
    th = np.linspace(-pi,pi)
    x=cos(th);y=sin(th)
    z=np.zeros(len(th))
    s=np.array([x,y,z,z,z,z]).T
    print s.shape
    X=[]; Y=[]
    for i in xrange(len(s)):
        ys, phi, dphi, d2phi = vm(s[i])
        X.append(ys[0])
        Y.append(ys[1])

    X=np.array(X)
    Y=np.array(Y)
    ax2.plot(X,Y,label='Yield locus')
    ## initial location of region A stress state
    ax2.plot(S6A[0][0],S6A[0][1],'r.',mfc='None',mec='r',label='A initial')
    ## final location of region A stress state
    ax2.plot(S6A[-1][0],S6A[-1][1],'rx',mfc='None',mec='r',label='A final')
    ## initial location of region B stress state
    ax2.plot(S6B[0][0],S6B[0][1],'g.',mfc='None',mec='g',label='B initial')
    ## final location of region B stress state
    ax2.plot(S6B[-1][0],S6B[-1][1],'gx',label='B final')

    A1=X*SA[-1];A2=Y*SA[-1]
    B1=X*SB[-1];B2=Y*SB[-1]
    ax3.plot(A1,A2,'-',label='Final Yield locus (A)')
    ax3.plot(B1,B2,'-',label='Final Yield locus (B)')
    # print 'A1'
    # print(A1)
    # print 'B1'
    # print(B1)
    ## initial location of region A stress state
    ax3.plot(S6A[0][0]*SA[0],S6A[0][1]*SA[0],'r.',mfc='None',mec='r',label='A initial')
    ## final location of region A stress state
    ax3.plot(S6A[-1][0]*SA[-1],S6A[-1][1]*SA[-1],'rx',mfc='None',mec='r',label='A final')
    ## initial location of region B stress state
    ax3.plot(S6B[0][0]*SB[0],S6B[0][1]*SB[0],'g.',label='B initial')
    ## final location of region B stress state
    ax3.plot(S6B[-1][0]*SB[-1],S6B[-1][1]*SB[-1],'gx',label='B final')
    ax2.legend();ax3.legend()

    fn='hist_plot.pdf'
    fig.tight_layout()
    fig.savefig(fn,bbox_inches='tight')
    print '%s has been saved'%fn

def pasapas(f0,S,tzero,yzero,ndds,dydx,xbb,f_hard,f_yld,verbose):
    """
    step by step integration

    f0  : initial inhomgeneity factor
    S   : stress state of region A
    tzero
    yzero
    ndds: dimension differential system
    dydx:
    xbb :
    f_hard : strain hardening function
    f_yld  : yield function
n    verbose

    Returns
    -------
    ynew
    Ahist
    Bhist
    absciss (T)
    xbb
    siga
    sa
    """
    import os
    absciss = tzero ## xcoordinate in the intergration
    yancien = np.zeros(ndds) ## y_old
    yancien = yzero[:]

    ## integration values
    nbpas  = 200000
    freq   = 100      # frequency of output (probably not needed for this)
    deltat = 0.001    # (stepsize)
    tmax   = 1.5      # (maximum effective strain upper limit (2.0?))
    class int_opt:
        def __init__(self):
            self.nbpas  = nbpas
            self.freq   = freq      # frequency of output (probably not needed for this)
            self.deltat = deltat    # (stepsize)
            self.tmax   = tmax      # (maximum effective strain upper limit (2.0?))

    k =-1
    t = tzero
    print 'nbpas:',nbpas
    print 'tmax:', tmax
    print 'dydx(1):',dydx[0]

    Ahist=[]; Bhist=[]
    time_used_in_syst=0.
    while(k<=nbpas and absciss<tmax and dydx[0]>=1e-1):
        k=k+1
        ## adjusting the incremental size size?
        if   dydx[0]<0.2: deltt = deltat/100.
        elif dydx[0]<0.5: deltt = deltat/10.
        else:             deltt = deltat
        t0=time.time()
        dydx, ynew, xbb, regA, regB, siga, SA = syst(
            deltt,t,f0,dydx,xbb,S,yancien,f_hard,f_yld,verbose)

        time_used_in_syst = time_used_in_syst + (time.time()-t0)
        Ahist.append(regA)
        Bhist.append(regB)
        k1   = deltt * dydx ## Y increments
        ynew = yancien + k1
        t    = t +deltt

        # np.set_printoptions(precision=6)
        # print '------------------------------------------------------------'
        # print ('%2i '+'%11.4e '*5)%(k, ynew[0],ynew[1],ynew[2],ynew[3],ynew[4])
        # np.set_printoptions(precision=3)

        absciss=t
        yancien[::] = ynew[::]

        # if k==3:
        #     import os
        #     print 'Stop in pasapas for debug'
        #     os._exit(1)

        pass

    uet(time_used_in_syst,'Total time used for iteration in syst')
    print
    print '-'*70
    print 'k,ynew resulting from pasapas'
    print ('%2i '+'%11.3f '*5)%(k, ynew[0],ynew[1],ynew[2],ynew[3],ynew[4])
    print '-'*70
    print 'Stop in paspas for debug'
    print 'End of loop in pasapas ----'
    # os._exit(1)
    return ynew,Ahist,Bhist,absciss,xbb,siga, SA

## differential system
def syst(deltt,t,f0,dydx,xbb,sa,y,f_hard,f_yld,verbose):
    """
    Arguments
    ---------
    t       : axis along the intergration occurs
    sa      : stress in region A
    y       : yancien defined in pasapas
    xbb     : [psi0,s1,s2,s6]
    f_hard  : strain (and strain rate) hardening function
    f_yld   : yield function
    verbose

    Returns
    -------
    dydx
    yancien
    siga     strain hardening flow stress, sig = hard(E)
    sa        stress of region a
    """
    import os
    xzero    = np.zeros(4)
    xzero[0] = dydx[0]*deltt
    xzero[1] = xbb[1]
    xzero[2] = xbb[2]
    xzero[3] = xbb[3]

    ndim  = 4
    ncase = 2

    bn    = np.zeros(20)
    bn[1] = f0
    bn[8] = deltt
    bn[9] = xbb[0]

    xfinal, fa, fb, bn, regionA, regionB, siga, sa\
        = new_raph_fld(
            T=t,ndim=ndim,ncase=2,xzero=xzero,y=y,
            b=bn,f_hard=f_hard,f_yld=f_yld,sa=sa,
            verbose=verbose)
    ##
    xbb[0] = bn[9]+bn[10]
    xbb[1] = xfinal[1]
    xbb[2] = xfinal[2]
    xbb[5] = xfinal[3]

    dydx[0] = xfinal[0]/deltt
    dydx[1] = fa[0]*dydx[0]
    dydx[2] = fa[1]*dydx[0]
    dydx[3] = fb[0]
    dydx[4] = fb[1]

    np.set_printoptions(precision=6)

    # print 'sa:'
    # print sa
    # print 'fa'
    # print(fa)
    # print 'fb'
    # print(fb)

    # print 'xbb:'
    # print(xbb)
    # print 'dydx'
    # print dydx[:5]

    np.set_printoptions(precision=3)

    # raise IOError
    return dydx, y, xbb, regionA[-1],regionB[-1], siga, sa

def new_raph_fld(
        T=None,ndim=None,ncase=None,
        xzero=None,y=None,b=None,f_hard=None,f_yld=None,
        sa=None,verbose=True):
    """
    Arguments
    ---------
    T
    ndim
    ncase
    xzero
    y
    b
    f_hard
    verbose=True

    Return
    ------
    """
    import os
    if ncase==2 and verbose:
        print 'sa:'
        print sa
        print 'T:'
        print T
        print 'b:'
        print b
        np.set_printoptions(precision=5)
        print 'y:'
        print y
        np.set_printoptions(precision=3)

    from for_lib import gauss,norme
    residu = 1.0
    # xn    = np.zeros(ndim)
    xn    = xzero[::]
    it    = 0
    itmax = 100
    eps   = 1e-10
    A_region=[]; B_region=[]

    totalTimeFunc = 0.
    while (residu>eps and it<itmax):
        it = it+1
        t0=time.time()
        if ncase==1:
            if verbose:
                print '-'*40
                print '%i ITERATION over func_fld1 in NR'%it
                print 'xn:'
                print xn
            F, J, fb        = func_fld1(ndim,b,xn,f_hard,f_yld,verbose)
            dt = time.time() - t0
        if ncase==2:
            if verbose:
                print '-'*40
                print '%i ITERATION over func_fld2 in NR'%it
            F, J, fa, fb, b, mat_A, mat_B, siga, sa \
                = func_fld2(ndim,T,sa,b,xn,y,f_hard,f_yld,verbose)
            dt = time.time() - t0
            A_region.append(mat_A)
            B_region.append(mat_B)

        totalTimeFunc = totalTimeFunc + dt

        if ncase==1 and verbose:
            print 'result prior to gauss'
            print 'F:'
            print F
            print 'J:'
            print J
        J,res = gauss(ndim=ndim,a=J,b=F)
        if ncase==1 and verbose:
            print 'result after gauss'
            print 'F:'
            print F
            print 'J:'
            print J
            print 'res:'
            print res
        elif ncase==2 and verbose:
            print 'F:'
            print F
            print 'J:'
            print J
            print 'res:'
            print res


        xn1=xn-res ## x^(n+1)
        residu = norme(ndim,res)
        xn=xn1[::] ## x^(n)
        if ncase==2 and verbose:
            np.set_printoptions(precision=5)
            print '%3i %13.4e'%(it,residu),xn
            np.set_printoptions(precision=3)
        pass

    if ncase==2 and verbose:
        print it,
        print 'b:'
        print b[:10]
        print 'xn in new_raph_fld ncase==2:'
        print xn

    # if no convergence
    if it>=itmax:
        print 'could not converge'
        return

    # if ncase==2:
    #     uet(totalTimeFunc, 'Total time spent for func_FLD in NR')
    # print

    if ncase==1: return xn1,fb
    if ncase==2: return xn1,fa,fb,b,A_region,B_region,siga,sa

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
    F,J,fa,fb,b,mat_A,mat_B, siga, sa
    """
    import os
    from lib import rot_6d

    class material:
        class H:
            pass
        class Y:
            pass

    mat_A = material()
    mat_B = material()

    f0     = b[1]
    deltat = b[8]

    ## stress in region B
    xb      = x[1]
    yb      = x[2]
    zb      = x[3]
    sb_dump = np.array([xb,yb,0,0,0,zb])

    ## Parameters in region A
    sa,phia,fa,f2a = f_yld(s)
    mat_A.Y.phi   = phia
    mat_A.Y.dphi  = fa
    mat_A.Y.d2phi = f2a
    mat_A.stress  = sa

    dpsi     = x[0] * (fa[0]-fa[1])*tan(b[9])/(1+tan(b[9])**2)
    b[10]    = dpsi
    psi_new  = b[9]+b[10]
    sa_rot   = rot_6d(sa,-psi_new)
    xa,ya,za = sa_rot[0], sa_rot[1],sa_rot[5]
    da_rot   = rot_6d(fa,-psi_new)
    na,ma,siga,dsiga,dma,qqa=f_hard(x[0]+yancien[0])
    mat_A.H.n    = na
    mat_A.H.m    = ma
    mat_A.H.eps  = x[0]+yancien[0]
    mat_A.H.sig  = siga
    mat_A.H.dsig = dsiga
    mat_A.H.dm   = dma
    mat_A.H.qq   = qqa
    mat_A.psi    = psi_new

    if verbose:
        print 'ma,siga,dsiga,dma,qqa'
        print ma,siga,dsiga,dma,qqa

    ## parameters in region B
    sb             = rot_6d(sb_dump,psi_new)
    sb,phib,fb,f2b = f_yld(sb)
    mat_B.Y.phi   = phib
    mat_B.Y.dphi  = fb
    mat_B.Y.d2phi = f2b
    mat_B.stress  = sb
    db_rot        = rot_6d(fb,-psi_new)
    nb,mb,sigb,dsigb,dmb,qqb = f_hard(T+deltat)
    mat_B.H.n = nb
    mat_B.H.m = mb
    mat_B.H.eps = T+deltat
    mat_B.H.sig = sigb
    mat_B.H.dsig = dsigb
    mat_B.H.dm   = dmb
    mat_B.H.qq   = qqb
    mat_B.psi    = psi_new
    if verbose:
        print 'mb,sigb,dsigb,dmb,qqb'
        print mb,sigb,dsigb,dmb,qqb
    E = -yancien[3] - yancien[4]-deltat* (fb[0]+fb[1])\
        + yancien[1]+yancien[2]\
        + x[0] * (fa[0]+fa[1])
    if verbose:
        print 'E:'
        print '%11.4e'%E

    d2fb=np.zeros((6,6))
    cp = cos(psi_new);  sp = sin(psi_new)
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
    if verbose:
        print 'd2fb:'
        print(d2fb)

    dxp    = np.zeros(6)
    dxp[0] = 2*sc*(yb-xb)-2*(c2-s2)*zb
    dxp[1] =-dxp[0]
    dxp[5] = (c2-s2)*(xb-yb)-4*sc*zb
    if verbose:
        print 'dxp:'
        print(dxp)

    dpe = (fa[0]-fa[1])*tan(psi_new)/(1+tan(psi_new)**2)
    if verbose:
        print 'dpe:'
        print(dpe)
    Q=x[0]/deltat
    F=np.zeros(4)

    ## conditions to satisfy
    F[0] = f0*np.exp(E)*sigb*xb-xa*siga*(Q**mb)*qqb**(ma-mb)
    F[1] = xb*za - zb*xa
    F[2] = phib - 1.0
    F[3] = deltat*db_rot[1] - x[0]*da_rot[1]

    if verbose:
        print 'F:'
        print(F)

    J=np.zeros((4,4))
    J[0,0] = (fa[0]+fa[1])*f0*np.exp(E)*sigb*xb
    if verbose:
        print J[0,0]

    J[0,0] = J[0,0] -xa*dsiga*(Q**mb)*qqb**(ma-mb)
    dum1=(mb/deltat)*(Q**(mb-1.))*qqb**(ma-mb)
    if verbose:
        print'Q:',Q
        print 'QQb:',qqb
        print 'dum1',dum1
    dum2=(dma/x[0])*np.log(qqb)*qqb*(ma-mb)
    dum=dum1+dum2

    if verbose:
        print 'dum2',dum2
        print 'dum',dum
    J[0,0] = J[0,0] -xa*siga*(dum)
    """J(1,1) = J(1,1) -XA*DSIGA*(Q**MB)*QQ**(MA-MB)
             -XA*SIGA*((MB/DELTAT)*(Q**(MB-1))*QQ**(MA-MB)
                       +(DMA/X(1))*DLOG(QQ)*QQ*(MA-MB))"""
    if verbose:
        print J[0,0]
        print qqb

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
    if verbose:
        print 'J:'
        print(J)
    return F, J, fa, fb, b, mat_A, mat_B, siga, sa

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
    s = np.array([s11,s22,0,0,0,s12])
    sb = rot_6d(s, psi0)

    if verbose:
        print 'psi0:',psi0
        print 'x:'
        print x[:4]
        print 'sb:'
        print sb
        pass

    sb,phib,fb,f2b=f_yld(sb)
    if verbose:
        print 'phib:',phib
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

## command line usage
if __name__=='__main__':
    main2(f0=0.996,psi0=5,th=0)
    # from lib import gen_tempfile
    # from mk_paths import returnPaths
    # DRD,PSRD,BBRD,BBTD,PSTD,DTD = returnPaths()
    # f0 = 0.996
    # logFN, tTime = main(f0,DRD)
    # print logFN
    # uet(tTime,'Time elapsed for the given path')
