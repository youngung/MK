"""
Adapted from FB's yld2000-2d subroutines for forming limit diagram predictions
"""
from numba import jit
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

def main_deprecated(f0=0.996,fGenPath=None,**kwargs):
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

@jit
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

def main(f0=0.996,psi0=0,th=0,logFileName=None):
    """
    Assumed proportional loadings

    Argument
    --------
    f0
    psi0         in degree
    th  (epsAng) in degree
    """
    # np.seterr(all='raise')
    # np.seterr(all='ignore')

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
    print('stressA:'+('%7.3f'*6)%(
        stressA[0],stressA[1],stressA[2],stressA[3],stressA[4],stressA[5]))
    print('strainA:'+('%7.3f'*6)%(
        dphi[0],dphi[1],dphi[2],dphi[3],dphi[4],dphi[5]))

    print('alpha: %7.4f'%alpha)
    print('rho  : %7.4f'%rho)

    if type(logFileName).__name__=='NoneType':
        logFileName = gen_tempfile(
            prefix='mk-f0%3.3i-th%4.4i-psi%2.2i'%(
                int(f0*1e3),int(th),int(psi0)),
            affix='log')
    logFile  = open(logFileName,'w')

    ## integrate for each path.
    absciss  = 1e3
    absciss0 = 1e3
    print('Iteration over the given psi angle')
    head = ('%8s'*9)%('epsRD','epsTD','psi0','psif','sigRD',
                      'sigTD','sigA','T','cmpT[s]\n')
    logFile.write(head)
    t0   = time.time()

    # ynew, Ahist, Bhist, absciss,xbb,siga,sx = onepath(
    ynew, absciss,xbb,siga,sx = onepath(
        f_yld=f_yld,sa=stressA,
        psi0=psi0*deg2rad,f0=f0,T=absciss)

    dTime = time.time() - t0

    psif1=xbb[0]
    cnt = ('%8.3f'*8+'%8i')%(
        ynew[1],ynew[2],
        psi0,
        psif1*rad2deg,
        sx[0],sx[1],
        siga,
        absciss,dTime)
    print(cnt)
    logFile.write(cnt+'\n')
    uet(dTime,'total time spent');print('')
    logFile.close()
    print('%s has been saved'%logFileName)
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
    print('siga,dsiga,ma,dma')
    print(siga,dsiga,ma,dma)
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


    print ('%7s'+'%11.3f'*6)%('fb   :',fb[0],fb[1],fb[2],
                              fb[3],fb[4],fb[5])
    print ('%7s'+'%20.12e'*5)%('yzero:',yzero[0],yzero[1],yzero[2],
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
    print ('%10.6f'*7)%(xbb[0],xbb[1],xbb[2],
                        xbb[3],xbb[4],xbb[5],xbb[6])

    t0=time.time()

    ## integrate function
    # ynew,Ahist,Bhist,absciss,xbb,siga, SA_fin = pasapas(
    ynew,absciss,xbb,siga, SA_fin = pasapas(
        f0,sa,tzero,yzero,ndds,dydx,xbb,f_hard,f_yld,verbose=False)
    psif = xbb[0]
    print ('%8.3f'*5)%(ynew[0],ynew[1],ynew[2],ynew[3],ynew[4])
    uet(time.time()-t0,'Elapsed time in step by step integration')
    print

    print 'After pasapas'
    print 'absciss:',absciss
    # hist_plot(f_yld,Ahist,Bhist)

    ## check the hardening curve?
    # return ynew,Ahist,Bhist,absciss,xbb,siga, SA_fin
    return ynew,absciss,xbb,siga, SA_fin

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
    verbose

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

    k =-1
    t = tzero
    Ahist=[]; Bhist=[]
    time_used_in_syst=0.
    while(k<=nbpas and absciss<tmax and dydx[0]>=1e-1):
        k=k+1
        ## adjusting the incremental size size
        if dydx[0]<0.5:
            deltt = deltat/10.
            if dydx[0]<0.2:
                deltt = deltat/100.
        else:
            deltt = deltat*1.0

        t0=time.time()
        dydx, ynew, xbb, siga, SA = syst(            
            deltt,t,f0,dydx,xbb,S,yancien,f_hard,f_yld,verbose)
        time_used_in_syst = time_used_in_syst + (time.time()-t0)
        k1     = deltt * dydx ## Y increments
        ynew   = yancien + k1
        t      = t +deltt
        absciss=t*1.
        yancien[::] = ynew[::]
        pass

    uet(time_used_in_syst,'Total time used for iteration in syst')
    return ynew,absciss,xbb,siga, SA

## differential system
def syst(deltt,t,f0,dydx,xbb,sa,y,f_hard,f_yld,verbose):
    """
    Arguments
    ---------
    t       : axis along the intergration occurs
    sa      : stress in region A
    y       : yancien defined in pasapas
    xbb     : [psi0,s1,s2,s3,s4,s5,s6] -- psi0 and stress state of region b
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
    xzero[3] = xbb[6]

    ndim  = 4
    ncase = 2

    bn    = np.zeros(20)
    bn[1] = f0
    bn[8] = deltt
    bn[9] = xbb[0]

    xfinal, fa, fb, bn, siga, sa\
        = new_raph_fld(
            T=t,ndim=ndim,ncase=2,xzero=xzero,y=y,
            b=bn,f_hard=f_hard,f_yld=f_yld,sa=sa,
            verbose=verbose)

    xbb[0] = bn[9]+bn[10]
    xbb[1] = xfinal[1]
    xbb[2] = xfinal[2]
    xbb[5] = xfinal[3]

    dydx[0] = xfinal[0]/deltt
    dydx[1] = fa[0]*dydx[0]
    dydx[2] = fa[1]*dydx[0]
    dydx[3] = fb[0]
    dydx[4] = fb[1]

    # return dydx, y, xbb, regionA[-1],regionB[-1], siga, sa
    return dydx, y, xbb, siga, sa

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
    xn1, fb (case 1)
    xn1, fa, fb, b, siga, sa (case 2)
    """
    import os
    from for_lib import gauss,norme
    from func_fld import func_fld1, func_fld2
    # from func_fld_cy import func_fld1, func_fld2 -- cython was not impressive...

    residu  = 1.0
    xn      = xzero[::]
    it      = 0
    itmax   = 100
    eps     = 1e-10

    totalTimeFunc = 0.
    # if ncase==2: verbose=True ##override
    while (residu>eps and it<itmax):
        it = it+1
        t0 = time.time()
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

            F, J, fa, fb, b,siga, sa \
                = func_fld2(ndim,T,sa,b,xn,y,f_hard,f_yld,verbose)
            dt = time.time() - t0

        totalTimeFunc = totalTimeFunc + dt
        J,res = gauss(ndim=ndim,a=J,b=F) ## f2py module
        xn1=xn-res ## x^(n+1)
        residu = norme(ndim,res)
        xn=xn1[::] ## x^(n)

    # if no convergence
    if it>=itmax:
        print 'could not converge'
        return

    if ncase==1: return xn1,fb
    if ncase==2: return xn1,fa,fb,b,siga,sa

## command line usage (may be called in mk_run for multi-threaded run)
if __name__=='__main__':
    from MP import progress_bar
    import argparse
    uet = progress_bar.update_elapsed_time
    #-------------------------------------------------------
    ## Arguments parsing
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fn',type=str,default='dummy-log-file-name',
        help='File name of the final result')
    parser.add_argument(
        '-f', type=float,default=0.995,
        help='Initial inhomogeneity factor')
    parser.add_argument(
        '-p', type=float,default=0.,
        help='Initial angle of the groove [degree]')
    parser.add_argument(
        '-t', type=float,default=0.,
        help='Angle [theta in degree] of'+\
            ' strain rate theta=atan2(eyy,exx) in [degree]')
    #-------------------------------------------------------
    args = parser.parse_args()
    f0   = args.f
    psi0 = args.p
    th   = args.t
    fn   = args.fn
    main(f0=f0,psi0=psi0,th=th,logFileName=fn)
