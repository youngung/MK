"""
Forming limig diagram predictive tool on the basis of
macro-mechanical constitutive models using anisotropic yield function
and hardening curve.

Adapted from FB's forming limit calculation subroutine/algorithm

Youngung Jeong
--------------
youngung.jeong@gmail.com
younguj@clemson.edu
--------------------------------------------
International Center for Automotive Research
Clemson University, Greenville, SC
"""
import matplotlib as mpl
mpl.use('Agg')            ## In case X-window is not available.
from numba import jit
from yf_for import vm
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

def main(
        f0=0.996,
        psi0=0,
        th=0,
        material=None,
        logFileName=None):
    """
    Run forming limit test for the given path
    using the given material (if none given, assume isotropic material)

    Arguments
    ---------
    f0           initial inhomogeneity factor
    psi0         [degree]
    th  (epsAng) [degree]
    material    = None (a material data from constitutive.Constitutive)
    logFileName = None
    """
    # np.seterr(all='raise')
    np.seterr(all='ignore')
    import os
    from mk.library.mk_lib   import findStressOnYS
    from mk.library.lib      import gen_tempfile, calcAlphaRho
    from mk_paths import constructBC,findCorrectPath
    import mk.materials.constitutive as constitutive
    import dill
    snapshot = constitutive.Snapshot()
    # from yf2 import wrapHill48

    print 'material:',material, type(material).__name__

    if type(material).__name__=='NoneType':
        print 'given material', material
        from materials.materials import IsoMat
        matA = IsoMat()
        matB = IsoMat()
    elif type(material).__name__=='str':
        with open(material,'rb') as fo:
            matA = dill.load(fo)
        with open(material,'rb') as fo:
            matB = dill.load(fo)
        matA.set_hrd()
        matA.set_yld()
        matB.set_hrd()
        matB.set_yld()
    else:
        raise IOError, 'Unexpected case'
        # ## Should work on here to allow
        # ## both A and B materials are described using the
        # ## same constitutive model
        # matA = material
        # matB = material

    rad2deg  = 180./np.pi
    deg2rad  =   1./rad2deg

    stressA_off, dum1, dum2 = constructBC(
        epsAng  = th,
        f_yld   = matA.f_yld,
        verbose = False)

    ## put the stress on the locus
    matA.update_yld(stressA_off)
    np.set_printoptions(precision=3)
    print('stressA:'+('%7.3f'*6)%(
        matA.stress[0],matA.stress[1],matA.stress[2],
        matA.stress[3],matA.stress[4],matA.stress[5]))
    print('strainA:'+('%7.3f'*6)%(
        matA.dphi[0],matA.dphi[1],matA.dphi[2],
        matA.dphi[3],matA.dphi[4],matA.dphi[5]))
    alpha,rho = calcAlphaRho(matA.stress,matA.dphi)
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
    nind = max([len(matA.logfn),len(matB.logfn)])+3
    print('Iteration over the given psi angle')
    head = (
        '%8s'*9+  ## variables
        ('%'+'%is'%nind)*2+ ## aLogFN and bLogFN
        '%'+'%is'%(len(snapshot.logfn)+3))%(
            'epsRD','epsTD','psi0','psif','sigRD',
            'sigTD','sigA','T','cmpt','aLogFN','bLogFN','ssFN')
    head = '%s\n'%head
    logFile.write(head)
    t0   = time.time()

    ynew, absciss, xbb= onepath(
        matA=matA,matB=matB,
        psi0=psi0*deg2rad,f0=f0,
        T=absciss,snapshot=snapshot)

    matA.recordCurrentStat()
    matB.recordCurrentStat()

    dTime = time.time() - t0
    psif1 = xbb[0]

    cnt = (
        '%8.3f'*8+
        '%8i'+
        ('%'+'%is'%nind)*2+
        '%'+'%is'%(len(snapshot.logfn)+3))%(
        ynew[1],ynew[2],psi0,
        psif1*rad2deg,
        matA.stress[0],matA.stress[1],
        matA.sig, ## hardening (effective stress)
        absciss,dTime,matA.logfn,matB.logfn,snapshot.logfn)
    print(cnt)
    logFile.write(cnt+'\n')
    uet(dTime,'total time spent');print('')
    logFile.close()
    print('%s has been saved'%logFileName)
    return logFileName,dTime, matA, matB

def onepath(matA,matB,psi0,f0,T,snapshot):
    """
    Run under the given condition that is
    characterized by the passed arguments

    Arguments
    ---------
    matA
    matB
    psi0
    f0
    T
    snapshot
    """
    import os
    from mk.library.lib import rot_6d
    from mk.materials.func_hard_for import return_swift

    ## A stress state referred in band axes
    sx = rot_6d(matA.stress,-psi0)
    # ## strain hardening can be passed
    # independently as the was f_yld is passed.

    matA.update_hrd(0.) ## initialize hardening parameters

    ## initial_conditions
    ndim = 4
    b    = np.zeros(20)
    b[0] = psi0
    b[1] = f0
    b[2] = matA.sig
    b[3] = matA.m
    b[4] = matA.qq
    ## stress state ratio within the band from
    ## region A stress state
    b[5] = sx[5]/sx[0]

    xzero = np.array([1,1,0,0])

    ## Determine the initial states
    ## x[:3]: stress of region b referred in the band axes
    ## x[:3] = [s11, s22, s12] of the region B referred in the band axes

    xfinal, fb = new_raph_fld(
        ndim=ndim,ncase=1,
        xzero=xzero,b=b,
        matA=matA,
        matB=matB,
        verbose=False)

    ## fb is the first derivative of region B yield function
    ## that gives the 'directions' of strain rate according to the AFR

    ## Initial values
    tzero = xfinal[3] ## d\labmda
    yzero = np.zeros(5)
    ## fb: first derivative in region B
    yzero[3] = tzero*fb[0] ## B strain 'increment' along RD
    yzero[4] = tzero*fb[1] ## B strain 'increment' along TD

    ndds    = 5 ## dimension differential system
    dydx    = np.zeros(ndds)
    dydx[0] = 1.

    ## xbb = [psi0,s1,s2,s3,s4,s5,s6]
    xbb    = np.zeros(7)
    xbb[0] = psi0
    ## caution: xbb[1:] corresponds to the stress components
    ## sx: A stress state referred in band axes
    xbb[1] = sx[0]
    xbb[2] = sx[1]
    xbb[6] = sx[5]

    t0=time.time()
    ## integrate through monotonic loading
    ynew,absciss,xbb\
        = integrateMono(
            f0,
            tzero,
            yzero,
            ndds,
            dydx,
            xbb,
            matA,
            matB,
            snapshot,
            verbose=False)

    psif = xbb[0]
    print ('%8.3f'*5)%(ynew[0],ynew[1],ynew[2],ynew[3],ynew[4])
    uet(time.time()-t0,'Elapsed time in step by step integration')
    print '\nAfter integrateMono'
    print 'absciss:',absciss

    ## check the hardening curve?
    return ynew,absciss,xbb

def integrateMono(
        f0,
        tzero,
        yzero,
        ndds,
        dydx,
        xbb,
        matA,
        matB,
        snapshot,
        verbose):
    """
    Step by step integration

    f0     : initial inhomgeneity factor
    S      : stress state of region A
    tzero
    yzero  :
        y[1] accumulative strain RD
        y[2] accumulative strain TD
    ndds   : dimension differential system
    dydx   :
    xbb    : [psi0, s1, s2, s3, s4, s5, s6]
    matA
    matB
    snapshot: snapshot object to record state variables to study
    verbose: flag to be or not to be verbose

    Returns
    -------
    ynew
    absciss (T)
    xbb
    """
    import os
    # S       = matA.stress
    absciss = tzero ## xcoordinate in the intergration
    yold = np.zeros(ndds) ## y_old
    yold = yzero[:]

    ## integration values
    nbpas  = 200000
    freq   = 100      # frequency of output (probably not needed for this)

    ## --- delta t
    deltat = 1e-3
    # (incremental stepsize)
    #
    tmax   = 3.0      # (maximum effective strain upper limit (2.0?))

    k =-1
    t = tzero
    time_used_in_syst=0.
    totalTimeFunc=0.
    while(k<=nbpas and absciss<tmax and dydx[0]>=1e-1):
        """
        dydx[0] = d\lambda^A /  d\lambda^B

        Forming limit criterion: if dydx<0.1

        i.e., the instant when the equivalent strain rate
        of region a becomes far less than that of region B.

        that implies that dydx = d\lambda^A / d\lamda^B
        """
        k=k+1
        ## adjusting the incremental size size
        if dydx[0]<0.5:
            deltt = deltat/10.
            if dydx[0]<0.2:
                deltt = deltat/100.
        else:
            deltt = deltat*1.0
        t0 = time.time()

        ## find solution at current deformation increment
        dydx, ynew, xbb, totalTimeFunc = syst(
            deltt,
            t,
            f0,
            dydx,
            xbb,
            yold,
            matA,
            matB,
            snapshot,
            verbose,totalTimeFunc)

        ## record current status of the two regions
        ## might need to reduce the writing frequency
        if np.mod(k,100)==0: ## every 10 steps
            matA.recordCurrentStat()
            matB.recordCurrentStat()
            snapshot.takeshot(
                k=k,              #0
                deltt=deltt,      #1

                dydx0=dydx[0],    #2  ## partials of d(ynew[0])/dx[0]
                dydx1=dydx[1],    #3  ## partials of d(ynew[0])/dx[1]
                dydx2=dydx[2],    #4  ## partials of
                dydx3=dydx[3],    #5
                dydx4=dydx[4],    #6

                xbb0 = xbb[0],    #7
                xbb1 = xbb[1],    #8
                xbb2 = xbb[2],    #9
                xbb3 = xbb[3],    #10
                xbb4 = xbb[4],    #11
                xbb5 = xbb[5],    #12
                xbb6 = xbb[6],    #13

                ynew0=ynew[0],    #14
                ynew1=ynew[1],    #15 ERD
                ynew2=ynew[2],    #16 ETD
                ynew3=ynew[3],    #17 ERD
                ynew4=ynew[4]     #18 ETD
            )

            snapshot.linebreak()

        time_used_in_syst = time_used_in_syst + (time.time()-t0)
        k1      = deltt * dydx ## Y increments
        ynew    = yold + k1
        t       = t +deltt ## accmulated equivalent strain
        absciss = t*1.
        yold[::]=ynew[::]

    uet(time_used_in_syst,'Total time used for iteration in syst')
    print
    uet(totalTimeFunc,'Total time elapsed in func')
    print

    return ynew,absciss,xbb

## differential system
def syst(
        deltt,
        t,
        f0,
        dydx,
        xbb,
        y,
        matA,
        matB,
        snapshot,
        verbose,
        totalTimeFunc=0.):
    """
    Arguments
    ---------
    deltt   : equivalent strain increment of region B
    t       : axis along the intergration occurs
    f0      : initial inhomogeneity
    dydx    :
    xbb     : [psi0,s1,s2,s3,s4,s5,s6]
              -- psi0 and stress state of region b
    y       : yold defined in integrateMono
        y[1]: accumulative strain RD
        y[2]: accumulative strain RD
    matA
    matB
    snapshot : (if not None, activated)
    verbose

    Returns
    -------
    dydx
    yold
    siga     strain hardening flow stress, sig = hard(E)
    totalTimeFunc
    """
    import os
    """
    xzero is the initial guesses
    xzero[0] = equivalent strain increment of region A
    """
    xzero    = np.zeros(4)
    xzero[0] = dydx[0]*deltt ## use \Delta\lambda^A * \Delta T as a guess
    xzero[1] = xbb[1] ## s1
    xzero[2] = xbb[2] ## s2
    xzero[3] = xbb[6] ## s6 of region B (stress referred in... )

    ndim  = 4
    ncase = 2

    bn    = np.zeros(20)
    bn[1] = f0
    bn[8] = deltt
    bn[9] = xbb[0] ## psi0

    xfinal, fa, fb, bn, totalTimeFunc\
        = new_raph_fld(
            T=t,
            ndim=ndim,
            ncase=2,
            xzero=xzero,
            y=y,
            b=bn,
            matA=matA,
            matB=matB,
            verbose=verbose,totalTimeFunc=totalTimeFunc)

    xbb[0] = bn[9]+bn[10]  ## psi^n + \delta psi
    xbb[1:3] = xfinal[1:3]
    # xbb[1] = xfinal[1]
    # xbb[2] = xfinal[2]
    xbb[5] = xfinal[3]

    dydx[0] = xfinal[0]/deltt ## delta lambda^A / lambda^B
    dydx[1] = fa[0]*dydx[0]
    dydx[2] = fa[1]*dydx[0]
    dydx[3] = fb[0]
    dydx[4] = fb[1]

    return dydx, y, xbb, totalTimeFunc

def new_raph_fld(
        T=None,ndim=None,ncase=None,
        xzero=None,y=None,b=None,#f_hard=None,f_yld=None,
        matA=None,matB=None,
        verbose=True,totalTimeFunc=0):
    """
    Find numerical solution using Newton-Raphson method.
    The relevant jacobian and objective functions are defined
    in func_fld.py.

    Iterative determination is done
    using the subroutine gauss in for.f

    Arguments
    ---------
    T
    ndim
    ncase
    xzero
    y
    b
    matA
    matB
    verbose=True
    totalTimeFunc

    Return
    ------
    xn1, fb (case 1)
    xn1, fa, fb, b
    """
    import os, time
    from yf_for import gauss,norme
    from func_fld import func_fld1, func_fld2

    residu  = 1.0
    xn      = xzero[::]
    it      = 0
    # itmax   = 200
    itmax   = 20
    # eps     = 1e-10
    eps = 1e-4

    # totalTimeFunc = 0.
    dt = 0.
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
            F, J, fb        = func_fld1(
                ndim,b,xn,matA,matB,verbose)
            dt = time.time() - t0
        if ncase==2:
            if verbose:
                print '-'*40
                print '%i ITERATION over func_fld2 in NR'%it
            F, J, fa, fb, b\
                = func_fld2(
                    ndim,
                    T,
                    b,
                    xn,
                    y,
                    matA,
                    matB,
                    verbose)
            dt = time.time() - t0
        totalTimeFunc = totalTimeFunc + dt ## to estimate compute performance
        F,J,res = gauss(ndim=ndim,a=J,b=F) ## f2py module
        xn1=xn-res ## x^(n+1)
        residu = norme(ndim,res)
        xn=xn1[::] ## x^(n)

    # uet(totalTimeFunc,'Time elapsed in new_raph_fld')

    # if no convergence
    if it>=itmax:
        print 'could not converge'
        return

    if ncase==1: return xn1,fb
    if ncase==2: return xn1,fa,fb,b,totalTimeFunc

## command line usage (may be called in mk_run for multi-threaded run)
"""
test case for command line usage:

$ python main.py --fn /tmp/dummy-log-file-name -f 0.995 -p 0 -t 0 --mat 0
"""
if __name__=='__main__':
    from MP import progress_bar
    import argparse, mk.materials.materials, mk.materials.constitutive,dill,pickle
    from mk.library.lib import gen_tempfile
    import mk.materials.func_hard_for
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
    parser.add_argument(
        '--mat', type=int, default=-1,
        help='Material card in materials.py (e.g., 0: IsoMat)')
    parser.add_argument(
        '--fnyld', type=str, default=None,
        help='Yield function (pickled) binary file name')
    parser.add_argument(
        '--fnhrd', type=str, default=None,
        help='Strain-hardening function (pickled) binary file name')
    #-------------------------------------------------------
    args = parser.parse_args()

    print 'args.mat:', args.mat
    print 'args.fnyld:', args.fnyld
    print 'args.fnhrd:', args.fnhrd

    ## Determine material cards  --
    if type(args.fnyld).__name__=='str' \
       and type(args.fnhrd).__name__=='str':

        with open(args.fnyld,'rb') as fo:
            yf_label=pickle.load(fo)
            p_yld = pickle.load(fo)
        with open(args.fnhrd,'rb') as fo:
            fhrd_type = dill.load(fo)
            p_hrd = dill.load(fo)

        matClass = mk.materials.constitutive.Constitutive(
            f_yld=None,f_hrd=None,
            params_yld=p_yld,label_yld=yf_label,
            params_hrd=p_hrd,label_hrd=fhrd_type)

        fn = gen_tempfile(prefix='mkmat',ext='dll')
        with open(fn,'wb') as fo:
            dill.dump(matClass,fo)
        mat = fn ## save file name to mat

    elif args.mat==-1:
        raise IOError, 'Unexpected case 1:   %s %s'%(
            type(args.fnyld).__name__,
            type(args.fnhrd).__name__)
    elif args.mat!=-1:
        print '-'*50
        print args.fnyld
        print args.fnhrd
        print '-'*50
        mat = mk.materials.materials.library(args.mat)
    else:
        raise IOError, 'Unexpected case 2:   %s %s %s'%(
            type(args.fnyld).__name__,
            type(args.fnhrd).__name__,
            type(args.mat).__name__)
    ## end of material card determination --

    if type(mat).__name__=='NoneType':
        raise IOError, 'None returned from the library'
    print 'mat:', mat
    print 'logfilename:',args.fn
    main(f0=args.f,psi0=args.p,th=args.t,logFileName=args.fn,material=mat)
    pass
