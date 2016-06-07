"""
Proportional paths considered in Marciniak-Kuczynski problems

The strain spaces are divided into DRD, PSRD, BBRD, BBTD, PSTD, and DTD.

Depending on the region, one might considered different 'guesses'
for the initial conditions used for numerical solver used in FLD
simulations. Each of functions listed in below gives
1. angles of psi
2. the values of rhos to be probed (in terms of strain rate vector)
3. left and rigt stresses that bounds the stress regions that
   would is likely to result from the given 'rhos'.


---------------------
DRD : drawing in RD
PSRD: plane-strain RD
BBRD: biaxial RD
BBTD: biaxial TD
PSTD: plane-strain TD
DTD : drawing in TD
"""
from yf_for import vm
import numpy as np

def DRD():
    """
    Drawing condition RD
    """
    ang1   =  0
    ang2   = 25
    anginc =  1
    angs   = [ang1,ang2,anginc]

    ## reference strain vector (pth) of rho
    ## from -0.6 to 0 with 6 number of probings
    pth    = np.zeros(4)
    pth[0] = 1.000
    pth[1] =-0.600
    pth[2] = 1.000
    pth[3] = 0.000
    npt    = 6

    ## Actual stress state should be bounded by
    ## below 'left' and 'right' stress points
    ## to be used to numerically find the actual stress
    ## corresponding to the strain vector given by pth

    stressRight=np.zeros(6)
    stressLeft=np.zeros(6)
    stressRight[0]= 0
    stressRight[1]=-1
    stressLeft[0] = 1
    stressLeft[1] = 1

    return angs,npt,pth,stressRight,stressLeft

def PSRD():
    """
    Near plane-strain RD
    """
    ang1   =  0
    ang2   = 25
    anginc =  5
    angs   = [ang1,ang2,anginc]

    ## reference strain vector (pth) of rho
    ## from 0.025 to 0.25
    pth    = np.zeros(4)
    pth[0] = 1.000
    pth[1] = 0.025
    pth[2] = 1.000
    pth[3] = 0.250
    npt    = 9

    stressRight=np.zeros(6)
    stressLeft=np.zeros(6)
    stressRight[0]=1
    stressRight[1]=0
    stressLeft[0]=0
    stressLeft[1]=1

    return angs,npt,pth,stressRight,stressLeft

def BBRD():
    """
    Near biaxial RD
    """
    ang1   =  0
    ang2   = 90
    anginc =  5
    angs   = [ang1,ang2,anginc]

    ## reference strain vector (pth) of rho
    pth    = np.zeros(4)
    pth[0] = 1.000
    pth[1] = 0.300
    pth[2] = 1.000
    pth[3] = 1.000
    npt    = 14

    stressRight=np.zeros(6)
    stressLeft=np.zeros(6)
    stressRight[0]=1
    stressRight[1]=0
    stressLeft[0]=0
    stressLeft[1]=1


    return angs,npt,pth,stressRight,stressLeft

def BBTD():
    """
    Near biaxial TD
    """
    ang1   = 90
    ang2   =  0
    anginc = -5
    angs   = [ang1,ang2,anginc]

    ## reference strain vector (pth) of rho
    pth    = np.zeros(4)
    pth[0] = 0.950
    pth[1] = 1.000
    pth[2] = 0.250
    pth[3] = 1.000
    npt    = 14

    stressRight=np.zeros(6)
    stressLeft=np.zeros(6)
    stressRight[0]=1
    stressRight[1]=0
    stressLeft[0]=0
    stressLeft[1]=1

    return angs,npt,pth,stressRight,stressLeft

def PSTD():
    """
    Near plane-strain TD
    """
    ang1   = 90
    ang2   = 60
    anginc = -5
    angs   = [ang1,ang2,anginc]

    ## reference strain vector (pth) of rho
    pth    = np.zeros(4)
    pth[0] = 0.225
    pth[1] = 1.000
    pth[2] = 0.000
    pth[3] = 1.000
    npt    = 9

    stressRight=np.zeros(6)
    stressLeft=np.zeros(6)
    stressRight[0]=1
    stressRight[1]=0
    stressLeft[0]=0
    stressLeft[1]=1

    return angs,npt,pth,stressRight,stressLeft

def DTD():
    """
    Drawing condition TD
    """
    ang1   = 90
    ang2   = 65
    anginc = -1
    angs   = [ang1,ang2,anginc]

    ## reference strain vector (pth) of rho
    pth    = np.zeros(4)
    pth[0] =-0.100
    pth[1] = 1.000
    pth[2] =-0.600
    pth[3] = 1.000
    npt    = 5

    stressLeft=np.zeros(6)
    stressRight=np.zeros(6)
    stressRight[0]=1
    stressRight[1]=1
    stressLeft[0]=-1
    stressLeft[1]=0

    return angs,npt,pth,stressRight,stressLeft

def returnPaths():
    return DRD, PSRD, BBRD, BBTD, PSTD, DTD

def findCorrectPsi(epsAng):
    """
    Recommend the most required (relevnt) psi anges
    for the given strain paths <epsAng>

    Arguments
    ---------
    epsAng in [degree]
    """
    path = findCorrectPath(epsAng)[0]

    print 'path name:',path.__name__

    rst = path()
    ang1,ang2,anginc = rst[0]

    ntotAngle = ((ang2-ang1)/anginc)+1
    rad2deg   = 180./np.pi
    deg2rad   =   1./rad2deg
    psi0s     = np.linspace(ang1,ang2,ntotAngle)*deg2rad
    return psi0s

def findCorrectPath(epsAng):
    """
    Argument
    --------
    epsAng [in degree]
    """
    import numpy as np
    pth = np.ones(2)
    e11 = np.cos(epsAng*np.pi/180.)
    e22 = np.sin(epsAng*np.pi/180.)

    if e11>e22:
        rho = e22/e11
    else:
        rho = e11/e22

    if epsAng<45:
        pth[1]=rho
    else:
        pth[0]=rho

    ## RD // Axis1
    if -45.<epsAng<0:
        return DRD, pth
    elif 0<=epsAng<30:
        return PSRD, pth
    elif 30<=epsAng<45:
        return BBRD, pth


    ## TD // Axis1
    elif 45<=epsAng<60:
        return BBTD, pth
    elif 60<=epsAng<90:
        return PSTD, pth
    elif 90<=epsAng<145:
        return DTD, pth
    else:
        raise IOError,'Could not find the corrent range'

def testAllPaths():
    paths = returnPaths()
    for i in xrange(len(paths)):
        testEachPath(paths[i])
        pass
    pass

def testEachPath(funcPath=None):
    from mk_lib import findStressOnYS
    if type(funcPath).__name__=='NoneType':
        angs,npt,pth,f_yld,stressR, stressL = DRD()
    elif type(funcPath).__name__=='function':
        angs,npt,pth,f_yld,stressR, stressL = funcPath()
    else:
        raise IOError, 'Error!'

    ang1,ang2,anginc = angs
    SR = stressR[::]
    SL = stressL[::]
    for npth in xrange(npt):
        pthr = pth[0] + npth*(pth[2]-pth[0])/npt
        ptht = pth[1] + npth*(pth[3]-pth[1])/npt

        strainVector = [pthr,ptht]

        s,dphi = findStressOnYS(
            f_yld,SR.copy(),SL.copy(),
            pth=strainVector,
            verbose=False)
        pass
    pass

def constructBC(epsAng,f_yld,verbose=False):
    """
    Arguments
    ---------
    epsAng
    f_yld
    verbose
    """
    stress, stressLowerAngle, stressUpperAngle = \
        calcStressWindow(epsAng,f_yld,verbose)
    stressLower = th2s6(stressLowerAngle)
    stressUpper = th2s6(stressUpperAngle)
    if verbose:
        print 'upper',stressUpperAngle*180./np.pi,stressUpper
        print 'lower',stressLowerAngle*180./np.pi,stressLower

    return stress, stressLower, stressUpper

def th2s6(th):
    from numpy import cos, sin
    x=cos(th)
    y=sin(th)
    return np.array([x,y,0,0,0,0])

def th2th(th):
    """
    Restrain th within (-pi,pi) range.
    """
    x,y = np.cos(th),np.sin(th)
    return np.arctan2(y,x)

def calcStressWindow(theta=0,f_yld=None,verbose=False):
    """
    Given the rho, provide the stress that most probably
    embraces the corresponding stress state in yield locus.

    Argument
    --------
    theta: the angle by (eyy,exx), degree
    f_yld: yield function (if not(None),
            use it as a the actual yield function
            if None, von Mises will be used as a default

    1. Calculates the corresponding von Mises stress bound first.
    2. Return an extended range to make sure that the stress is embraced within

    Returns
    -------
    s, stressLowerAngle, stressUpperAngle  (all in radian)
    """
    from numpy import cos, sin, tan, pi
    d2r = pi/180.
    r2d = 1./d2r
    th = theta * d2r
    e1 = cos(th); e2 = sin(th)

    ## vm
    s=np.array([1,1,0,0,0,0],dtype='float') ## initial guess

    tol  = 1e-11
    diff = 1.
    dx   = 1e-13
    maxIter = 1000

    th_sig = np.arctan2(s[1],s[0])
    th_sig_ = th2th(th_sig)
    it = 0
    if verbose:
        print ('%4s %8s %8s %8s %11s %11s')%(
            'it','th_sig','th_eps','th','diff','jac')
    while abs(diff)>tol:
        it = it + 1
        if it>maxIter:
            raise IOError, 'could not find the stress ...'

        x0     = th_sig
        x1     = th_sig+dx
        F1     = objf(x1,f_yld)-th
        F0     = objf(x0,f_yld)-th
        jacob  = (F1-F0)/(x1-x0)
        th_sig = th_sig - F1 / jacob
        th_eps = objf(th_sig,f_yld)
        diff   = th_eps  - th

        th_sig_ = th2th(th_sig)
        th_eps_ = th2th(th_eps)
        if verbose: print '%4i %8.2f %8.2f %8.2f %11.3e %11.3e'%(
            it,th_sig_*r2d,th_eps_*r2d,th*r2d,diff*r2d,jacob)

    ## the final stress estimate:
    s[0] = cos(th_sig_)
    s[1] = sin(th_sig_)

    window = 40*d2r
    return s, th_sig_-window,  th_sig_+window

def objf(th=None,f_yld=None):
    """
    Given stress th, return strain rate angle (-pi,pi)
    th = arc2tan(s2,s1)
    """
    from numpy import cos, sin, tan, pi

    if type(f_yld).__name__=='NoneType':
        f_yld = vm
    # f2py modules are named 'fortran'
    elif type(f_yld).__name__=='fortran':
        pass
    elif type(f_yld).__name__=='function': ## just in case
        pass
    else:
        raise IOError, 'unexpected choice of f_yld'

    s1=cos(th); s2=sin(th)
    s=np.array([s1,s2,0,0,0,0])
    s,phi,dphi,d2phi = f_yld(s)
    th_eps = np.arctan2(dphi[1],dphi[0])
    return th_eps

if __name__=='__main__':
    ## test.
    testAllPaths() ## testing all pre-defined path functions

    ## test the calcStressWindow
    epsAngles = np.linspace(-45,135,30)

    for i in xrange(len(epsAngles)):
        th = epsAngles[i]
        calcStressWindow(theta=th,verbose=True)
