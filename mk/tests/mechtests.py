"""
Various mechanical tests...

<inplaneTension>
  - uniaxial tension tests to reveal the anisotropy in terms
    of r-values and uniaxial yield stresses at various angles
<locus>
  - Conduct in-plane biaxial stress tests to reavel yield
    locus in plane-stress space.
"""
import numpy as np
from mk.library.lib import rot_6d
import mk.yieldFunction.yf2

def inplaneTension(fYLD,iopt=1,**kwargs):
    """
    inplane tesion tests provides
    R-value and yield stress variations along difference
    direction from the rolling direction of sheet metals

    fYLD: yield function (a function of stress state only)

    Returns
    -------
    psis (rotation angles)
    rvs  (r-values)
    phis (yield stresses)
    """
    psis = np.linspace(0,+np.pi/2.,100)

    ## stress state for inplane tension in the lab axes
    ## is of course uniaxial stress state:
    sLab=np.zeros(6)
    sLab[0]=1.

    phis=[];rvs=[]
    ## rotate this stress state to 'material' axes for
    ## each psi angles and collect results
    for i in xrange(len(psis)):
        sMat = rot_6d(sLab, psis[i])
        ysMat, Phi, dPhiMat, d2PhiMat = fYLD(s=sMat,**kwargs)
        ysLab = rot_6d(ysMat,  -psis[i])
        deLab = rot_6d(dPhiMat,-psis[i])
        if iopt==0:
            phis.append(Phi)
        elif iopt==1:
            phis.append(ysLab[0])
        e1,e2,e3 = deLab[0],deLab[1],- deLab[0]-deLab[1]
        rvs.append(e2/e3)
    return psis, rvs, phis

def locus(func,nth=100):
    """
    in-plane biaxial locus

    Arguments
    ---------
    func       (yield function)
    nth = 100
    """
    import numpy as np
    pi = np.pi
    cos=np.cos
    sin=np.sin
    th=np.linspace(-pi,pi,nth)
    x=cos(th); y=sin(th)
    z=np.zeros(len(th))
    s=np.array([x,y,z,z,z,z]).T
    X=[]; Y=[]
    for i in xrange(len(s)):
        rst = func(s[i])
        ys = rst[0]
        X.append(ys[0])
        Y.append(ys[1])
    return np.array(X),np.array(Y)

def locus_shear(func,nth=100):
    """
    S11/S12 locus

    Arguments
    ---------
    func       (yield function)
    nth = 100
    """
    pi = np.pi
    cos=np.cos
    sin=np.sin
    th=np.linspace(-pi,pi,nth)
    x=cos(th); y=sin(th)
    z=np.zeros(len(th))
    s=np.array([x,z,z,z,z,y]).T
    X=[]; Y=[]
    for i in xrange(len(s)):
        rst = func(s[i])
        ys = rst[0]
        X.append(ys[0])
        Y.append(ys[5])
    return np.array(X),np.array(Y)

def defineObjf(th,yfunc,s12):
    def obj(rad):
        s11,s22=rad*np.cos(th),rad*np.sin(th)
        s=[s11,s22,0,0,0,s12]
        snew,phi,dphi,d2phi=yfunc(s)
        snew = snew/phi
        return abs(snew[5] - s12)
    return obj

def surface(func):
    """"
    3D surface of (s11,s22,s12)
    Return (S11,S22) contours with a fixed level of s12.

    Arguments
    --------
    func : yield fuction
    """
    import time
    from scipy.optimize import minimize

    t0=time.time()
    ## find s12 yield stress
    snew, phi, dphi, d2phi = func([0,0,0,0,0,1.])
    y_s12 = snew[-1]

    nlev_s12 = 20
    s12Cnts = np.linspace(0,y_s12, nlev_s12)[:-1] ## exlude the last level

    nths = 100
    ths = np.linspace(-np.pi,np.pi,100) ## azimuth
    yieldSurface = np.zeros((nlev_s12,nths,3)) ## altitude

    rad = 1.
    for iz in xrange(nlev_s12-1):
        s12 = s12Cnts[iz] ## level, s12 should be met by this value.
        for i in xrange(len(ths)):
            objf = defineObjf(ths[i],func,s12Cnts[iz])
            res = minimize(
                fun=objf,x0=rad,
                tol=1e-20,method='BFGS',jac=False)
            rad = res.x
            s11,s22 = rad*np.cos(ths[i]), rad*np.sin(ths[i])
            snew,phi,dphi,d2phi = func([s11,s22,0,0,0,s12])
            yieldSurface[iz,i,:]=snew[0],snew[1],snew[5]
    print 'Elapsed time in mk.tests.mechtests.surface:', time.time()-t0
    return yieldSurface, y_s12

def testSurface(
        yfunc =None): #mk.yieldFunction.yf2.wrapHill48R([1.5,2.4,2.9])):
    """
    Calculate the 3D space yield surface in (s11,s22,s12)

    Arguments
    ---------
    yfunc
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import MP.lib.mpl_lib
    import matplotlib as mpl

    ## Test Hill48
    # yfunc = mk.yieldFunction.yf2.wrapYLD(r=[1.5,2.4,2.9,1.],y=[1,1,1,1.3])
    if type(yfunc).__name__=='int':
        raise IOError

    yieldSurface, y_s12 = surface(yfunc)
    cmap = MP.lib.mpl_lib.color_map(mn=0,mx=y_s12,cmap='viridis')

    fig1 = plt.figure()
    ax1  = fig1.add_subplot(111)
    fig2 = plt.figure(figsize=(6,5))
    ax2  = fig2.add_subplot(111,projection='3d')
    ax2.view_init(elev=25., azim=45.)

    for i in xrange(len(yieldSurface)):
        x=yieldSurface[i,:,0]
        y=yieldSurface[i,:,1]
        z=yieldSurface[i,:,2]

        c = cmap(z[0])
        ax1.plot(x,y,color=c)
        ax2.plot(x,y,z,color=c)
        pass

    ax1.set_aspect('equal'); ax2.set_aspect('equal')
    for ax in [ax1,ax2]:
        ax.set_xlabel(r'$s_{11}$')
        ax.set_ylabel(r'$s_{22}$')
    ax2.set_zlabel(r'$s_{12}$')
    ax2.set_xticks([-1,0,1])
    ax2.set_yticks([-1,0,1])
    ax2.set_zticks([0,0.5,1])

    fig1.tight_layout(); fig2.tight_layout()
    fig1.savefig('yieldsurface_2d.pdf',bbox_to_inches='tight')
    fig2.savefig('yieldsurface_3d.pdf',bbox_to_inches='tight')

    return yieldSurface, y_s12

if __name__=='__main__':
    testSurface()
