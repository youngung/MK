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
