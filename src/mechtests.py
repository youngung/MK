"""
various mechanical tests using
"""
import numpy as np
from lib import rot_6d

def inplaneTension(fYLD):
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
        ysMat, Phi, dPhiMat, d2PhiMat = fYLD(s=sMat)
        ysLab = rot_6d(ysMat,  -psis[i])
        deLab = rot_6d(dPhiMat,-psis[i])
        phis.append(Phi)
        e1,e2,e3 = deLab[0],deLab[1],- deLab[0]-deLab[1]
        rvs.append(e2/e3)

    return psis, rvs, phis
