"""
Proportional paths considered in marciniak-kuczynski problems
"""

from for_lib import vm
import numpy as np


def DRD():
    """
    Drawing condition RD
    """
    f_yld  = vm ## material card...

    ang1   = 0.
    ang2   = 25.
    anginc = 1.
    angs   = [ang1,ang2,anginc]

    ## reference strain vector (pth) of rho
    pth    = np.zeros(4)
    pth[0] = 1.
    pth[1] = -0.6
    pth[2] = 1.
    pth[3] = 0.
    npt    = 6

    s1=np.zeros(6); s2=np.zeros(6)
    s1[0]=0; s1[1]=-1; s2[0]=1; s2[1]=1

    return angs,npt,pth,f_yld,s1,s2

def PSRD():
    """
    Near plane-strain RD
    """
    f_yld = vm ## material card...

    ang1   = 0.
    ang2   = 25.
    anginc = 5.
    angs   = [ang1,ang2,anginc]

    ## reference strain vector (pth) of rho
    pth    = np.zeros(4)
    pth[0] = 1.
    pth[1] = 0.025
    pth[2] = 1.
    pth[3] = 0.25
    npt    = 6

    s1=np.zeros(6); s2=np.zeros(6)
    s1[0]=1; s1[1]=0; s2[0]=0; s2[1]=1

    return angs,npt,pth,f_yld,s1,s2
