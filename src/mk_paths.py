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
from for_lib import vm
import numpy as np

def DRD():
    """
    Drawing condition RD
    """
    f_yld  = vm ## material card...
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

    return angs,npt,pth,f_yld,stressRight,stressLeft

def PSRD():
    """
    Near plane-strain RD
    """
    f_yld  = vm ## material card... should replace this.
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

    return angs,npt,pth,f_yld,stressRight,stressLeft

def BBRD():
    """
    Near biaxial RD
    """
    f_yld  = vm ## material card...
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


    return angs,npt,pth,f_yld,stressRight,stressLeft

def BBTD():
    """
    Near biaxial TD
    """
    f_yld  = vm ## material card...
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

    return angs,npt,pth,f_yld,stressRight,stressLeft

def PSTD():
    """
    Near plane-strain TD
    """
    f_yld  = vm ## material card...
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

    return angs,npt,pth,f_yld,stressRight,stressLeft

def DTD():
    """
    Drawing condition TD
    """
    f_yld  = vm ## material card...
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

    return angs,npt,pth,f_yld,stressRight,stressLeft

def returnPaths():
    return DRD, PSRD, BBRD, BBTD, PSTD, DTD

def testAllPaths():
    paths = returnPaths()
    for i in xrange(len(paths)):
        testEachPath(paths[i])

def testEachPath(funcPath=None):
    from mk_lib import findStressOnYS
    if type(funcPath).__name__=='NoneType':
        angs,npt,pth,f_yld,stressR, stressL = DRD()
    elif type(funcPath).__name__=='function':
        angs,npt,pth,f_yld,stressR, stressL = funcPath()
    else:
        raise IOError, 'Error!'

    ang1,ang2,anginc = angs
    for npth in xrange(npt):
        SR = stressR[::]
        SL = stressL[::]
        pthr = pth[0] + npth*(pth[2]-pth[0])/npt
        ptht = pth[1] + npth*(pth[3]-pth[1])/npt

        strainVector = [pthr,ptht]

        s,dphi = findStressOnYS(
            f_yld,SR.copy(),SL.copy(),
            pth=strainVector,
            verbose=False)

if __name__=='__main__':
    ## test.
    test()
