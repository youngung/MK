### yf version that relies on fortran version of yield functions
from for_lib import vm, hqe
from lib import c6p, rot_tensor_r, c2s6 ## cauchy stress to principal stresses
import numpy as np

def VonMises(s):
    snew,phi,dphi,d2phi = vm(s)
    return snew,phi,dphi,d2phi

def HillQuad(s=None,r0=2.0,r90=2.3):
    """
    Arguments
    ---------
    s    : 6d stress
    r0   : r-value along 0
    r90  : r-value along 90
    """
    princStress,rotMatrix = c6p(s)
    sPrinc,phi,dphi,d2phi = hqe(s=princStress,r0=r0,r90=r90) ## principal values

    s1 = rotMatrix[:,0] * sPrinc[0]
    s2 = rotMatrix[:,1] * sPrinc[1]
    s3 = rotMatrix[:,2] * sPrinc[2]
    stressPrinc = np.zeros((3,3))
    stressPrinc[:,0] = s1[::]
    stressPrinc[:,1] = s2[::]
    stressPrinc[:,2] = s3[::]

    # stressOrig = np.zeros((3,3))

    ## snew, dphi, d2phi all are referred to the principal space.
    ## should rotate back to the original space... and how?
    ## rotMatrix rotates from te principal space to the original

    ## Below has not been thoroughly checked yet
    # snew = np.dot(rotMatrix,snew)
    s33  = rot_tensor_r(stressPrinc,rotMatrix)
    snew = c2s6(s33)

    # dphi = np.dot(rotMatrix,dphi)


    # d2phi = np.dot(rotMatrix,dphi)  - need to define d2phi in 'for.f' first.
    return snew,phi,dphi,d2phi

def test1():
    s1  = [1.,0.,0.,0.,0.,0.] ## uniaxial s11
    s2  = [0.,1.,0.,0.,0.,0.] ## uniaxial s22
    s3  = [0.,0.,0.,0.,0.,1.] ## uniaxial s12
    s4  = [1.,1.,0.,0.,0.,0.] ## uniaxial balanced biaxial

    ss = [s1,s2,s3,s4]
    r0 =1.; r90=1.
    xs=[];ys=[];zs=[]
    for i in xrange(len(ss)):
        snew,phi,dphi,d2phi = HillQuad(s=ss[i],r0=r0,r90=1.)
        s11,s22,s12 = snew[0],snew[1],snew[5]
        xs.append(s11)
        ys.append(s22)
        zs.append(s12)

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax  = fig.add_subplot(111,projection='3d')
    ax.scatter(xs,ys,zs)
    ax.set_xlabel(r'$\sigma_{11}$')
    ax.set_ylabel(r'$\sigma_{22}$')
    ax.set_zlabel(r'$\sigma_{12}$')

    fn ='yf2_test.pdf'
    fig.savefig(fn)
    print '%s has been saved'%fn

if __name__=='__main__':
    test1()
