### yf version that relies on fortran version of yield functions
from for_lib import vm, hqe
from lib import c6p, s62c, rot_tensor_r, c2s6, rotPrincOrig ## cauchy stress to principal stresses
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

    s33    = rotPrincOrig(rotMatrix,sPrinc)
    snew   = c2s6(s33)
    dphi33 = rotPrincOrig(rotMatrix,dphi)
    dphi   = c2s6(dphi33)

    # d2phi = np.dot(rotMatrix,dphi)  - need to define d2phi in 'for.f' first.
    return snew,phi,dphi,d2phi

def test2():
    from lib import rot
    import numpy as np
    psis = np.linspace(-np.pi,+np.pi)
    s=[1.,0.,0.,0.,0.,0.]
    s33 = s62c(s)

    yl=[]
    for i in xrange(len(psis)):
        psi = psis[i]
        # r = rot(psi)
        # s33new = rot_tensor_r(s33,r) ## rotated material axes
        # s6new = c2s6(s33new)

        s6=np.zeros(6)
        s6[0] = (np.cos(psi))**2
        s6[1] = (np.sin(psi))**2
        s6[2] = 0.0
        s6[3] = 0.0
        s6[4] = 0.0
        s6[5] = np.cos(psi)*np.sin(psi)

        snew, phi, dphi, d2phi = HillQuad(s=s6,r0=2,r90=2.3)
        strs = snew[0]*(np.cos(psi))**2 + snew[1]*(np.sin(psi))**2 + 2*snew[5]*np.sin(psi)*np.cos(psi)

        # RVAL = (E(1)*(DSIN(AL))**2 + E(2)*(DCOS(AL))**2 - 2*E(6)*DSIN(AL)*DCOS(AL))/(-E(1)-E(2))

        yl.append(strs)


    import matplotlib.pyplot as plt
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(psis*180./np.pi,yl)
    fn = 'yf2_test2.pdf'
    fig.savefig(fn)
    print '%s has been saved'%fn


def test1():
    s1  = [1.,0.,0.,0.,0.,0.] ## uniaxial s11
    s2  = [0.,1.,0.,0.,0.,0.] ## uniaxial s22
    s3  = [0.,0.,0.,0.,0.,1.] ## uniaxial s12
    s4  = [1.,1.,0.,0.,0.,0.] ## uniaxial balanced biaxial
    s5  = [1.,-1.,0.,0.,0.,0.] ## pure shear1
    s6  = [-1.,1.,0.,0.,0.,0.] ## pure shear2

    ss = [s1,s2,s3,s4,s5,s6]
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

    ax.set_aspect('equal')

    fn ='yf2_test.pdf'
    fig.savefig(fn)
    print '%s has been saved'%fn

if __name__=='__main__':
    # test1()
    test2()
