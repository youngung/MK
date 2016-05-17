### yf version that relies on fortran version of yield functions
import matplotlib as mpl
mpl.use('Agg') ## In case X-window is not available.
from for_lib import vm, hqe, hill48
from lib import c6p, s62c, rot_6d,rot_tensor_r, c2s6, rotPrincOrig ## cauchy stress to principal stresses
import numpy as np

def VonMises(s):
    snew,phi,dphi,d2phi = vm(s)
    return snew,phi,dphi,d2phi

def wrapHill48(r0,r90):
    import tuneH48
    Hill48params = tuneH48.main(r0=r0,r90=r90)
    f,g,h,n = Hill48params
    def func(s):
        return Hill48(s,f,g,h,n)
    return func

def Hill48(s,f,g,h,n):
    """
    Arguments
    ---------
    s    :6d-stress
    """
    snew,phi,dphi,d2phi = hill48(s,f,g,h,n)
    return snew, phi, dphi, d2phi

def wrapHillQuad(r0,r90):
    def func(s):
        return HillQuad(s,r0,r90)
    return func

def HillQuad(s=None,r0=2.0,r90=2.3):
    """
    Arguments
    ---------
    s    : 6d stress   (it should be referenced by material axes) axis1//RD, axi2//TD
    r0   : r-value along 0
    r90  : r-value along 90
    """
    ## incomplete here - axis1 // rd, it should be such that axis 2 // td
    ## This, in the current form, may not enforce axis 1 to e aligned with RD
    ##
    srd = s[0]; std = s[1]
    if srd>=std: pdir = 'rd'
    elif srd<std: pdir='td'

    princStress,rotMatrix = c6p(s)
    sPrinc,phi,dphi,d2phi = hqe(s=princStress,r0=r0,r90=r90) ## principal values

    s33    = rotPrincOrig(rotMatrix,sPrinc)
    snew   = c2s6(s33)
    dphi33 = rotPrincOrig(rotMatrix,dphi)
    dphi   = c2s6(dphi33)

    # d2phi = np.dot(rotMatrix,dphi)  - need to define d2phi in 'for.f' first.
    return snew,phi,dphi,d2phi

def test3():
    """
    """
    import vm_check
    vm_check.main()

def test2(r0=2.1,r90=2.7):
    """
    uniaxial tension tests

    Arguments
    ---------
    r0  =2.1
    r90 =2.7
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from lib import rot
    import numpy as np

    h48 = wrapHill48(r0=r0,r90=r90)

    psis = np.linspace(0,+np.pi/2.,100)
    ## uniaxial tension stress state in the laboratory axes
    sUniaxial=np.zeros(6)
    sUniaxial[0]=1.
    phis=[];x=[];y=[];z=[];rvs=[]

    print '%6s %6s %6s %6s %6s %6s %6s %6s'%(
        'th','phi','s11','s22','s12','s1','s2','s3')

    for i in xrange(len(psis)):
        ## rotate uniaxial stress in the lab axes to material axes
        sMaterial                     = rot_6d(sUniaxial,psis[i])
        ysMat, Phi, dPhi, d2Phi       = h48(s=sMaterial)
        ysLab = rot_6d(ysMat, -psis[i])
        dpLab = rot_6d(dPhi,  -psis[i])
        # print '%6f %6f'%(psis[i],Phi)
        phis.append(Phi)

        et = -dpLab[0]-dpLab[1]
        r  = dpLab[1] / et
        rvs.append(r)

        # sMaterial, phi, dphi, d2phi = HillQuad(s=sMaterial,r0=r0,r90=r90)
        # sLab = rot_6d(sMaterial, -psis[i])
        # print 'sLab:',sLab
        # dum1, phi, dum2, dum3 = HillQuad(s=sLab)
        # print 'phi:',phi


    fig = plt.figure(figsize=(12,3))
    ax1  = fig.add_subplot(141,projection='3d')
    ax2  = fig.add_subplot(142)
    ax3  = fig.add_subplot(143)
    ax4  = fig.add_subplot(144)
    import vm_check
    # ax1.plot(x,y,z)
    ax2.plot(psis*180/np.pi,rvs)
    # ax2.set_xlabel(r'$\theta$ [Ang] rotation from $S_1$')
    # ax2.set_xticks(np.arange(0,90.01,45))
    ax3.plot(psis*180/np.pi,phis)
    vm_check.main(ax=ax4)
    ax1.set_xlabel(r'$\sigma_{11}$')
    ax1.set_ylabel(r'$\sigma_{22}$')
    ax1.set_zlabel(r'$\sigma_{12}$')
    ax1.set_aspect('equal')

    ax1.set_ylabel(r'R value')
    ax3.set_ylabel(r'$\Phi$')
    fig.tight_layout()

    for ax in [ax2,ax3]:
        ax.set_xlabel(r'$\theta$ [deg] from RD')
        ax.set_xticks(np.arange(0,90.001,30))

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
    test2() # - uniaxial tension along various directions - not working yet
    #test3()
