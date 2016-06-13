"""
Tune-up parameters of Hill48.

Hill48 parameters can be tuned by r-values or yield stresses
obtained by a series of uniaxial tension tests.
"""
#from for_lib import vm
from yf_for import vm
from yf2 import HillQuad,Hill48, wrapHill48Gen
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from vpscyld.lib_dat import xy2rt
from scipy import interpolate
from mk.tests.mechtests import inplaneTension

pi  = np.pi
sin = np.sin
cos = np.cos
nth = 800

def tuneGenY(y):
    """
    Place-holder to tune H48 yield function parameters
    by fitting with the in-plane variation of uniaxial
    yield stresses <y>

    If only three elements are present,
    h,g,f,n parameters in Hill48 yield function
    are analytically obtained.

    Argument
    --------
    y   - yield stress list
    """
    if len(y)==3:
        pass
    elif len(y)<3:
        raise IOError, 'At least 3 parameters are necessary.'
    else:
        pass


def tuneGenR(r=[2.2,2.0,2.9]):
    """
    Tune based on r-values

    r-values are assumed to be in the form
    of an array with the first element being
    associated with RD and the last TD.
    If only three elements are present,
    h,g,f,n parameters in Hill48 yield function
    are analytically obtained.

    Arguments
    ---------
    r    - r value list
    """
    if len(r)==3:
        r0,r45,r90 = r
        ## Eq 3 in Dasappa et al. IJSS, vol 49, (2012)
        h = r0/(r0+1)
        g = 1 - h
        f = g * r0/r90
        n = (r45+0.5)*(r0/r90+1)*g
    elif len(r)<2:
        raise IOError, 'At least 3 parameters are necessary.'
    elif len(r)==2:
        print 'Warning: only two r-values are given'
        print 'Tune Hill48 prameters by using'
        print "Hill's quadratic plane-stress yield locus"
        f,g,h,n = tuneR2(r0=r[0],r90=r[2])
    else:
        ## case that many more r-values are available
        ## numerically determine h,g,f,n that 'best'
        ## fits the r-value profile.

        ## Initial guess is approximated by
        ## the three r-values.
        r0  = r[0]
        r45 = r[int(len(r)/2.)]
        r90 = r[-1]

        h   = r0/(r0+1)
        g   = 1-h
        f   = g * r0/r90
        n   = (r45+0.5)*(r0/r90+1)*g

        x0 =[h,g,f,n]

        objf = returnObjRV(rv=r,fYLD=Hill48)
        res = minimize(fun=objf, x0=x0,method='BFGS',
                       jac=False,tol=1e-10,options=dict(maxiter=20))

        popt = res.x
        n_it = res.nit
        fopt = res.fun

        f,g,h,n=popt
        y=(g+h)
        params = np.array([f,g,h,n])
        params = params / y
        f,g,h,n = params

    print 'Hill48 parameter tuning in tuneH48.tuneGenR'
    print ('%7s'*4)%('f','g','h','n')
    print ('%7.3f'*4)%(f,g,h,n)

    return f,g,h,n

def tuneR2(r0=1.,r90=1.):
    """
    Tune Hill48 using Hill Quadratic plane-strss yield locus
    that is characterized by two r-values (r0 and r90)

    Arguments
    ---------
    r0
    r90
    """
    th=np.linspace(-pi,+pi,nth)
    x=cos(th); y=sin(th)
    z=np.zeros(len(th))
    s=np.array([x,y,z,z,z,z]).T
    X=[]; Y=[]

    for i in xrange(len(s)):
        rst = HillQuad(s[i],r0=r0,r90=r90)
        ys = rst[0]
        X.append(ys[0])
        Y.append(ys[1])

    ysHQ = np.array([X,Y])
    objf = returnObjYS(ref=ysHQ,fYLD=Hill48)

    res = minimize(
        fun=objf,
        x0=[0.277,0.357,0.749,1.464],
        bounds=((0.2,0.8),(0.2,0.8),(0.2,2),(0.5,1.5)),

        # method='BFGS',
        method='Nelder-Mead',

        jac=False,tol=1e-20,
        options=dict(maxiter=400))

    popt = res.x
    n_it = res.nit
    fopt = res.fun

    # print 'popt:',popt
    f,g,h,n=popt

    y=(g+h)
    params = np.array([f,g,h,n])
    params = params / y
    f,g,h,n = params

    print 'Hill48 parameter tuning in tuneH48.tuneR2'
    print ('%6s'*4)%('f','g','h','n')
    print ('%6.3f'*4)%(f,g,h,n)
    return f,g,h,n

## Generate objective function that compares yield locus
## in the plane-stress space.
def returnObjYS(ref=None,fYLD=Hill48):
    """
    Given the reference data <ref>,
    provide an objective function that is
    subject to optimized when a proper
    set of parameter is obtained.

    Arguments
    ---------
    ref
    fYLD
    """
    def objf(xs):
        """
        The objective function
        generated in mk.yieldFunction.tuneH48.returnObjYS

        Argument
        --------
        xs -- <f,g,h,n> the four Hill48 parameters
        in the plane-stress space.
        """
        th=np.linspace(-pi,+pi,nth)
        x=cos(th); y=sin(th)
        z=np.zeros(len(th))
        s=np.array([x,y,z,z,z,z]).T
        f,g,h,n = xs

        X=[]; Y=[]
        for i in xrange(len(s)):
            rst=fYLD(s[i],f,g,h,n)
            ys = rst[0]
            X.append(ys[0])
            Y.append(ys[1])

        X=np.array(X)
        Y=np.array(Y)

        ## rst to th-r coordinate
        r,th = xy2rt(X,Y)
        R,TH = xy2rt(ref[0],ref[1])
        diff = ((r-R)**2).sum()/(len(th)-1)
        return diff
    return objf

def returnObjRV(rv,fYLD=Hill48):
    """
    Arguments
    ---------
    rv
    fYLD
    """
    def objf(xs):
        nth  = len(rv)
        psis_ref = np.linspace(0,np.pi/2.,nth)
        f, g, h, n = xs
        psis, rvs, phis = inplaneTension(fYLD=fYLD,f=f,g=g,h=h,n=n)
        funcINT = interpolate.interp1d(x=psis,y=rvs)
        rvs_ref = funcINT(psis_ref)
        diff = np.sqrt(((rvs_ref - rv)**2).sum())/(nth-1.)
        return diff
    return objf
