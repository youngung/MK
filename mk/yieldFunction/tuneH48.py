"""
Tune-up parameters of Hill48 using Hill Quad
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
    Place-holder to tune based on yield functions
    """
    if len(y)==3:
        pass
    pass

def tuneGenR(r=[2.2,2.0,2.9]):
    """
    Tune based on r-values

    Arguments
    ---------
    r    - r value list
    fYLD - yield function
    """
    if len(r)==3:
        r0,r45,r90 = r
        # ## original formula in R. Hill, JMPS, V38, 1990
        # h = 2*r0 / (2*r0+1)
        # g = 1 - h
        # f = h/(2.*r90)
        # n = (2*r45+1) * (g+f) / 2.

        ## Eq 3 in Dasappa et al. IJSS, vol 49, (2012)
        h = r0/(r0+1)
        g = 1 - h
        f = g * r0/r90
        n = (r45+0.5)*(r0/r90+1)*g
    else:
        ## approximate...
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
        #               jac=False,tol=1e-20,options=dict(maxiter=400))

        popt = res.x
        n_it = res.nit
        fopt = res.fun

        # print 'popt:',popt
        # print 'fopt:',fopt

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
    Arguments
    ---------
    ref
    fYLD
    """
    def objf(xs):
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
