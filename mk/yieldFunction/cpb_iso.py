"""
Example of cpb iso Cazacu et al., IJP (2006) 1171-1194
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import mk.library.lib

def eq4(S,a,k):
    ## Find principal values of sigma1, sigma2
    sig,dum = mk.library.lib.convert_6sig_princ(S)
    sig=np.sort(sig)[::-1] ## descending order (large to small)

    ## homogeneous function of degree m
    ##
    """
    Eq (4) in Cazacu et al.
    F=(|S_1|-k*S_1)**a + (|S_2|-k*S_2)**a +(|S_3|-k*S_3)**a
    with S_i, i.e., the principal values of the stress deviator.

    F: the size of the yield locus (should be a constant for any
      stress state that meets the given yield condition)

    F is the homogeneous function of degree a.
    """
    f = 0.
    for i in xrange(len(sig)):
        f = f + (np.abs(sig[i])-k*sig[i])**a
    return f

def eq5a(a,k):
    p1=2./3.
    p2=1./3.

    rst =(   (  (p1*(1+k))**a+2*(p2*(1-k))**a  ) /  ( (p1*(1-k))**a+2*(p2*(1+k))**a   )   ) **(1./a)
    return rst

def eq5b(x,a):
    """
    Convert tension/compression ratio (x) to k
    with the given value of exponent a
    """
    h = eq5c(x,a)
    k = (1-h)/(1+h)
    return k

def eq5c(x,a):
    """
    Convert tension/compression ratio (x) to
    obtain h(value)
    """
    return ((2.**a - 2.*(x**a)) / ( (2*x)**a - 2   )) **(1./a)

def deviator(a):
    if type(a).__name__=='ndarray':
        S=a.copy()
    else:
        S=np.array(a,dtype='float')
    # make sure that the given stress <sigma> is a deviator
    p = S[:3].sum()
    S[:3] = S[:3] - p/3.
    return S

def main(s=[1,0,0,0,0,0]):
    """
    Argument
    --------
    s=[1,0,0,0,0,0]
    """
    if type(s).__name__=='ndarray':
        S=s.copy()
    else:
        S=np.array(s,dtype='float')

    Sdev = deviator(S)

    a = 4.
    k = 2./3.

    f=eq4(S=Sdev.copy(),a=a,k=k)
    f1=(f)**(1./a)
    Sdev2 = deviator(S.copy()/f1)
    f2=eq4(Sdev2.copy(),a=a,k=k)
    print 'f2:',f2

    return S.copy()/f1

def locus(nth=100):
    """
    in-plane biaxial locus

    Argument
    --------
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
    X=[]; Y=[]; fs=[]; fsinv=[]
    sigma=[]
    for i in xrange(len(s)):
        newStress  = main(s[i].copy())
        X.append(newStress[0])
        Y.append(newStress[1])

    fs=np.array(fs)
    fsinv=np.array(fsinv)

    fig=plt.figure()
    ax1=fig.add_subplot(121)
    ax2=fig.add_subplot(122)
    #ax1.plot(X,Y)
    ax1.plot(X,Y)
    ft=dict(fontsize=18)
    ax1.set_xlabel(r'$\bar{\Sigma}_{RD}$',ft)
    ax1.set_ylabel(r'$\bar{\Sigma}_{TD}$',ft)
    ax1.grid('on')
    ax1.set_aspect('equal')

    ax1.set_xlim(-2.0,3.0); ax1.set_ylim(-2.0,3.0)
    ticks=np.linspace(-2,3,6)
    ax1.set_xticks(ticks)
    ax1.set_yticks(ticks)

    fig.savefig('cpb_iso_locus.pdf',bbox_inches='tight')

def inplaneTension():
    psis = np.linspace(0,+np.pi/2.,100)

    ## stress state for inplane tension in the lab axes
    ## is of course uniaxial stress state:
    sLab=np.zeros(6)
    sLab[0]=1.

    phis=[];rvs=[]
    ## rotate this stress state to 'material' axes for
    ## each psi angles and collect results
    for i in xrange(len(psis)):
        sMat = mk.library.lib.rot_6d(sLab, psis[i])
        ysMat = main(s=sMat)
        ysLab = mk.library.lib.rot_6d(ysMat,  -psis[i])
        phis.append(ysLab[0])

    fig=plt.figure()
    ax1=fig.add_subplot(111)
    ax1.plot(psis,phis)
    fig.savefig('cpb_iso_inplanetest.pdf')

if __name__=='__main__':
    locus()
    # inplaneTension()
