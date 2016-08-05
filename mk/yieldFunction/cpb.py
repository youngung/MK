"""
Testing Plunkett, Cazacu, Barlat yield function
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import mk.library.lib

def data(iopt=0):
    """
    Data library
    """
    if iopt==0:
        pass
    elif iopt==1:
        ## 2090-T3 Al alloy (Plunkett et al., 2008, IJP vol 24 pp847-866)
        c11  =  0.453
        c12  = -0.841
        c13  = -1.248
        c22  = -1.058
        c23  = -2.284
        c33  = -3.201
        c44  =  0.000
        c55  =  0.000
        c66  =  1.026

        c11_ =  0.453
        c12_ = -0.705
        c13_ =  1.148
        c22_ =  0.139
        c23_ = -0.519
        c33_ =  0.878
        c44_ =  0.000
        c55_ =  0.000
        c66_ =  1.978
        k =0.054
        k_=0.027
        a =2.
    elif iopt==2:
        ## AZ31B Mg
        c11  =  0.771
        c12  = -0.914
        c13  = -1.312
        c22  =  0.654
        c23  = -1.209
        c33  =  0.293
        c44  =  0.000
        c55  =  0.000
        c66  =  1.617

        c11_ =  0.771
        c12_ =  0.337
        c13_ = -0.764
        c22_ =  0.646
        c23_ = -0.401
        c33_ = -0.022
        c44_ =  0.000
        c55_ =  0.000
        c66_ =  0.499
        k =0.295
        k_=0.900
        a =2.
    else:
        raise IOError, 'unexpected option given'

    c  = [[ c11,  c12,  c13,   0.,   0.,   0.],
          [ c12,  c22,  c23,   0.,   0.,   0.],
          [ c13,  c23,  c33,   0.,   0.,   0.],
          [  0.,   0.,   0.,  c44,   0.,   0.],
          [  0.,   0.,   0.,   0.,  c55,   0.],
          [  0.,   0.,   0.,   0.,   0.,  c66]]
    c_ = [[c11_, c12_, c13_,   0.,   0.,   0.],
          [c12_, c22_, c23_,   0.,   0.,   0.],
          [c13_, c23_, c33_,   0.,   0.,   0.],
          [  0.,   0.,   0., c44_,   0.,   0.],
          [  0.,   0.,   0.,   0., c55_,   0.],
          [  0.,   0.,   0.,   0.,   0., c66_]]
    c  = np.array(c)
    c_ = np.array(c_)
    return c, c_, k, k_, a


def main(s=[1,0,0,0,0,0]):
    """
    Arguments
    ---------
    s = [1,0,0,0,0,0]
    """
    c1,c2,k1,k2,a=data(1)

    ## transform the given stress s to sigmas
    sigma1=np.tensordot(c1,s,axes=([1,0]))
    sigma2=np.tensordot(c2,s,axes=([1,0]))

    ## Find principal values of sigma1, sigma2
    sig1,rot1 = mk.library.lib.convert_6sig_princ(sigma1)
    sig2,rot2 = mk.library.lib.convert_6sig_princ(sigma2)

    # phi
    f = 0.
    for i in xrange(3):
        f = f + (np.abs(sig1[i])-k1*sig1[i])**a
        f = f + (np.abs(sig2[i])-k2*sig2[i])**a

    print 'sigma1:',sigma1
    print 'sigma2:',sigma2
    print 'sig1:',sig1
    print 'sig2:',sig2
    print 'f:', f

    ct = (f/2.) ** (1./a)
    s=s*ct
    return s

def locus(nth=100,iplot=False):
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
    X=[]; Y=[]
    for i in xrange(len(s)):
        ys = main(s[i])
        X.append(ys[0])
        Y.append(ys[1])

    if iplot:
        fig=plt.figure()
        ax1=fig.add_subplot(111)
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
        fn='cpb_locus.pdf'
        fig.savefig(fn,bbox_inches='tight')
        print '%s has been saved'%fn
        pass
    return X,Y


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
    fig.savefig('cpb_inplanetest.pdf')


if __name__=='__main__':
    locus()
    # inplaneTension()
