import matplotlib as mpl
mpl.use('Agg') ## In case X-window is not available.
from for_lib import vm
from yf2 import wrapHill48, wrapHillQuad
import numpy as np
import matplotlib.pyplot as plt
pi=np.pi
sin=np.sin
cos=np.cos

def locus(func):
    """
    in-plane biaxial locus
    """
    th=np.linspace(-pi,pi,100)
    x=cos(th); y=sin(th)
    z=np.zeros(len(th))
    s=np.array([x,y,z,z,z,z]).T
    X=[]; Y=[]
    for i in xrange(len(s)):
        rst = func(s[i])
        ys = rst[0]
        X.append(ys[0])
        Y.append(ys[1])
    return X,Y

def main(ax=None):
    if type(ax).__name__=='NoneType':
        fig = plt.figure()
        ax  = fig.add_subplot(111)

    ## two parameters from Jeong et al. Acta Mat 112, 2016
    ## for the interstitial-free steel
    r0  = 2.1; r90 = 2.7
    hq = wrapHillQuad(r0=r0,r90=r90)
    h48= wrapHill48(r0=r0,r90=r90) ## tune f,g,h by HillQuad

    funcs = vm, hq, h48
    labs  = ['von Mises','Hill Quad','Hill48']
    ls    = ['-','--',':']
    for i in xrange(len(funcs)):
        x,y = locus(funcs[i])
        ax.plot(x,y,label=labs[i],ls=ls[i])

    ax.legend(loc='best')
    fn='vm_check.pdf'
    print '%s has been saved'%fn
    ax.set_xlim(-0.25,1.75)
    ax.set_ylim(-0.25,1.75)
    ax.set_aspect('equal')

    try:
        fig.savefig(fn)
    except:
        pass

if __name__=='__main__':
    main()
