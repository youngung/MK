import matplotlib as mpl
mpl.use('Agg') ## In case X-window is not available.
from for_lib import vm
from yf2 import wrapHill48, wrapHill48R, wrapHillQuad, wrapYLD
import numpy as np
import matplotlib.pyplot as plt
pi=np.pi
sin=np.sin
cos=np.cos

def locus(func,nth=100):
    """
    in-plane biaxial locus
    """
    th=np.linspace(-pi,pi,nth)
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

    import time
    t0=time.time()

    ## two parameters from Jeong et al. Acta Mat 112, 2016
    ## for the interstitial-free steel
    r0   = 2.20
    r45  = 2.0
    r90  = 2.9
    rb   = 1.

    y0   =1
    y45  =1
    y90  =1
    yb   =1.

    # hq   = wrapHillQuad(r0=r0,r90=r90)
    # h48  = wrapHill48(r0=r0,r90=r90) ## tune f,g,h by HillQuad
    h48g = wrapHill48R([r0,r45,r90])
    yld2000 = wrapYLD(r=[r0,r45,r90,rb],y=[y0,y45,y90,yb],m=6,k=2)

    funcs =    vm,    h48g,        yld2000
    labs  = ['von Mises', 'Hill48R','yld2000']
    ls    = ['-','--','-.',':']
    xs=[];ys=[]
    t_indv=[]
    for i in xrange(len(funcs)):
        t_indv0= time.time()
        x,y = locus(funcs[i])
        t_indv.append(time.time()-t_indv0)
        xs.append(x)
        ys.append(y)

    for i in xrange(len(funcs)):
        print 'time %s :'%labs[i], t_indv[i]
    print 'total time:',time.time()-t0
    for i in xrange(len(funcs)):
        ax.plot(xs[i],ys[i],label=labs[i],ls=ls[i])

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
