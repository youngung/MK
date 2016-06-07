import matplotlib as mpl
mpl.use('Agg') ## In case X-window is not available.
from for_lib import vm
from yf2 import HillQuad, Hill48
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

    funcs=(vm,HillQuad, Hill48)
    labs =['von Mises','HillQuadratic','Hill48']
    for i in xrange(len(funcs)):
        x,y = locus(funcs[i])
        ax.plot(x,y,label=labs[i])

    ax.legend(loc='best')
    fn='vm_check.pdf'
    print '%s has been saved'%fn
    ax.set_xlim(-1.5,1.5)
    ax.set_ylim(-1.5,1.5)
    ax.set_aspect('equal')

    try:
        fig.savefig(fn)
    except:
        pass

if __name__=='__main__':
    main()
