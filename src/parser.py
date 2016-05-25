## parsing various files generated in MK
import numpy as np

def plotMat(fn,ax,**kwargs):
    dat=readMatFile(fn)
    eps=dat[0]
    sig=dat[1]
    stress=dat[2:8]
    dphi=dat[8:14]
    s1=sig[::]*stress[0]
    s2=sig[::]*stress[1]
    d1=dphi[0]
    d2=dphi[1]
    rad=20.
    for i in xrange(len(s1)):
        th = np.arctan2(d2[i],d1[i])
        x,y = s1[i],s2[i]
        ax.plot(y,x,'k.',alpha=0.3)
        dx,dy = rad*np.cos(th),rad*np.sin(th)
        ax.arrow(y,x,dy,dx,**kwargs)
        
def readMatFile(fn):
    """
    """
    dat = np.loadtxt(fn).T
    return dat


