## parsing various files generated in MK
import numpy as np

def plotMat(fn,ax,**kwargs):
    """
    Plot loading history of regions A and B

    Arguments
    ---------
    fn
    ax
    **kwargs  will be arguments for matplotlib plot and arrow
    """
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
        ax.plot(y,x,'.',alpha=0.3,**kwargs)
        dx,dy = rad*np.cos(th),rad*np.sin(th)
        
        ## may be just ticks instead of arrows
        # ax.arrow(y,x,dy,dx,**kwargs)
        ax.plot([x,x+dx],[y,y+dy],**kwargs)

def readMatFile(fn):
    """
    material log file <fn> is parsed to data.

    Argument
    --------
    fn  <file name> of the material history log file
    """
    ## No head is appended thus np.loadtxt is good enough
    dat = np.loadtxt(fn).T
    return dat
