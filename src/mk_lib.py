from for_lib import vm
import numpy as np
import time
from MP import progress_bar
uet = progress_bar.update_elapsed_time
cos=np.cos
sin=np.sin
tan=np.tan
log=np.log
atan2=np.arctan2
sqrt=np.sqrt

def findStressOnYS(f_yld,stressLeft,stressRight,pth,verbose):
    """
    Given the yield function <f_yld> find the stress
    that corresponding to a certain strain rate <pth>
    on the basis of the associatd flow rule.

    Arguments
    ---------
    f_yld:
        yield function that is supposed to provide derivatives
        to obtain strain rate vectors on the basis of the associated
        flow rule

    (stressLeft, stressRight): 
        bounds that approximately embrace the 
        presumable stress location to improve
        numerical stability

    pth      : strain rate vector

    verbose 
    """
    tol = 1e-10  # tolerance

    ## normalize pth just in case.
    pth  = np.array(pth)/sqrt((np.array(pth)**2).sum())

    
    if pth[0]>pth[1]:
        majorDir = 'RD'
        _rho_ = pth[1]/pth[0]
    else:
        majorDir = 'TD'
        _rho_ = pth[0]/pth[1]

    print 'strain rate vector:',\
        ('%6.3f '*2)%(pth[0],pth[1]), 'rho: %6.3f  %s'%(
            _rho_,majorDir)

    diff = tol*100 ## sufficiently large initial tolerance

    ##
    if verbose: print (9*'%6s ')%(
        'it','phi','sigL','sigR','SL0','SL1',
        'SR0','SR1','diff')
    it = 0
    while diff>tol:
        it = it + 1
        s  = (stressLeft[::]+stressRight[::])/2.
        s,phi,dphi,d2phi = f_yld(s) # s on the locus
        rac     = sqrt(dphi[0]**2 + dphi[1]**2)
        dphi[0] = dphi[0]/rac
        dphi[1] = dphi[1]/rac

        ## narrow down stressLeft-stressRight bounds
        ## project the given edot and the edot from yield function
        ## see if the angle is acute or obtuse

        if dphi[0]*pth[1]-dphi[1]*pth[0]>=0:
            stressLeft[:] = s[:] ## bring up stressLeft
        else:
            stressRight[:] = s[:] ## bring down stressRight

        if verbose:
            print ('%6i '+7*'%6.3f '+'%9.3e')%(
                it,phi,s[0],s[1],stressLeft[0],stressLeft[1],
                stressRight[0],stressRight[1],diff)

        if (it>100): raise IOError, 'Could not find the proper s'
        diff = sqrt(((stressLeft-stressRight)**2).sum())

    rho   = dphi[1]/dphi[0]
    alpha = s[1]   /s[0]

    if verbose: 
        print '# of iterations:',it
        print ('%8s %5.2f')%('rho  :',rho)
        print ('%8s %5.2f')%('alpha:',alpha)

    return s,dphi
