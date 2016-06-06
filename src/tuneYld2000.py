"""
Cases for tuning-up YLD2000-2D parameters
"""
from yf2 import wrapYLD
from vm_check import locus
from vpscyld.lib_dat import xy2rt, rad_xy2
import numpy as np
import matplotlib.pyplot as plt

def returnObjf_FullyYS(YS):
    """
    Tune purely based on YS
    YS; the given ys locus in (th,r) coordinates
    """
    def returnDiff(params):
        """
        Tune using the given in-plane locus

        Arguments
        =========
        params [:4]  : r0,r45,r90,rb
        params [4:7] :    y45,y90,yb
        """
        ## Find 8 parameters...
        paramR = np.ones(4)
        paramY = np.ones(4)

        paramR[:4]  = params[:4]
        paramY[1:4] = params[4:7]
        k=2; m=8 ## the exponent is fixed?

        yldFunc = wrapYLD(r=paramR,y=paramY,m=m,k=k)
        if type(yldFunc).__name__=='int':
            if yldFunc==-1:
                return 0.5

        yldLocus = locus(yldFunc,1000) ##
        ## convert yldLocus in to (th,r) coordinates
        r,th=xy2rt(np.array(yldLocus[0]),np.array(yldLocus[1]))

        th0,r0 = YS## reference
        radIntp = rad_xy2(r,th,th0)
        radIntp = np.array(radIntp)

        return np.sqrt(((radIntp-r0)**2).sum())/ (len(r)-1)
    return returnDiff


def returnObjf_YSRV(YS,rv):
    """
    Tune based on YS and R-values
    YS; the given ys locus in (th,r) coordinates
    rv: r0, r45, r90
    """
    def returnDiff(params):
        """
        Tune using the given in-plane locus

        Arguments
        =========
        params [0]  :          rb
        params [1:4]:  y45,y90,yb
        """
        ## Find 4 parameters...
        paramR = np.ones(4)
        paramY = np.ones(4)

        paramR[0:4]  = rv[0],rv[1],rv[2],params[0]
        paramY[1:4] = params[1:4]
        k=2; m=8 ## the exponent is fixed?

        yldFunc = wrapYLD(r=paramR,y=paramY,m=m,k=k)
        if type(yldFunc).__name__=='int':
            if yldFunc==-1:
                return 1

        yldLocus = locus(yldFunc,1000) ##
        ## convert yldLocus in to (th,r) coordinates
        r,th=xy2rt(np.array(yldLocus[0]),np.array(yldLocus[1]))

        th0,r0 = YS## reference
        radIntp = rad_xy2(r,th,th0)
        radIntp = np.array(radIntp)

        return np.sqrt(((radIntp-r0)**2).sum())/ (len(r)-1)
    return returnDiff

def extractParamsFromYS(locus_planestress,psi_uni,rvs_uni,ys_uni):
    """
    Given plane-stress yield locus, extract parameters.
    """
    from vpscyld import lib_dat, calc_normal
    ## normalize the locus by stress along axis 1
    ysx, ysy =lib_dat.nys_th(np.array(locus_planestress[0]),np.array(locus_planestress[1]),t=0)
    ysr, yst = lib_dat.xy2rt(ysx,ysy)

    yb = lib_dat.rad_xy(ysr,yst,45*np.pi/180.)  ## balanced biaxial yield stress...
    yb = yb / np.sqrt(2)
    y0 = lib_dat.rad_xy(ysr,yst,0)  ## balanced biaxial yield stress...
    y90 = lib_dat.rad_xy(ysr,yst,90*np.pi/180.)  ## balanced biaxial yield stress...

    th = calc_normal.main_var(locus_xy=[ysr,yst],psi_ref=[45.])
    rb = 1./np.tan(th[0])

    y0_uni  = ys_uni[0]
    y90_uni = ys_uni[-1]
    r0_uni  = rvs_uni[0]
    r90_uni = rvs_uni[-1]
    r45_uni = np.interp(45*np.pi/180., psi_uni, rvs_uni)
    y45_uni = np.interp(45*np.pi/180., psi_uni, ys_uni)


    if abs(y0_uni-y0)>1e-3:
        print 'y0 diff'
        print y0_uni,y0
        # raise IOError, 'Consider using finer step for yield locus'
    if abs(y90_uni-y90)>1e-3:
        print 'y90 diff'
        print y90_uni, y90
        # raise IOError, 'Consider using finer step for yield locus'

    return y0,y90,yb,rb, y45_uni, r45_uni

def case1(YS):
    """
    Given, normalized in-plane yield locus <YS>
    Find yld2000-2D parameters that fits <YS>
    """
    method = 'Nelder-Mead'
    # method = 'BFGS'
    from scipy.optimize import minimize
    objf = returnObjf_FullyYS(YS)
    res = minimize(objf,x0=[1,1,1,1,1,1,1,1],method=method,jac=False,
                   tol=1e-10,options=dict(maxiter=4000))
    popt = res.x
    n_it = res.nit
    fopt = res.fun

    print 'popt:',popt
    print 'fopt:',fopt

    return popt

def case2(YS,rv):
    """
    Given, normalized in-plane yield locus <YS> and three r-values,
    i.e., r0, r45, r90,
    Find yld2000-2D parameters (rb, y45, y90, yb) that fits <YS>
    """
    method = 'Nelder-Mead'
    # method = 'BFGS'
    from scipy.optimize import minimize
    objf = returnObjf_YSRV(YS,rv)
    res = minimize(objf,x0=[1,1,1,1.35],method=method,jac=False,
                   tol=1e-10,options=dict(maxiter=4000))
    popt = res.x
    n_it = res.nit
    fopt = res.fun

    print 'n_it:',n_it
    print 'popt:',popt
    print 'fopt:',fopt

    return popt

def ex2():
    """
    Characterize through in-plane tension tests, and rb, yb
    """
    import tuneH48, yf2, mechtests
    rv=[2.2,2.0,2.9]
    f,g,h,n = tuneH48.tuneGenR(r=rv)
    # f,g,h,n = 0.237,0.313,0.687,1.374
    yfunc_H48 = yf2.wrapHill48Gen(f,g,h,n)
    locus_H48 = locus(yfunc_H48,30000)

    psi_uni_H48, rvs_uni_H48, ys_uni_H48 \
        = mechtests.inplaneTension(yfunc_H48,1)
    y0, y90, yb, rb, y45, r45 = extractParamsFromYS(
        locus_H48,psi_uni_H48,rvs_uni_H48,ys_uni_H48)
    print y0, y90, yb, rb, y45, r45

    r_yld     = [rv[0],     rv[1],  rv[2], rb]
    y_yld     = [ys_uni_H48[0], y45,   ys_uni_H48[-1], yb]
    yfunc_yld = yf2.wrapYLD(r=r_yld, y=y_yld, m=6.)

    if type(yfunc_yld).__name__=='int':
        raise IOError

    locus_yld = locus(yfunc_yld,30000)
    psi_uni_yld, rvs_uni_yld, ys_uni_yld \
        = mechtests.inplaneTension(yfunc_yld,1)

    fig =plt.figure(figsize=(11,3.));
    ax1=fig.add_subplot(131);ax2=fig.add_subplot(132)
    ax3=fig.add_subplot(133)

    ax1.plot(psi_uni_H48*180/np.pi,ys_uni_H48,label='Hill48')
    ax2.plot(psi_uni_H48*180/np.pi,rvs_uni_H48,label='Hill48')
    ax1.plot(psi_uni_yld*180/np.pi,ys_uni_yld,'--',label='YLD2000')
    ax2.plot(psi_uni_yld*180/np.pi,rvs_uni_yld,'--',label='YLD2000')
    ax3.plot(locus_H48[0],locus_H48[1],label='Hill48')
    ax3.plot(locus_yld[0],locus_yld[1],'--',label='YLD2000')
    ax1.legend()
    ax2.legend()
    ax3.legend()
    ax3.set_xlim(0.,)
    ax3.set_ylim(0.,)
    ax3.set_aspect('equal')
    fig.tight_layout()
    fig.savefig('in-plane-tension.pdf',bbox_to_inches='tight')

def ex1():
    """
    Excercise - Tune H48 using three r-values and calculate YS.
    Using the three given r-values used in H48
    find rb, y0, y45, y90, and yb by fitting YS of H48 with
    that of yld2000-2d
    """
    import tuneH48
    import yf2
    rv=[2.2,2.0,2.9]
    f,g,h,n = tuneH48.tuneGenR(r=rv)
    yfunc_H48 = yf2.wrapHill48Gen(f,g,h,n)
    YS_H48X, YS_H48Y = locus(yfunc_H48,1000)
    r, th = xy2rt(np.array(YS_H48X), np.array(YS_H48Y))

    ## popt = [rv, y45, y90, yb]

    popt = case2([r,th], rv=rv)
    rv.append(popt[0])
    ys = np.ones(4)
    ys[1:] = popt[1:]

    ##
    print 'rv:', rv
    print 'ys:', ys
    yfunc_yld2000 = wrapYLD(r=rv,y=ys,m=8)

    import matplotlib.pyplot as plt
    fig = plt.figure(); ax=fig.add_subplot(111)

    YS_YLDX, YS_YLDY = locus(yfunc_yld2000,1000)
    ax.plot(YS_H48X,YS_H48Y,'-',label='Hill48')
    ax.plot(YS_YLDX,YS_YLDY,'-',label='YLD2000')
    ax.set_xlim(0.,)
    ax.set_ylim(0.,)
    ax.set_aspect('equal')
    ax.legend(loc='best')
    # print YS_YLDX
    # print YS_YLDY
    fig.savefig('Yld2000-Hill48.pdf')


def test():
    # ex1()
    ex2()

if __name__=='__main__':
    test()
