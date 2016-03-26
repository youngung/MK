"""
Common libraries
"""
import os
import numpy as np
pi = np.pi
from MP.mat import voigt
ijv=voigt.ijv
vij=voigt.vij
from MP.lib import mpl_lib
ticks_bins = mpl_lib.ticks_bins

def draw_guide(ax,r_line = [-0.5,0. ,1],max_r=2,
               ls='--',color='k',alpha=0.5):
    """
    Maximum should be a radius...
    """
    # guide lines for probed paths
    xlim=ax.get_xlim(); ylim=ax.get_ylim()
    for i in xrange(len(r_line)):
        r = r_line[i]
        mx=max_r
        mx = mx/np.sqrt(1.+r**2)
        ys = np.linspace(0.,mx)
        xs = r * ys

        ax.plot(xs,ys,ls=ls,color=color,alpha=alpha)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

def rot(psi):
    psi = psi * np.pi/180.
    r = np.zeros((3,3)); c = np.cos(psi); s = np.sin(psi)
    r[0,0]= c;  r[0,1]=-s
    r[1,0]= s;  r[1,1]= c
    r[2,2]= 1.
    return r

def rot_vec(r,vect):
    """
    v[i] = r[i,j]*vect[j]
    """
    v=np.zeros(3)
    for i in xrange(3):
        for j in xrange(3):
            v[i]=v[i]+r[i,j]*vect[j]
    return v

def rot_tensor(a,psi):
    """
    Arguments
    ---------
    a
    psi (in degree)
    """
    a=np.array(a)
    b=np.zeros((3,3))
    r=rot(psi)
    for i in xrange(3):
        for j in xrange(3):
            b[i,j] = 0.
            for k in xrange(3):
                for l in xrange(3):
                    b[i,j] = b[i,j] + r[i,k] * a[k,l] * r[j,l]
    return b

def th_2_planestress(th):
    """
    Given theta, return the stress.

    Argument
    --------
    th
    """
    sigma=np.zeros(6)
    sigma[0]=np.cos(th)
    sigma[1]=np.sin(th)
    return sigma

def th_planestress(th,yfunc,**kwargs):
    """
    Return stress tensors that gives the same
    size (value) of phi

    Argument
    --------
    th
    yfunc
    **kwargs for the given yfunc

    - yfunc is assumed to take stress and **kwargs
    """
    Sigma = th_2_planestress(th)
    y     = yfunc(Sigma,**kwargs)
    return Sigma / y

def th_planestress_c(th,yfunc):
    """
    Return stress tensors that gives the same
    size (value) of phi

    Argument
    --------
    th
    yfunc
    **kwargs for the given yfunc

    - yfunc is assumed to take stress and **kwargs
    """
    Sigma = th_2_planestress(th)
    y     = yfunc(Sigma)
    return Sigma / y


def convert_6sig_princ(s6):
    """
    Convert 6D stress vectors to
    Eigen vetors and values

    Argument
    --------
    s6

    Returns
    -------
    w : Eigen values
    v : Eigen vectors
    """
    sig33=np.zeros((3,3))
    for k in xrange(6):
        i,j = ijv[:,k]
        sig33[i,j] = s6[k]
        if i!=j: sig33[j,i]=s6[k]

    w,v = np.linalg.eig(sig33)
    return w,v

def convert_sig33_sig6(sig33):
    s6=np.zeros(6)
    for k in xrange(6):
        i,j = ijv[:,k]
        s6[k] = sig33[i,j]
    return s6

def convert_sig6_sig33(sig6):
    s33=np.zeros((3,3))
    for k in xrange(6):
        i,j = ijv[:,k]
        s33[i,j] = sig6[k]
        if i!=j: s33[j,i] = sig6[k]
    return s33

def proj_sig33_nt(sig33,n,t):
    """
    Arguments
    ---------
    sig33 [3,3]
    n [3]
    t [3]
    """
    snt = np.dot(t,np.dot(sig33,n))
    snn = np.dot(n,np.dot(sig33,n))
    return snn, snt


## alias
c6p  = convert_6sig_princ
c2s6 = convert_sig33_sig6
s62c = convert_sig6_sig33

def ys_temp(ax):
    """
    Plane stress space (11,22)
    """
    ax.grid()
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\Sigma_\mathrm{11}$',fontsize=17)
    ax.set_ylabel(r'$\Sigma_\mathrm{22}$',fontsize=17)
    ticks_bins(ax,axis='x',n=4)
    ticks_bins(ax,axis='y',n=4)

def ys_tempr(ax):
    """
    Plane stress space (11,22)
    """
    ax.grid()
    ax.set_aspect('equal')
    ax.set_ylabel(r'$\Sigma_\mathrm{11}$',fontsize=17)
    ax.set_xlabel(r'$\Sigma_\mathrm{22}$',fontsize=17)
    ticks_bins(ax,axis='x',n=4)
    ticks_bins(ax,axis='y',n=4)
    draw_guide(ax,r_line=[0,1,2],max_r=500.)

def es_temp(ax):
    """
    2D strain space (e11, e22)
    """
    ax.grid()
    ax.set_aspect('equal')
    ax.set_xlabel(r'$\mathrm{E_{11}}$',fontsize=17)
    ax.set_ylabel(r'$\mathrm{E_{22}}$',fontsize=17)
    ticks_bins(ax,axis='x',n=4)
    ticks_bins(ax,axis='y',n=4)


def es_tempr(ax):
    """
    2D strain space (e11, e22)
    """
    ax.grid()
    ax.set_aspect('equal')
    ax.set_ylabel(r'$\mathrm{E_{11}}$',fontsize=17)
    ax.set_xlabel(r'$\mathrm{E_{22}}$',fontsize=17)
    draw_guide(ax,r_line=[-0.5,0,1,2,2.5])
    ticks_bins(ax,axis='x',n=4)
    ticks_bins(ax,axis='y',n=4)

def ss_temp(ax):
    """
    Hardening curve
    """
    ax.grid()
    ax.set_xlabel(r'$\bar{\epsilon}$',fontsize=17)
    ax.set_ylabel(r'$\bar{\sigma}$',fontsize=17)
    ticks_bins(ax,axis='x',n=4)
    ticks_bins(ax,axis='y',n=4)

def pi_proj(sd):
    """
    Deviatoric stress to pi-plane projection
    """
    sq6 = np.sqrt(6.)
    sq2 = np.sqrt(2.)
    x = 2.*sd[0]/sq6 - sd[1] / sq6 - sd[2] / sq6
    y =                sd[1] / sq2 - sd[2] / sq2
    return x,y

def devit(x,p=0.):
    """
    Convert it to a deviatoric space with
    hydrostatic stress of <p>

    Argument
    --------
    x  : Conver stress <x> to a deviatric
    p  : pressure (optional, if non-zero,
                   translate the stress along hydrostatic axis)
    """
    x=np.array(x,dtype='float')
    m=x[:3].sum()/3.
    x[:3]=x[:3]-m
    if p!=0: x[:3]=x[:3]+p/3.
    return x

def y_locus(nths,yfunc,**kwargs):
    ths = np.linspace(-pi,pi,nths)
    locus_ps = np.zeros((2,nths))
    xy=[]
    for i in xrange(len(ths)):
        ys = th_planestress(ths[i],yfunc,**kwargs)
        locus_ps[:,i]=ys[0],ys[1]
        sd = np.zeros(6)
        sd[0],sd[1] = ys[0],ys[1]
        sd = devit(sd)
        x, y = pi_proj(sd)
        xy.append([x,y])
        locus_pi=np.array(xy).T
    return locus_ps, locus_pi

def y_locus_c(nths,yfunc):
    ths = np.linspace(-pi,pi,nths)
    locus_ps = np.zeros((2,nths))
    xy=[]
    for i in xrange(len(ths)):
        ys = th_planestress_c(ths[i],yfunc)
        locus_ps[:,i]=ys[0],ys[1]
        sd = np.zeros(6)
        sd[0],sd[1] = ys[0],ys[1]
        sd = devit(sd)
        x, y = pi_proj(sd)
        xy.append([x,y])
        locus_pi=np.array(xy).T
    return locus_ps, locus_pi

def norm_vec(a):
    return np.sqrt((a**2).sum())

def assoc_flow(s6,lamb,yfunc,**kwargs):
    """
    Argument
    --------
    s6  (6D cauchy stress)
    lambda (proportional factor in
            the associated flow rule equation)
    yfunc  (yield function)
    **kwargs (key-worded arguments for yfunc

    Returns
    -------
    edot in 6D (strain rate vector)
    """
    dlt = 1e-10
    s1  = np.zeros(6)
    s2  = np.zeros(6)
    dki = np.identity(6)
    e_k = np.zeros(6)

    ## S should be 'deviatoric'
    for k in xrange(6):
        dum=0.
        s1=np.zeros(6);
        s2=np.zeros(6);
        for i in xrange(6):
            s1[i] = s6[i] + dki[k,i] * dlt
            s2[i] = s6[i] - dki[k,i] * dlt
        e_k[k] = lamb*(yfunc(s1,**kwargs)
                       - yfunc(s2,**kwargs))/(2*dlt)

    ## to by-pass the issue with Quad-Hill
    ## only the direction matters here.
    e_k[2] = - e_k[:2].sum()
    return e_k

def assoc_flow_c(s6,lamb,yfunc):
    """
    Argument
    --------
    s6  (6D cauchy stress)
    lambda (proportional factor in
            the associated flow rule equation)
    yfunc  (yield function)

    Returns
    -------
    edot in 6D (strain rate vector)
    """
    dlt = 1e-10
    s1  = np.zeros(6)
    s2  = np.zeros(6)
    dki = np.identity(6)
    e_k = np.zeros(6)

    ## S should be 'deviatoric'
    for k in xrange(6):
        dum=0.
        s1=np.zeros(6);
        s2=np.zeros(6);
        for i in xrange(6):
            s1[i] = s6[i] + dki[k,i] * dlt
            s2[i] = s6[i] - dki[k,i] * dlt
        e_k[k] = lamb*(yfunc(s1)
                       - yfunc(s2))/(2*dlt)

    ## to by-pass the issue with Quad-Hill
    ## only the direction matters here.
    e_k[2] = - e_k[:2].sum()
    return e_k


def alph2sig(alpha,beta):
    if alpha<=1.:
        return alph2sig1(alpha,beta)
    if alpha>1.:
        return alph2sig2(2-alpha,beta)

def alph2sig1(alpha,beta):
    """
            |   1   beta   0 |
    sigma = | beta  alpha  0 |
            |   0     0    0 |
    """
    sigma = np.zeros((3,3))
    sigma[0,0] = 1.
    sigma[1,1] = alpha
    sigma[0,1] = beta
    sigma[1,0] = beta
    return sigma

def alph2sig2(alpha,beta):
    """
            |alpha  beta   0 |
    sigma = | beta    1    0 |
            |   0     0    0 |
    """
    sigma = np.zeros((3,3))
    sigma[1,1] = 1.
    sigma[0,0] = alpha
    sigma[0,1] = beta
    sigma[1,0] = beta
    return sigma


def alph2sig6(alpha,beta):
    """
    (alpha,beta) to sigma6
    """
    s33 = alph2sig(alpha,beta)
    s6  = c2s6(s33)
    return s6

def alph2eps(alpha,beta,potential,**kwargs):
    """
    """
    ## 6D stress
    cs6 = alph2sig6(alpha,beta)
    ## 6D strain rate
    de6 = assoc_flow(cs6,1.,potential,**kwargs)
    vol = de6[:3].sum()
    de6[:3] = de6[:3] - vol/3.
    return de6

def alph2eps_c(alpha,beta,potential):
    """
    based on ready characterizied potential
    """
    ## 6D stress
    cs6 = alph2sig6(alpha,beta)
    ## 6D strain rate
    de6 = assoc_flow_c(cs6,1.,potential)
    vol = de6[:3].sum()
    de6[:3] = de6[:3] - vol/3.
    return de6

def get_stime():
    """
    date/time format for naming convention in VPSC-FLD
    """
    import time
    t=time.asctime()
    date = time.strftime('%Y%m%d')
    dat,month,dum,time,year = t.split()
    hr,mn,sec = time.split(':')
    return date,hr+mn+sec

def gen_hash_code(nchar=6):
    import hashlib, time
    ## -------------------------------------------------------
    ## Gen HASH code
    m = hashlib.md5()
    m.update(get_stime()[0]+get_stime()[1])
    m.update(time.asctime())
    m.update(str(time.time()))
    return m.hexdigest()[:nchar]
