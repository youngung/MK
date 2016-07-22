"""
Strain hardening functions


All stress/strain terms are 'equivalent' scholar values of the
respective tensors.
"""
import numpy as np
def func_swift(eps,k,eps_0,n):
    """
    sigma = k * (eps+eps_0) **n

    Arguments
    ---------
    eps
    k
    eps_0
    n

    Returns
    -------
    k*(eps+eps_0)**n
    """
    return k * (eps+eps_0) **n

def func_voce(eps,a,b0,c,b1):
    """
    sigma = a-b0*exp(-c*eps) + b1 *eps

    Arguments
    ---------
    eps
    a
    b0
    c
    b1

    Returns
    -------
    a-b0*exp(-c*eps) + b1 *eps
    """
    return a-b0*np.exp(-c*eps) + b1 *eps

# def func_swift_d(eps,k,eps_0,n):
#     """
#      sigma  =     k * (eps+eps_0) ** n
#     dsigma  = n * k * (eps+eps_0) ** (n-1)

#     Arguments
#     ---------
#     eps,k,eps_0,n
#     """
#     return k * (eps+eps_0) **n

def func_hollomon(eps,k,n):
    """
    Hollomon 1945
    sigma = K * eps **n

    Arguments
    ---------
    eps,k,n

    Returns
    -------
    k * eps**n
    """
    return k * eps**n

def func_ludwik(eps,k,sig0,n):
    """
    Ludwik, 1909
    sigma = sig0 + k*eps**n

    Arguments
    ---------
    eps,k,sig0,n

    Returns
    -------
    sig0 + k*eps**n
    """
    return sig0 + k*eps**n

def func_hs(eps,sig0,k,eps0,n):
    """
    Hartley and Srinivasa, 1983

    sigma = sig0 + K (eps0 + eps) **n

    Arguments
    ---------
    eps,sig0,k,eps0,n

    Returns
    -------
    sig0 + K * (eps0 + eps)**n
    """
    return sig0 + K * (eps0 + eps)**n

def func_ludwigson(eps,k1,n1,k2,n2,):
    """
    Ludwigson, 1971

    sigma = k1 * eps**n1 + exp(k2+n2*eps)

    Arguments
    --------
    eps
    k1
    n1
    k2
    n2

    Returns
    -------
    k1*eps**n1+np.exp(k2+n2*eps)
    """
    return k1*eps**n1+np.exp(k2+n2*eps)

def func_baragar(eps,sig0,c,d,e):
    """
    Baragar, 1987

    sig = sig0 + c*eps**0.4 + d*eps**0.8+e*eps**1.2
    """
    return sig0 + c*eps**0.4+ d*eps**0.8+e*eps**1.2


## place-holder for flow stress that accounts for strain-rate
## sensitivity
def flowStress(f_hard,eps,d,q,a):
    """
    Flow stress as function of strain-hardening curve
    multiplied by a factor that accounts for strain rate
    sensitivity. Here, the strain rate sensitivity
    is described by a power law such that


    sigma = (f).(a)^q

    with f is the strain-hardening curve and (a)^q describes
    the strain rate dependence.

    Arguments
    ---------
    f_hard
    eps
    d
    q
    a
    """

    pass


def c_G(iopt,**kwargs):
    """
    Return characterized hardening function
    The returned hardening functions would be useful
    when the repeated use of hardening function is required.

    Returned functions are already characterized by **kwargs.
    Returned functions are now functions of only equivalent strain (eps_eq)

    Argumemtns
    ==========
    iopt
    **kwargs
    """
    def func(eps_eq):
        if iopt==0:
            return func_swift(eps_eq,**kwargs)
        elif iopt==1:
            return func_hollomon(eps_eq,**kwargs)
        elif iopt==2:
            return func_ludwik(eps_eq,**kwargs)
        elif iopt==3:
            return func_hs(eps_eq,**kwargs)
        elif iopt==4:
            return func_ludwigson(eps_eq,**kwargs)
        elif iopt==5:
            return func_baragar(eps_eq,**kwargs)

    # def func_d(eps_eq):
    #     if iopt==0:
    #         return func_swift_d(eps_eq,**kwargs)
    #     elif iopt==1:
    #         return np.nan
    #     elif iopt==2:
    #         return np.nan
    #     elif iopt==3:
    #         return np.nan
    #     elif iopt==4:
    #         return np.nan
    #     elif iopt==5:
    #         return np.nan

    return func, func_d
