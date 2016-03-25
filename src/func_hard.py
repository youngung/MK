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
    eps,k,eps_0,n
    """
    return k * (eps+eps_0) **n


def func_hollomon(eps,k,n):
    """
    Hollomon 1945
    sigma = K * eps **n

    Arguments
    ---------
    eps,k,n
    """
    return k * eps**n

def func_ludwik(eps,k,sig0,n):
    """
    Ludwik, 1909
    sigma = sig0 + k*eps**n

    Arguments
    ---------
    eps,k,sig0,n
    """
    return sig0 + k*eps**n


def func_hs(eps,sig0,k,eps0,n):
    """
    Hartley and Srinivasa, 1983

    sigma = sig0 + K (eps0 + eps) **n

    Arguments
    ---------
    eps,sig0,k,eps0,n
    """
    return sig0 + K * (eps0 + eps)**n

def func_ludwigson(eps,k1,n1,k2,n2,):
    """
    Ludwigson, 1971

    sigma = k1 * eps**n1 + exp(k2+n2*eps)
    """
    return k1*eps**n1+np.exp(k2+n2*eps)

def func_baragar(eps,sig0,c,d,e):
    """
    Baragar, 1987

    sig = sig0 + c*eps**0.4 + d*eps**0.8+e*eps**1.2
    """
    return sig0 + c*eps**0.4+ d*eps**0.8+e*eps**1.2


# def func_swift_sr(eps_dot,eps,k,eps_0,n,m):
#     """
#     sigma = k * (eps + eps_0)**n * (eps_dot)**m
#     """
#     return k * (eps + eps_0)**n * (eps_dot)**m
