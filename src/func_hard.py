"""
"""
def func_swift(eps,k,eps_0,n):
    """
    sigma = k * (eps+eps_0) **n

    Arguments
    ---------
    eps,k,eps_0,n
    """
    return k * (eps+eps_0) **n

def func_swift_sr(eps_dot,eps,k,eps_0,n,m):
    """
    sigma = k * (eps + eps_0)**n * (eps_dot)**m
    """
    return k * (eps + eps_0)**n * (eps_dot)**m
