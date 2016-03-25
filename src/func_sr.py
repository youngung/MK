"""
Strain rate sensitivity functions


All strain rate terms are 'equivalent' scholar values of 
the full tensorial strain rate
"""
import numpy as np

def func_power(ed,ed0,m):
    """
    Power law model

    Kleemola and Ranta-Ekola 1979
    Hosford and Caddell, 1983

    sigma = (ed/ed0)**m
    """
    return (ed/ed0)**m

def func_jc(ed,e0,m):
    """
    Johnson-Cook, 1983
    sigma = (1+m*np.log(ed/ed0))
    """
    return 1+m*np.log(ed/ed0)

def c_F(iopt,**kwargs):
    """
    Return characterized strain rate function
    The returned strain rate functions would be useful 
    when the repeated use of strain rate function is required.

    Returned functions are already characterized by **kwargs.
    Returned functions are now functions of only strain rate (edot_eq)

    Argumemtns
    ==========
    iopt
    edot_eq
    **kwargs
    """
    def func(edot_eq):
        if iopt==0:
            return func_power(edot_eq,**kwargs)
        elif iopt==1:
            return func_jc(edot_eq,**kwargs)
        else:
            raise IOError, 'Unexpected case'
    return func
