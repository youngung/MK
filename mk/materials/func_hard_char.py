# ### characterize hardening functions
import numpy as np
from scipy.optimize import curve_fit

def wrapper(func,*args):
    """
    Hardening function wrapper

    Arguments
    ---------
    func
    *args

    Returns
    -------
    func(x,*args) that is a function of only strain (x).
    """
    def f_hard_char(x):
        """
        Argument
        --------
        x
        """
        return func(x,*args)
    return f_hard_char

def main(exp_dat,f_hard,params):
    """
    Arguments
    ---------
    exp_dat
    f_hard
    params (initial guess)
    """
    x,y = exp_dat
    # bounds --
    # print 'params:', params
    popt, pcov = curve_fit(f_hard,x,y,p0=params)

    return wrapper(f_hard,*popt), popt, pcov

def test1():
    from func_hard import func_swift
    popt_guess = (518.968, 0.0007648, 0.28985) ## ks, e0, n
    x=np.linspace(0,0.2,1000)
    y=func_swift(x,*popt_guess)
    exp_dat= (x,y)
    func = main(exp_dat, func_swift, popt_guess)

if __name__=='__main__':
    test1()
