# ### characterize hardening functions
import numpy as np
from scipy.optimize import curve_fit

def main(exp_dat,f_hard,params):
    """
    Arguments
    ---------
    exp_dat
    f_hard
    x0
    """
    x,y = exp_dat
    # bounds --
    popt, pcov = curve_fit(f_hard,x,y,p0=params)
    def wrapper(func,*args):
        def f_hard_char(x):
            return func(x,*args)
        return f_hard_char
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
