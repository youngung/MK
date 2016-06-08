# ### characterize hardening functions
import numpy as np
def char_objf(f_hard,xs,ys):
    """
    objf1:
    f_hard should be a function of only args
    -- It is assumed that strain <epsilon>
    are characterized...

    objf2:
    f_hard should be a function of only strain.
    """
    def objf1(args):
        yf = f_hard(xs,args)
        diff = np.sqrt(((yf - ys)**2).sum())/(len(xs)-1)
        return diff
    def objf2(args):
        def f(xs):
            print args
            return f_hard(xs,args)
        return f
    return objf1, objf2

def main(exp_dat,f_hard,x0):
    """
    Given the strain hardening function *wrapper*,
    identify the proper hardening parameters by fitting
    with the experiment data.

    Once the parameters are identified, return
    the characterized strain hardening function that is
    returned by the wrapper.

    Arguments
    ---------
    exp_dat  :  exp_dat experimental data in the form of array
    f_hard   :  strain hardening function
    x0
    """
    from scipy.optimize import minimize
    x,y = exp_dat
    xs = np.linspace(x[0],x[1],100)
    ys = np.interp(xs,x,y)
    objf1, objf2 = char_objf(f_hard=f_hard,xs=xs,ys=ys)
    res = minimize(objf1,x0=x0,method='Powell')
    func_hard = objf2(res.x)
    return func_hard

def charFunc(func,xs,ys,popt_guess):
    """
    Arguments
    ---------
    func
    xs
    ys
    popt_guess
    """
    exp_dat = (xs,ys)
    f_char_hard = main(
        exp_dat = exp_dat,
        f_hard = func, x0=popt_guess)
    return f_char_hard

def test1():
    import numpy as np
    from func_hard import func_swift

    popt_guess = (518.968, 0.0007648, 0.28985) ## ks, e0, n
    x=np.linspace(0,0.2,1000)
    y=func_swift(x,popt_guess)
    exp_dat= (x,y)

    f_hard = main(
        exp_dat=exp_dat,
        f_hard = func_swift,x0=popt_guess)

    return f_hard

if __name__=='__main__':
    test1()
