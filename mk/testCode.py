"""
Various tests
"""

def test1():
    print 'Removed feature'
    # import fldb
    # fldb.test_save_a()

def test2():
    print 'Removed feature'
    # import fldb
    # fldb.main()

def test3():
    import numpy as np
    import matplotlib as mpl
    mpl.use('agg')
    import matplotlib.pyplot as plt
    import mk.tests.mechtests, mk.yieldFunction.yf2
    fyld = mk.yieldFunction.yf2.wrapHill48R([1.,2,1.5])
    psis, rvs, phis = mk.tests.mechtests.inplaneTension(fyld)

    ##
    r2d = np.pi/180.
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    ax1.plot(psis*r2d, rvs)
    fn='testCodeTest3.pdf'
    print '%s has been saved'%fn
    fig.savefig(fn)

def test4():
    """
    MK module compute time benchmark
    for the below condition:
    f0    = 0.995 (initial inhomogeneity factor)
    psi0  = 0     (orientation of the groove)
    th    = 0 (plane-strain tension RD condition)
    """
    import mk, time
    ## isotropic material that follows von Mises constitutive law
    from materials import IsoMat
    print 'TEST4 in test.py for estimating the compute time'
    from MP import progress_bar
    uet = progress_bar.update_elapsed_time
    t0 = time.time()
    mk.main(f0=0.995, psi0=0, th=0,logFileName='/tmp/mk-test4.log')
    dt = time.time() - t0
    uet(dt,'Total elapsed time for benchmark (test4)');print

def test5():
    import numpy as np
    import mk
    from lib import rho2th
    th = rho2th(-0.6)*180./np.pi
    mk.main(f0=0.996,psi0=3, th=th,logFileName='/tmp/dum.log')

def test6():
    """
    a short FLD calcultion
    """
    cmd = 'python mk_run.py --f0 0.995 --r0 -0.5 --r1 1. --nr 4 --mat 0'
    import os
    print 'cmd:',cmd
    os.system(cmd)

if __name__=='__main__':
    import argparse
    ## Arguments parsing
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', type=int, help='test number')

    args = parser.parse_args()
    ind = args.i

    if ind==1:
        test1()
    elif ind==2:
        test2()
    elif ind==3:
        test3()
    elif ind==4:
        test4()
    elif ind==5:
        test5()
    elif ind==6: ## FLD actual calculations
        test6()
    else:
        raise IOError, 'Could not find the test'
