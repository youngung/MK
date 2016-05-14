"""
Various tests
"""

def test1():
    import fldb
    fldb.test_save_a()

def test2():
    import fldb
    fldb.main()

def test3():
    import materials
    materials.tension_tests()

def test4():
    """
    MK module compute time benchmakr
    """
    import mk, time
    print 'TEST4 in test.py for estimating the compute time'
    from MP import progress_bar
    uet=progress_bar.update_elapsed_time

    t0 = time.time()
    mk.main(f0=0.995, psi0=0, th=0,logFileName='/tmp/mk-test4.log')
    dt = time.time() - t0

    uet(dt,'Total elapsed time');print

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
    else:
        raise IOError, 'Could not find the test'

