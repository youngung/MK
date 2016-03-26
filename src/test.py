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
    else:
        raise IOError, 'Could not find the test'

