## f2py modules for hardening functions
def return_swift(n,m,ks,e0,qq):
    """
    Swift requires
    ks, e0, n, m (m being the rate sensitivity)

    Arguments
    ---------
    n
    m
    ks
    e0
    qq
    """
    from yf_for import swift
    def f_swift(e):
        return swift(e,ks,n,e0,m,qq)
    return f_swift


def test():
    """
    test various functions.
    """
    ## test swift function
    n=0.1
    m=1.
    ks=2.
    e0=3.
    qq=4.
    f_hard      = return_swift(n,m,ks,e0,qq)
    a,b,c,d,e,f = f_hard(0.)
    print a    # (n)
    print b    # (m)
    print c    # sig
    print d    # dsig
    print e    # dm
    print f    # (qq)
