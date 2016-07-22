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
        return swift(e=e,ks=ks,n=n,e0=e0,m=m,qq=qq)
    return f_swift

def return_voce(a,b0,c,b1,m,qq):
    """
    Voce requires
    a, b0, c, b1, m (m being the rate sensitivity)

    Arguments
    ---------
    a
    b0
    c
    b1
    m
    qq
    """
    from yf_for import voce
    def f_voce(e):
        return voce(e=e,a=a,b0=b0,c=c,b1=b1,m=m,qq=qq)
    return f_voce

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
    f_hard      = return_swift(n=n,m=m,ks=ks,e0=e0,qq=qq)
    a,b,c,d,e,f = f_hard(0.)

    print 'm:', a,m    # (m)
    print 'sig:',b    # sig
    print 'dsig:',c    # dsig
    print 'dm:',d    # dm
    print 'qq:',e    # (qq)
    print 'n:',f    # (n)
