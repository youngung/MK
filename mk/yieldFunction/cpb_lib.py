import numpy as np

def deviator(a):
    """
    Calculate deviator of 6-component stress tensor
    """
    if type(a).__name__=='ndarray':
        S=a.copy()
    else:
        S=np.array(a,dtype='float')
    # make sure that the given stress <sigma> is a deviator
    p = S[:3].sum()
    S[:3] = S[:3] - p/3.
    return S

def calcC(c11=1.,c22=1.,c33=1.,c44=1.,c55=1.,c66=0.,c23=0.,c13=0.,c12=0.):
    """
    Calculate C (6x6) matrix
    The default value of each individual diagonal components is unity.
    The default value of each individual off-diagonal components is zero.

    Arguments
    ---------
    c11
    c22
    c33
    c44
    c55
    c66
    c23
    c13
    c12
    """
    c  = [[ c11,  c12,  c13,   0.,   0.,   0.],
          [ c12,  c22,  c23,   0.,   0.,   0.],
          [ c13,  c23,  c33,   0.,   0.,   0.],
          [  0.,   0.,   0.,  c44,   0.,   0.],
          [  0.,   0.,   0.,   0.,  c55,   0.],
          [  0.,   0.,   0.,   0.,   0.,  c66]]
    return np.array(c)

def C2comp(C):
    """
    Easy converter to break-down C(6x6) matrix to
    individual components.

    Arguments
    ---------
    C (6x6) matrix
    """
    ## diagonal
    c11=C[0,0]
    c22=C[1,1]
    c33=C[2,2]
    c44=C[3,3]
    c55=C[4,4]
    c66=C[5,5]

    ## off-diagonal
    c23=C[1,2]
    c13=C[0,2]
    c12=C[0,1]

    return c11,c22,c33,c44,c55,c66,c23,c13,c12

def returnT():
    """
    4th-order deviatoric projection that transforms a 2nd-order
    tensor to its devaitor.

    Return
    ------
    T
    """
    T=np.zeros((6,6))
    T[0, :3] =  2., -1., -1.
    T[1, :3] = -1.,  2., -1.
    T[2, :3] = -1., -1.,  2.

    for i in [3,4,5]:
        T[i,i] = 3.
        pass

    T = T /3.
    return T
