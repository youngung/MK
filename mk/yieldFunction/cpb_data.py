"""
Collection of data
"""
import cpb_lib

def tab1_cpb(iopt):
    """
    Mg-Th coefficients
         k         c12    c13    c22    c23    c33
    1%   0.3539 0.4802 0.2592 0.9517 0.2071 0.4654
    5%   0.2763 0.3750 0.0858 0.9894 0.0659 0.1238
    10%  0.0598 0.6336 0.2332 1.4018 0.5614 0.7484

    iopt=0: returns  1%
    iopt=1: returns  5%
    iopt=2: returns 10%
    """
    k    = 0.
    c11  = 0.
    c12  = 0.
    c13  = 0.
    c22  = 0.
    c23  = 0.
    c33  = 0.
    c44  = 0.
    c55  = 0.
    c66  = 0.

    if iopt==0:
        k,c12,c13,c22,c23,c33=0.3539,0.4802,0.2592,0.9517,0.2071,0.4654
    elif iopt==1:
        k,c12,c13,c22,c23,c33=0.2763,0.3750,0.0858,0.9894,0.0659,0.1238
    elif iopt==2:
        k,c12,c13,c22,c23,c33=0.0598,0.6336,0.2332,1.4018,0.5614,0.7484
    else:
        raise IOError, 'Unexpected case.'
    c = cpb_lib.calcC(c12=c12,c13=c13,c22=c22,c23=c23,c33=c33)
    return c, k

def tab2_cpb(iopt):
    """
    Mg-Li coefficients
         k         c12    c13    c22    c23    c33
    1%   0.2026 0.5871 0.6975 0.9783 0.2840 0.1497
    5%   0.2982 0.6103 0.8056 1.0940 0.5745 0.1764
    10%  0.1763 0.5324 0.8602 1.0437 0.8404 0.2946

    iopt=0: returns  1%
    iopt=1: returns  5%
    iopt=2: returns 10%
    """
    k    = 0.
    c11  = 0.
    c12  = 0.
    c13  = 0.
    c22  = 0.
    c23  = 0.
    c33  = 0.
    c44  = 0.
    c55  = 0.
    c66  = 0.

    if iopt==0:
        k,c12,c13,c22,c23,c33=0.2026,0.5871,0.6975,0.9783,0.2840,0.1497
    elif iopt==1:
        k,c12,c13,c22,c23,c33=0.2982,0.6103,0.8056,1.0940,0.5745,0.1764
    elif iopt==2:
        k,c12,c13,c22,c23,c33=0.1763,0.5324,0.8602,1.0437,0.8404,0.2946
    else:
        raise IOError, 'Unexpected case.'
    c = cpb_lib.calcC(c12=c12,c13=c13,c22=c22,c23=c23,c33=c33)
    return c, k
