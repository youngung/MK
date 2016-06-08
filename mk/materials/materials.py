## material library in the form of constitutive models
## with an assumption of 'isotropic' hardening
## Collection of materials
from constitutive import Constitutive
from func_hard_for import return_swift
from mk.yieldFunction.yf2 import wrapHill48R, VonMises, wrapYLD

def IFsteel_yld2000_case1():
    """
    Return constitutive decription for
    IF steel parameters based on YLD2000-2D
    """
    f_hrd = return_swift(n=0.28985,m=5e-2, ks=518.968, e0=0.0007648, qq=1e3)
    f_yld = tuneYld2000.H48toYld(rv=[2.2,2.0,2.9,1.0],m=6)
    return Constitutive(f_yld=f_yld, f_hrd=f_hrd)

def IFsteel():
    """
    Return constitutive description for
    IF steel parameters
    """
    ## hardening model
    f_hrd = return_swift(n=0.28985,m=5e-2, ks=518.968, e0=0.0007648, qq=1e3)
    ## yield function characterized by three r-values
    f_yld = wrapHill48R([2.2, 2.0, 2.9]) # r0, r45, r90
    return Constitutive(f_yld=f_yld, f_hrd=f_hrd)


def IFsteel_Hill48R_20160608():
    """
    IF steel parameters tuned for Hill48R

    Hardening parameters were tuned for bulge test data
    """
    ## hardening model
    f_hrd = return_swift(n=0.255392940,m=5e-2,ks= 6.00162331e+02, e0=4.23405410e-04, qq=1e3)
    f_yld = wrapHill48R([2.092113652699876, 1.8999691982523326, 2.8779441147053473])
    return Constitutive(f_yld=f_yld,f_hrd=f_hrd)


def IFsteel_yld2000_2d_1_20160608():
    """
    IF steel parameters tuned for yld2000_2d
    -- yld2000-2d model tuned by experimental r0, r45, r90.
    rb, y0, y45, y90, and yb were tuned using Hill48R(r0,r45,r90)

    Hardening parameters were tuned for bulge test data
    """
    ## hardening model
    f_hrd = return_swift(n=0.255392940,m=5e-2,ks= 6.00162331e+02, e0=4.23405410e-04, qq=1e3)
    import mk.yieldFunction.tuneYld2000
    f_yld = mk.yieldFunction.tuneYld2000.H48toYld(
        rv=[2.092113652699876, 1.8999691982523326, 2.8779441147053473],m=6)
    return Constitutive(f_yld=f_yld,f_hrd=f_hrd)


def IFsteel_yld2000_2d_2_20160608():
    """
    IF steel parameters tuned for yld2000_2d
    -- yld2000-2d model tuned by experimental r0, r45, r90, y0, y45, y90
    rb yb were tuned using Hill48R(r0,r45,r90)

    Hardening parameters were tuned for bulge test data
    """
    ## hardening model
    f_hrd = return_swift(n=0.255392940,m=5e-2,ks= 6.00162331e+02, e0=4.23405410e-04, qq=1e3)
    import mk.yieldFunction.tuneYld2000
    f_yld = mk.yieldFunction.tuneYld2000.H48toYld_withYS(
        rv=[2.092113652699876, 1.8999691982523326, 2.8779441147053473],
        ys=[1.        ,  1.02781532,  0.98441769],m=6)
    return Constitutive(f_yld=f_yld,f_hrd=f_hrd)

def IsoMat():
    """
    Return constitutive description for
    Isotropic von Mises material

    Hardening parameters were tuned for bulge test data (of the IF steel)
    """
    ## hardening model
    f_hrd = return_swift(n=0.255392940,m=5e-2,ks= 6.00162331e+02, e0=4.23405410e-04, qq=1e3)
    ## yield function
    f_yld = VonMises
    return Constitutive(f_yld=f_yld, f_hrd=f_hrd)

def library(iopt):
    if iopt==0:
        return IsoMat
    elif iopt==1:
        return IFsteel
    elif iopt==2:
        return IFsteel_yld2000_case1
    elif iopt==3:
        return IFsteel_Hill48R_20160608
    elif iopt==4:
        return IFsteel_yld2000_2d_1_20160608
    elif iopt==5:
        return IFsteel_yld2000_2d_2_20160608
    else:
        raise IOError,'not quite ready for this material'
