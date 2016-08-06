"""
Collection of materials
"""
## material library in the form of constitutive models
## with an assumption of 'isotropic' hardening
## Collection of materials
from constitutive import Constitutive
from func_hard_for import return_swift, return_voce
import mk.yieldFunction.yf2
from mk.yieldFunction.yf2 import wrapHill48R, VonMises, wrapYLD

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

def IFsteel_yld2000_case1():
    """
    Return constitutive decription for
    IF steel parameters based on YLD2000-2D
    """
    import mk.yieldFunction.tuneYld2000
    f_hrd = return_swift(n=0.28985,m=5e-2, ks=518.968, e0=0.0007648, qq=1e3)
    f_yld = mk.yieldFunction.tuneYld2000.H48toYld(rv=[2.2,2.0,2.9,1.0],m=6)
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


## 20160722
r_exp    = [2.0112440248052055, 1.8235555903215106, 2.669009844748163]
y_exp    = [1.                , 1.02572944        , 0.99007923]
rv3_vpsc = [1.8451121         , 1.64375413        , 2.27514493]##6k grains     #[1.42453631        , 1.43288609        , 1.94065753]##2k grains
ys3_vpsc = [1.                , 1.02364121        , 1.01127312]##6k grains    # [1.                , 1.01883562        , 1.01582869]##2k grains
rb_vpsc = 1.
yb_vpsc = 1.10682318938

def IFsteel_H48R_20160722():
    f_hrd = return_voce(a=479.00408591,b0=339.71480479,c=7.68395984,b1=70.86783572,m=5e-2,qq=1e3)
    f_yld = mk.yieldFunction.yf2.wrapHill48R(rvs=[r_exp[0],r_exp[1],r_exp[2]])
    return Constitutive(f_yld=f_yld, f_hrd=f_hrd)
def IFsteel_H48YR0_20160722():
    f_hrd = return_voce(a=479.00408591,b0=339.71480479,c=7.68395984,b1=70.86783572,m=5e-2,qq=1e3)
    f_yld = mk.yieldFunction.yf2.wrapHill48Y(ys=[y_exp[0],y_exp[1],y_exp[2]],r0=r_exp[0])
    return Constitutive(f_yld=f_yld, f_hrd=f_hrd)
def IFsteel_H48YR90_20160722():
    f_hrd = return_voce(a=479.00408591,b0=339.71480479,c=7.68395984,b1=70.86783572,m=5e-2,qq=1e3)
    y_yld = mk.yieldFunction.yf2.wrapHill48Y(ys=[y_exp[0],y_exp[1],y_exp[2]],r90=r_exp[2])
    return Constitutive(f_yld=f_yld, f_hrd=f_hrd)
def IFsteel_H48R_vpsc_20160722():
    f_hrd = return_voce(a=479.00408591,b0=339.71480479,c=7.68395984,b1=70.86783572,m=5e-2,qq=1e3)
    f_yld = mk.yieldFunction.yf2.wrapHill48R(rvs=rv3_vpsc)
    return Constitutive(f_yld=f_yld, f_hrd=f_hrd)
def IFsteel_H48YR0_vpsc_20160722():
    f_hrd = return_voce(a=479.00408591,b0=339.71480479,c=7.68395984,b1=70.86783572,m=5e-2,qq=1e3)
    f_yld = mk.yieldFunction.yf2.wrapHill48Y(r0=rv3_vpsc[0], ys=ys3_vpsc)
    return Constitutive(f_yld=f_yld, f_hrd=f_hrd)
def IFsteel_H48YR90_vpsc_20160722():
    f_hrd = return_voce(a=479.00408591,b0=339.71480479,c=7.68395984,b1=70.86783572,m=5e-2,qq=1e3)
    f_yld = mk.yieldFunction.yf2.wrapHill48Y(r90=rv3_vpsc[0], ys=ys3_vpsc)
    return Constitutive(f_yld=f_yld, f_hrd=f_hrd)
def IFsteel_yld1_20160722():
    f_hrd = return_voce(a=479.00408591,b0=339.71480479,c=7.68395984,b1=70.86783572,m=5e-2,qq=1e3)
    f_yld = mk.yieldFunction.tuneYld2000.H48toYld_withYS(rv=[r_exp[0],r_exp[1],r_exp[2]],ys=y_exp,m=6)
    return Constitutive(f_yld=f_yld, f_hrd=f_hrd)
def IFsteel_yld2_20160722():
    f_hrd = return_voce(a=479.00408591,b0=339.71480479,c=7.68395984,b1=70.86783572,m=5e-2,qq=1e3)
    f_yld = mk.yieldFunction.tuneYld2000.wrapYLD(r=[rv3_vpsc[0],rv3_vpsc[1],rv3_vpsc[2],rb],y=[ys3_vpsc[0],ys3_vpsc[1],ys3_vpsc[2],yb_vpsc],m=6)
    return Constitutive(f_yld=f_yld, f_hrd=f_hrd)
## End of 20160722 block

def library(iopt):
    if iopt==0:
        return IsoMat
    elif iopt==1:
        return IFsteel
    elif iopt==2:
        return IFsteel_yld2000_case1
    ## 20160608 starts here
    elif iopt==3:
        return IFsteel_Hill48R_20160608
    elif iopt==4:
        return IFsteel_yld2000_2d_1_20160608
    elif iopt==5:
        return IFsteel_yld2000_2d_2_20160608
    ## 20160722 starts here
    elif iopt==6:
        return IFsteel_H48R_20160722
    elif iopt==7:
        return IFsteel_H48YR0_20160722
    elif iopt==8:
        return IFsteel_H48YR90_20160722
    elif iopt==9:
        return IFsteel_H48R_vpsc_20160722
    elif iopt==10:
        return IFsteel_H48YR0_vpsc_20160722
    elif iopt==11:
        return IFsteel_H48YR90_vpsc_20160722
    elif iopt==12:
        return IFsteel_yld1_20160722
    elif iopt==13:
        return IFsteel_yld2_20160722
    else:
        print 'given iopt:', iopt
        print 'type(iopt).__name__:', type(iopt).__name__
        raise IOError,'not quite ready for this material'
