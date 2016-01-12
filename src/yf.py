"""
Yield functions
"""
import numpy as np
from scipy import optimize
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import os


##
import lib;reload(lib)
from lib import *
##
from MP.mat import voigt
ijv=voigt.ijv
vij=voigt.vij


## load lib_dat for template
p_home = os.getcwd()
p_load = '/Users/yj/repo/vpsc/vpsc-dev-yld/src/py_pack/yld_hah/'
os.chdir(p_load)
import lib_dat
os.chdir(p_home)
##----------------------------------------

def QuadHill(cauchy6,R0,R90):
    """
    Arguments
    ---------
    Cauchy6: 6D Cauchy stress
    R0     : R-value along 0.
    R90    : R-value along 90.
    """
    ## Principal stress space...
    w, v = c6p(cauchy6)
    s1, s2, s3 = w
    phi = s1**2 + R0 * (1+R90) / (R90*(1+R0))*s2**2 -\
          (2*R0)/(1+R0)*s1*s2
    return phi**(1./2.)

def VonMises(cauchy6):
    """
    Von Mises yield locus general

    Arguments
    ---------
    Cauchy6: 6D Cauchy stress

    Returns phi^VM
    ----------
    cauchy6 <6 dimensional>
    """
    s11,s22,s33,s31,s23,s12\
        =cauchy6[0],cauchy6[1],cauchy6[2],\
        cauchy6[3],cauchy6[4],cauchy6[5]
    y = (s11-s22)**2.+(s22-s33)**2.+(s33-s11)**2.
    y = y + 6.*(s12**2+s23**2+s31**2)
    return np.sqrt(0.5*y)

def Hosford(cauchy6,a):
    """
    Arguments
    ---------
    Cauchy6: 6D Cauchy stress
    a      : exponent
    """
    w,v = c6p(cauchy6)
    s1,s2,s3 = w
    phi=0.5*(np.abs(s1-s2)**a + np.abs(s1-s3)**a + np.abs(s2-s3)**a)
    return phi**(1./a)

def illu():
    """
    To illustrate the difference between yield functions
    """
    nths=100
    yld_funcs=[VonMises,QuadHill,Hosford]
    nyld = len(yld_funcs)

    fig = plt.figure(figsize=(3.5*2,3.))
    ax1=fig.add_subplot(121)
    ax2=fig.add_subplot(122)


    yfunc = VonMises
    loc1, loc2 = y_locus(nths,yfunc)
    ax1.plot(loc1[0],loc1[1],label='Von Mises')
    ax2.plot(loc2[0],loc2[1])

    yfunc = QuadHill
    loc1, loc2 = y_locus(nths,yfunc,R0=1.5,R90=2.)
    ax1.plot(loc1[0],loc1[1],label='Quad Hill')
    ax2.plot(loc2[0],loc2[1])

    yfunc = Hosford
    loc1, loc2 = y_locus(nths,yfunc,a=6)
    ax1.plot(loc1[0],loc1[1],label='Hosford')
    ax2.plot(loc2[0],loc2[1])

    ax1.legend()
    ys_temp(ax1)
    lib_dat.pi_rad(ax2,rot=150,rel=1.5)
    fig.savefig('ys_illu.pdf')
