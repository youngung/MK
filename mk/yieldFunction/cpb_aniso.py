"""
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import mk.library.lib,cpb_data,cpb_iso
import cpb_lib

def eq9(S,a,k):
    """
    Orthotropic criterion:

    (|S_1| - k S_1)**a + (|S_2| - k S_2)**a + (|S_3|-k S_3)**a = F

    Arguments
    ---------
    S
    a
    k

    Returns
    -------
    F
    """
    if len(S)!=3:
        print 'Error, I expect three eigen values'
    F=0.
    for i in xrange(len(S)):
        F=F+(np.abs(S[i]) - k * S[i])**a
    return F

def main(s=[1,0,0,0,0,0],iopt=0,inorm=False):
    """
    Arguments
    ---------
    s = [1,0,0,0,0,0] - Cauchy stress
    iopt (used in cpb_data.tab1_cpb)
    inorm = False (if True, normalize the yield locus by
            tension stress along RD)
    """
    if type(s).__name__=='ndarray':
        S=s.copy()
    else:
        s=np.array(s,dtype='float')
        S=s.copy()
    Sdev = cpb_lib.deviator(S)

    C, k = cpb_data.tab1_cpb(iopt=iopt)
    # C, k = cpb_data.tab2_cpb(iopt=iopt)

    # ## find principal values
    Sig = princ(s, C)
    Sig = np.sort(Sig)

    a = 2.
    F = eq9(S=Sig, a=a, k=k)
    f1 = (F)**(1./a)

    if inorm:
        """
        normarlize w.r.t. uniaxial tension yield stress along RD
        """
        s0t, s0c = eq12(C,a,k)
        f1=f1*s0t
    return S.copy()/f1

def locus(yfunc=main,nth=100,iplot=False,**kwargs):
    """
    Arguments
    ---------
    yfunc
    nth
    iplot
    **kwargs

    Returns
    -------
    X,Y coordinates of polygonal yield locus in plane stress space
    of (s11, s22) with s_ij being cauchy stress compoents.
    """
    x,y = cpb_iso.locus(yfunc,nth,iplot,**kwargs)
    return x, y

def princ2(cauchyStress,C):
    """
    Calculate the principal values using
    Sig = C:T:s form

    Arguments
    ---------
    cauchyStress - the linearly transformed stress Sigma
    C (6x6) tensor for linear transformation

    Returns
    -------
    S1
    S2
    S3
    """
    sqrt=np.sqrt
    s = cauchyStress.copy()
    T = cpb_lib.returnT()
    Sig = np.tensordot(C,np.tensordot(T,s,axes=(1,0)),axes=(1,0))

    Sxx,Syy,Szz,Sxy = Sig[0], Sig[1], Sig[2], Sig[5]

    S1 = 0.5*(Sxx + Syy + sqrt( (Sxx-Syy)**2 + 4*Sxy**2))
    S2 = 0.5*(Sxx + Syy - sqrt( (Sxx-Syy)**2 + 4*Sxy**2))
    S3 = Szz
    return S1, S2, S3

def princ(cauchyStress,C):
    """
    Arguments
    ---------
    cauchyStress - the linearly transformed stress Sigma
    C (6x6) tensor for linear transformation

    Returns
    -------
    S1
    S2
    S3
    """
    sqrt=np.sqrt
    Sxx,Syy,Szz,Sxy = eq11_aux(cauchyStress,C)

    S1 = 0.5*(Sxx + Syy + sqrt( (Sxx-Syy)**2 + 4*Sxy**2))
    S2 = 0.5*(Sxx + Syy - sqrt( (Sxx-Syy)**2 + 4*Sxy**2))
    S3 = Szz
    return S1, S2, S3

def eq11_aux(cauchyStress,C):
    """
    Use the given Cauch stress, calculate
    the componets required to obtain
    principal values of 'Sigma'

    Arguments
    ---------
    cauchyStress
    C (6x6) tensor for linear transformation

    Returns
    -------
    Sxx
    Syy
    Szz
    Sxy
    """
    s=cauchyStress[::]

    # make sure that the stress state is planar.
    if (s[2:5]!=0.).any():
        raise IOError, 'The given state is not plane stress'
    sxx, syy, sxy = s[0], s[1], s[5]

    ## to help write equation
    c11,c22,c33,c44,c55,c66,c23,c13,c12 = cpb_lib.C2comp(C)

    p1=1./3.
    p2=2./3.

    Sxx = (p2*c11 - p1*c12 - p1*c13)*sxx + (-p1*c11 + p2*c12 - p1*c13)*syy
    Syy = (p2*c12 - p1*c22 - p1*c23)*sxx + (-p1*c12 + p2*c22 - p1*c23)*syy
    Szz = (p2*c13 - p1*c23 - p1*c33)*sxx + (-p1*c13 + p2*c23 - p1*c33)*syy
    Sxy = c66*sxy

    return Sxx, Syy, Szz, Sxy

def calcPhis(C):
    """
    Returns Phi1, Phi2, Phi3 for Eq12

    Argument
    --------
    C (6x6) matrix
    """
    p1   = 2./3.
    p2   = 1./3.
    c11,c22,c33,c44,c55,c66,c23,c13,c12 = cpb_lib.C2comp(C)

    Phi1 = (p1*c11-p2*c12-p2*c13)
    Phi2 = (p1*c12-p2*c22-p2*c23)
    Phi3 = (p1*c13-p2*c23-p2*c33)
    return Phi1, Phi2, Phi3

def eq12(C,a,k):
    """
    Calculate tension and compression yield stress
    at phi = 0.

    Arguments
    ---------
    C (6x6)
    a : exponent
    k : k parameter
    """
    p1,p2,p3 = calcPhis(C)
    Phis = np.array([p1,p2,p3])
    s0t_ = 0.
    s0c_ = 0.
    for i in xrange(len(Phis)):
        s0t_ = s0t_ + (np.abs(Phis[i])-k*Phis[i])**a
        s0c_ = s0c_ + (np.abs(Phis[i])+k*Phis[i])**a
    s0t = (1./s0t_) ** (1./a)
    s0c = (1./s0c_) ** (1./a)
    return s0t, s0c

def calcPsis(C):
    """
    Returns Psi1, Psi2, Psi3 for eq13

    Argument
    --------
    C (6x6) matrix
    """
    p1   = 1./3.
    p2   = 2./3.
    c11,c22,c33,c44,c55,c66,c23,c13,c12 = cpb_lib.C2comp(C)

    Psi1 = -p1*c11+p2*c12-p1*c13
    Psi2 = -p1*c12+p2*c22-p1*c23
    Psi3 = -p1*c13+p2*c23-p1*c33
    return Psi1, Psi2, Psi3

def eq13(C,a,k):
    """
    Calculate tension and compression yield stress
    at phi = 0.

    Arguments
    ---------
    C (6x6)
    a : exponent
    k : k parameter
    """
    p1,p2,p3 = calcPsis(C)
    Psis = np.array([p1,p2,p3])
    s90t_=0
    s90c_=0.
    for i in xrange(len(Psis)):
        s90t_ = s90t_ + (np.abs(Psis[i])-k*Psis[i])**a
        s90c_ = s90c_ + (np.abs(Psis[i])+k*Psis[i])**a
    s90t = (1./s90t_)**(1./a)
    s90c = (1./s90c_)**(1./a)
    return s90t, s90c
