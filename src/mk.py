from for_lib import vm
import numpy as np
cos=np.cos
sin=np.sin
tan=np.tan
log=np.log
atan2=np.arctan2
sqrt=np.sqrt

def main(iverbose=0):
    ## given rho:
    rho = -0.6
    ## reference strain vector (pth) of rho
    pth = np.zeros(6)
    pth[0] = 1.
    pth[1] = rho

    ## Find the corresponding stress direction on the yield locus while maintaining (phi=1.)
    ## corresponding two locations of stress that embraces the presumable stress that gives the rho. - correct stress will be found between these two stresses
    s1=np.zeros(6); s2=np.zeros(6)
    
    ## s1 should be far left from uniaxial; s2 should be far right from uniaxial
    s1[0]=0; s1[1]=-1
    s2[0]=1; s2[1]= 1

    ## find the correct stress that givesn rho of -0.6
    if iverbose>=3: print (8*'%6s ')%('phi','s1','s2','sa1','sa2','sb1','sb2','diff')
    diff = np.sqrt(((s1-s2)**2).sum())
    it=0
    while diff>1e-10:
        s       = (s1[::]+s2[::])/2.
        phi,dphi,d2phi = vm(s)
        rac     = sqrt(dphi[0]*dphi[0] + dphi[1]*dphi[1])
        dphi[0] = dphi[0]/rac
        dphi[1] = dphi[1]/rac

        ## narrow down s1-s2 bounds
        if dphi[0]*pth[1]-dphi[1]*pth[0]>=0:
            s1[:] = s[:]
        else:
            s2[:] = s[:]

        it=it+1
        if (it>100): raise IOError, 'Could not find the proper s'
        diff = np.sqrt(((s1-s2)**2).sum())
        if iverbose>=3: print (7*'%6.3f '+'%9.3e')%(phi,s[0],s[1],s1[0],s1[1],s2[0],s2[1],diff)

    rho_ = dphi[1]/dphi[0]
    alf_ = s[1]/s[0]
    if iverbose>=3: print 'it:',it
    if iverbose>=3: print 'rho, rho_, alf_',rho,rho_,alf_

    onepath(sa=s,psi0=0.,f0=0.999)

def return_swift(m,ks,e0,qq):
    from for_lib import swift
    def f_swift(e):
        return swift(e,ks,m,e0,qq)
    return f_swift

def test():
    m=1.
    ks=2.
    e0=3.
    qq=4.
    f_hard=return_swift(m,ks,e0,qq)
    a,b,c,d,e=f_hard(0.)
    print a    # (m)
    print b    # sig
    print c    # dsig
    print d    # dm
    print e    # (qq)

def onepath(sa,psi0,f0):
    from lib import rot_6d
    # from for_lib import swift
    phia,fa,f2a=vm(sa)

    # ## debug
    # sa=np.array([1,2,0,0,0,0])
    # psi0=15.*np.pi/180.

    sx=rot_6d(sa,-psi0)
    print 'sa:,',sa
    print 'xa0,ya0,za0'
    print (3*'%6.3f ')%(sx[0],sx[1],sx[5])
    # return

    ks = 500
    ma  = 5e-1
    e0 = 1e-5
    qq = 1000. ##Strain rate ratio (E_a.dot / E_0.dot)
    f_hard = return_swift(ma,ks,e0,qq)

    print ks,ma,e0



    ma,siga,dsiga,dma,qq = f_hard(0.)
    print 'siga,dsiga,ma,dma'
    print siga,dsiga,ma,dma


    ## initial_conditions
    dydx = np.zeros(500)
    dydx[0] = 1.
    ndim = 4
    ncase=1
    b=np.zeros(20)
    b[0]= psi0
    b[1]= f0
    b[2]= siga
    b[3]= ma
    b[4]= qq
    b[5]= sx[5]/sx[0]

    xzero = np.array([1,1,0,0])


    for i in xrange(6):
        print 'B%i'%(i+1),'%7.2f'%b[i]
        pass

    xfinal=new_raph_fld(ndim=ndim,ncase=ncase,xzero=xzero,b=b,f_hard=f_hard)

    print 'xfinal:', xfinal
    return

def new_raph_fld(ndim,ncase,xzero,b,f_hard):
    from for_lib import gauss,norme
    xn=np.zeros(ndim)
    xn=xzero[::]
    residu = 1.0
    it = 0
    itmax=100
    eps=1e-10

    while (residu>eps and it<itmax):
        it = it+1
        F, J = fonc_fld(ndim,ncase,b,xn,f_hard)
        res = gauss(ndim,J,F)
        print res
        xn1=xn-res
        residu = norme(ndim,res)
        xn=xn1[::]
        print it,xn
        pass


    if it>=itmax: raise IOError

    return xn1



def fonc_fld(ndim,ncase,b,x,f_hard):
    from lib import rot_6d

    psi0, f0, siga, ma, qq, b6 = b[:6]
    # psi0=15*np.pi/180.
    # x[0]=1.
    # x[1]=2.
    # x[3]=0.

    s11,s22,s12 = x[:3]
    s = np.array([s11,s22,0,0,0,s12])
    sb = rot_6d(s, psi0)

    print 'psi0:',psi0
    print 'x:'
    print x[:3]
    print 'sb:'
    print sb

    phib,fb,f2b=vm(sb)
    print 'phib:',phib
    print 'f2b:'
    for i in xrange(6):
        for j in xrange(6):
            print '%6.2f'%f2b[i,j],
        print

    b[6]=x[3]*fb[0]
    b[7]=x[3]*fb[1]

    db = rot_6d(fb,-psi0)
    print 'db:',db

    f2xb = np.zeros((6,6))
    cp = cos(psi0)
    sp = sin(psi0)
    c2 = cp*cp
    s2 = sp*sp
    sc = sp*cp
    f2xb[0,0] = f2b[0,0]*c2 + f2b[0,1]*s2 + f2b[0,5]*sc
    f2xb[1,0] = f2b[1,0]*c2 + f2b[1,1]*s2 + f2b[1,5]*sc
    f2xb[5,0] =(f2b[5,0]*c2 + f2b[5,1]*s2 + f2b[5,5]*sc)/2
    f2xb[0,1] = f2b[0,0]*s2 + f2b[0,1]*c2 - f2b[0,5]*sc
    f2xb[1,1] = f2b[1,0]*s2 + f2b[1,1]*c2 - f2b[1,5]*sc
    f2xb[5,1] =(f2b[5,0]*s2 + f2b[5,1]*c2 - f2b[5,5]*sc)/2
    f2xb[0,5] = 2*(f2b[0,1]-f2b[0,0])*sc + f2b[0,5]*(c2-s2)
    f2xb[1,5] = 2*(f2b[1,1]-f2b[1,0])*sc + f2b[1,5]*(c2-s2)
    f2xb[5,5] =(2*(f2b[5,1]-f2b[5,0])*sc + f2b[5,5]*(c2-s2))/2
    
    for i in xrange(6):
        for j in xrange(6):
            print '%6.2f'%f2xb[i,j],
        print

    mb,sigb,dsigb,dmb,qq=f_hard(x[3])
    print 'x(4):',x[3]
    f=np.zeros(4)
    f[0] = db[1]
    f[1] = x[2] - x[0]*b[5]
    f[2] = phib - 1.
    f[3] =-log(b[1]*sigb/b[2])+(b[3]-mb)*log(b[4])-x[3]*db[0]
    print '--'
    print 'f:',f
    # print '-log(b[1]*sigb/b[2])',-log(b[1]*sigb/b[2])
    # raise IOError
    # print '(b[3]-mb)*log(b[4])-x[3]*db[0]',(b[3]-mb)*log(b[4])-x[3]*db[0]


    J=np.zeros((4,4))

    J      = np.zeros((4,4))
    J[0,0] = s2  * f2xb[0,0] + c2 * f2xb[1,0] - 2*sc*f2xb[5,0]
    J[0,1] = s2  * f2xb[0,1] + c2 * f2xb[1,1] - 2*sc*f2xb[5,1]
    J[0,2] = s2  * f2xb[0,5] + c2 * f2xb[1,5] - 2*sc*f2xb[5,5]
    J[1,0] = -b[5]
    J[1,2] = 1.
    J[2,0] = db[0]
    J[2,1] = db[1]
    J[2,2] = 2*db[5]
    J[3,3] = -dsigb / sigb - dmb*log(b[4]) - db[0]
    for i in xrange(4):
        for j in xrange(4):
            print '%12.2f '%J[i,j],
        print
    return f, J

if __name__=='__main__':
    main(iverbose=3)
