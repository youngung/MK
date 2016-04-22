from fldb import *
import modules; reload(modules)
from scipy import optimize
import os, glob
from MP.mat import mech
from MP import progress_bar
import lib
uet=progress_bar.update_elapsed_time
sin = np.sin
cos = np.cos

class Path:
    """
    Collection of loading history

    s: 6D sigma (stress)
    e: 6D strain accumulative strain
    d: 6D strian rate
    t     time stamps
    """
    def __init__(self,s,e,d,t,dt,ebar,sbar,debar):
        self.s=s.T
        self.e=e.T
        self.d=d.T
        self.t=t
        self.dt=dt
        self.ebar=ebar
        self.sbar=sbar
        self.debar=debar
        self.nstp = len(t)

    def values_mk(self,istp,psi):
        """
        Arguments
        ---------
        istp
        yfunc
        """
        s = self.s[:,istp]
        e = self.e[:,istp]
        d = self.d[:,istp]
        debar = self.debar[istp]
        sbar  = self.sbar[istp]

        self.alpha = s[1]/s[0]
        self.rho   = d[1]/d[0]
        self.beta  = s[5]/s[0]
        self.gamma = d[5]/d[0]
        self.eta   = s[0] / sbar
        self.k1    = 1. + self.alpha * self.rho + 2.*self.beta*self.gamma
        self.k2    = func_k2(psi,self.rho,  self.gamma)
        self.k3    = func_k3(psi,self.alpha,self.beta)

    # def rot_istp(self,istp,psi):
    #     """
    #     Apply coorindate transformation for tensors
    #     (The transformation tensor is defined
    #     by in-plane rotation quantified by psi angle)

    #     Arguments
    #     ---------
    #     istp
    #     psi
    #     """
    #     s = lib.rot_6d(self.s[:,istp],psi)
    #     e = lib.rot_6d(self.e[:,istp],psi)
    #     d = lib.rot_6d(self.d[:,istp],psi)
    #     return s,e,d

class Region:
    def __init__(self,mat):
        ## Yield functions and the derivatives
        self.func_yd=mat.func_yd
        self.func_yd1=mat.func_yd1
        self.func_yd2=mat.func_yd2

        self.func_hd  =mat.func_hd
        self.func_hd_d=mat.func_hd_d

        pass
    def assign(self,s,d,e):
        self.s=s; self.d=d; self.e=e
        self.calc_etc()
        pass
    def calc_etc(self):
        self.wrate = np.dot(self.s,self.d)
        self.s_ = self.func_yd(self.s)
        self.d_ = self.wrate/self.s_
        self.w  = func_hd(self.e) ## flow stress (hardening)
    def calc_s_grv(self,psi):
        self.s_grv=lib.rot_6d(self.s,psi)
        self.b6 = self.s_grv[5]/self.s_grv[0]
        pass

def calc_grv(a,psi):
    """
    """
    return lib.rot_6d(a,psi)

class PathCollection:
    def __init__(self):
        self.paths=[]
        pass
    def add_path(self,p):
        self.paths.append(p)
        self.npath = len(self.paths)
    def draw(self):
        import matplotlib.pyplot as plt
        fig=plt.figure(figsize=(9,7))
        ax1=fig.add_subplot(221)
        ax2=fig.add_subplot(222)
        ax3=fig.add_subplot(223)
        ax4=fig.add_subplot(224)

        for i in xrange(self.npath):
            p = self.paths[i]
            ax1.plot(p.e[1],p.e[0],'k-')
            ax2.plot(p.s[1],p.s[0],'k-')
        return fig

def load_a_fc(fn=None,hash_code=None):
    """
    Arguments
    ---------
    fn = None
    hash_code = None
    """
    p=os.popen('ls FLDA_MACRO_*_*_%s'%hash_code)
    elems = p.readlines()
    if len(elems)!=1:
        raise IOError, 'Could not find the file'
    fn = elems[0].split('\n')[0]
    with open(fn,'rb') as fi:
        data_collection_FLDA = pickle.load(fi)
    print '%s was restored'%fn

    P = PathCollection()
    for i in xrange(len(data_collection_FLDA)):
        e6,s6,ebar,sbar,d6,t,dt,debar = data_collection_FLDA[i]
        p = Path(s6,e6,d6,t,dt,ebar,sbar,debar)
        P.add_path(p)

    ## draw
    fig = P.draw();fn='regiona.pdf'
    fig.savefig(fn)
    print '%s has been saved'%fn
    return P

def func_k2(psi, rho, gamma):
    return sin(psi)**2 + rho * cos(psi)**2 - 2*gamma*sin(psi)*cos(psi)
def func_k3(psi,alpha,beta):
    return cos(psi)**2 + alpha*sin(psi)**2 + 2*beta*sin(psi)*cos(psi)
def eq17(alpha_a,beta_a,alpha_b,psi):
    aa = alpha_a
    ba = beta_a
    ab = alpha_b
    c = cos(psi)
    s = sin(psi)
    ## return beta^b
    return (ba* (ab*s**2 - c**2) + s*c*(ab-aa)) / (aa*s**2-c**2)

def eq20(da, eta_a,k1_a,rho_a,rho_b,k2_a,k2_b,D):
    dD = da / (eta_a*k1_a) * (1+rho_a-(1+rho_b)*k2_a/k2_b)*D
    return dD

def eq18(eta_a,eta_b,k3_a,k3_b,D,wa,wb):
    if eta_b==0 or k3_b==0 or wb==0:
        # print eta_b, k3_b, wb
        raise IOError
    return eta_a  * k3_a / (eta_b*k3_b) - D * wb / wa

def update_psi(psi,A,da_):
    """
    psi
    A
    da_
    """
    f = np.tan(psi)*(A.eta*A.k1+da_) / (A.eta*A.k1+A.rho*da_)
    new = np.arctan(f)
    dpsi = new - psi
    return dpsi

def Barlat_objf(P,psi,sa,da,f,af,yfunc,wa,wb):
    def objf(alpha_b,guess_stt):
        """
        Guess alpha_b
        """
        #
        sa_grv = lib.rot_6d(sa.copy(),psi)
        sa_nn,sa_nt = sa_grv[0],sa_grv[5]
        sb_nn = sa_nn/f
        sb_nt = sa_nt/f

        if type(guess_stt).__name__=='NoneType':
            x = sa_grv[1] ## guess as sa_tt
        else:
            x=guess_stt

        nit = 1;mx=50
        dx  = 1.e-7
        tol = 1.e-9

        #verbose=True
        verbose=False

        # fn=os.path.join(lib.find_tmp(),'MKobjf_%s.log'lib.gen_hash_code(6))
        # fo=open(fn,'w')

        hist_iter = ''
        while (True):
            sb_grv0 = np.array([sb_nn,x,0,0,0,sb_nt])
            sb0 = lib.rot_6d(sb_grv0,-psi)
            alpha0 = sb0[1]/sb0[0]
            beta0  = sb0[5]/sb0[0]
            beta_b0 = eq17(P.alpha,P.beta,alpha0,psi)
            F0      = alpha0 - alpha_b
            F0_beta = beta0 - beta_b0

            sb_grv1 = np.array([sb_nn,x+dx,0,0,0,sb_nt])
            sb1     = lib.rot_6d(sb_grv1,-psi)
            alpha1  = sb1[1]/sb1[0]
            beta1   = sb1[5]/sb1[0]
            beta_b1 = eq17(P.alpha,P.beta,alpha1,psi)
            F1      = alpha1 - alpha_b
            F1_beta = beta1 - beta_b1

            head = ('%10s'*11)%('alpha0','beta0','alpha_b','alpha1',
                                 'beta1','beta_b0','beta_b1',
                                 'F0','F1','F0_beta','F1_beta')
            fmt = '%10.2e'*11
            cnt = fmt%(alpha0, beta0, alpha_b,beta_b0, alpha1,beta1,beta_b1,\
                       F0,F1,F0_beta,F1_beta)
            if nit==1:
                hist_iter=hist_iter+head
            else:
                hist_iter=hist_iter+'\n'+cnt

            # if nit==1: print head
            # print cnt


            if (abs(F0)<tol):
                if verbose:
                    if nit==1: print '%s'%'nit',  '%5s'%'x', \
                       '%8s'%'F0','%8s'%'alph_b','%8s'%'alph_0',\
                       ('%7s'*3)%('s1','s2','s6')
                    print '%2.2i'%nit, '%5.1f'%x, '%8.1e'%F0,\
                        '%8.1e'%alpha_b, '%8.1e'%alpha0, \
                        ('%7.1f'*3)%(sb0[0],sb0[1],sb0[5])
                break
            if verbose and nit==1: print '***'*30
            if (nit>mx):
                print hist_iter
                raise IOError,'Cannot find solution'

            jac_new = (F1-F0)/dx
            if jac_new==0:
                if nit==1: raise IOError
                else: jac_new = jac ## use the previous jacobian
            else: jac = jac_new
            x = x - F0/jac

            if verbose:
                if nit==1: print head
                print cnt

            nit=nit+1

        #print ('%5s'*3)%('s1','s2','s6')
        #print ('%9.5f'*3)%(sb0[0],sb0[1],sb0[5])

        # if verbose: print '***'*30
        sb = sb0[::]
        db = af(sb) ## is not multiplied by lambda yet
        alpha_b = sb[1]/sb[0]
        beta_b  = sb[5]/sb[0]
        rho_b   = db[1]/db[0]
        gamma_b = db[5]/db[0]
        eta_b   = sb[0] / yfunc(sb)
        k2_b = func_k2(psi,rho_b,gamma_b)
        k3_b = func_k3(psi,alpha_b,beta_b)

        # da_grv = lib.rot_6d(da,psi)
        # db_grv = lib.rot_6d(db,psi)

        # # print da_grv
        # # print db_grv

        # dum = da_grv[1] / db_grv[1]
        # db_grv_new = db_grv * dum
        # db = lib.rot_6d(db_grv_new,-psi)

        class B: pass
        B.s = sb
        B.d = db
        B.alpha = alpha_b
        B.beta = beta_b
        B.rho = rho_b
        B.gamma = gamma_b
        B.k1 = 1 + alpha_b * rho_b + 2*beta_b * gamma_b
        B.k2 = k2_b
        B.k3 = k3_b
        B.eta = eta_b

        return eq18(P.eta,B.eta,P.k3,B.k3,f,wa,wb),B
    return objf

def main(FLDA_pck_name=None,
         psi0=0,f0=0.990):
    """
    Arguments
    ---------
    psi0 =0.
    f0   =0.990
    """
    t0 = time.time()
    ## Choice of material / loading condition
    mat = materials.iso_metal_vm()
    # mat = materials.iso_metal_hf8()
    uet(time.time() - t0,'Time for MAT');print
    """
    mat.func_hd as a function of ebar
    mat.func_sr as a function of ebar dot
    mat.func_yd as a function of cauchy stress
    """
    t0 = time.time()
    bnd = materials.prop_loading_long()
    # bnd = materials.bb() ## balanced biaxial
    uet(time.time() - t0,'Time for BND');print
    # bnd = materials.prop_loading_refine()

    if type(FLDA_pck_name).__name__=='NoneType':
        ## Save FLDa
        t0 = time.time()
        FLDA_pck_name = save_a(mat, bnd)
        uet(time.time() - t0,'Time for FLA');print

    hash_a = FLDA_pck_name[::-1][:6][::-1]
    ## Load from pickles.
    Paths = load_a_fc(hash_code=hash_a)

    print 'npath:', Paths.npath, '\n'

    ## Boundary condition
    beta    = bnd.beta
    sr_eq   = bnd.sr_eq
    dbar    = bnd.dbar
    ebar_mx = bnd.ebar_mx

    head= '%4s'%'step'+('%9s'+'%12s'+'%7s'*4)%(
        'Da','Db','Ea','Eb','Wa','Wb')+\
          ('%10s'*4)%('D1','D2','D3','D6')+\
          ('%10s'*4)%('E1','E2','E3','E6')+\
          ('%7s'*3)%('S1','S2','S6')+\
          '|'+\
          ('%10s'*4)%('D1','D2','D3','D6')+\
          ('%10s'*4)%('E1','E2','E3','E6')+\
          ('%7s'*3)%('S1','S2','S6')+\
          '%7s'%'f'+'%5s'%'psi'+'%8s'%'ratio'

    head=head+'\n'+'-'*276

    for ipath in xrange(Paths.npath):
        P   = Paths.paths[ipath]
        psi = psi0 * np.pi/180.
        f   = f0 * 1.0

        ## accumulative 6D strain for region B
        eb  = np.zeros(6)

        ## Equivalent quantifies for stress, strain and strain rate
        sb_ = 0.; eb_ = 0.; db_ = 0.

        ## history
        class matA: pass
        class matB: pass
        matA.s=[]; matA.e=[]; matA.d=[]
        matB.s=[]; matB.e=[]; matB.d=[]

        ## Region A, Region B (current step)
        #class A: pass; class B: pass
        A=Region(mat); B=Region(mat); A.ebar=0; B.ebar=0

        ## log file
        fn_log = os.path.join(lib.find_tmp(),'%s.log'%lib.gen_hash_code(6))
        fo_log = open(fn_log,'w')
        fo_log.write(head+'\n')

        t0 = time.time()
        for istp in xrange(P.nstp):
            t  = P.t[istp]
            dt = P.dt[istp]

            ## assign properties to A
            A.assign(P.s[:,istp],P.d[:,istp],P.e[:,istp])
            A.calc_s_grv(psi=psi) # A.s_grv

            ## Guess properties of B
            B.s = A.s[::]; B.d = A.d[::]

            ## flow stress from hardening

            ## guess 4 initial values for region B
            B.x =[1.,1.,0.,0.] ## X_b, Y_b, Z_b, E0_b
            newton_raph_fld(A,B,psi,f)
            B.x # X_b,Y_b,Z_b,E_b
            print B.x
            raise IOError

def newton_raph_fld(A,B,psi,f):
    """
    Arguments
    ---------
    A
    B
    psi
    f
    """
    from mk_for import gauss
    residu = 1.0; it = 0; itmax=100
    eps = 1e-10 ## convergence criterion
    while (residu>eps and it<itmax):
        F, J   = fonc_fld(A,B,psi,f) ##
        ndim   = len(F)
        res    = gauss(ndim=ndim,a=J,b=F)

        print it, res, B.x

        if np.isnan(res).any():
            raise IOError
        B.x    = B.x - res
        residu = np.dot(res,res)
        it = it + 1

def fonc_fld(A,B,psi,f):
    """
    Function F and matrix J (jacobian)

    Arguments
    ---------
    A
    B
    psi
    f

    Return
    ------
    F, J
    """
    C=np.cos(psi); S=np.sin(psi); SC=C*S; C2=C**2; S2=S**2

    ## guess stress of B in groove
    # B.x[0] sxx
    # B.x[1] syy
    # B.x[2] sxy
    # B.x[3] ebar (equivalent strain)

    ## guess (tilde) stress
    B.sx     = np.array([B.x[0],B.x[1],0., 0., 0., B.x[2]]) ## grv
    B.s      = calc_grv(B.sx,-psi) # B stress in the lab axes

    B.phi    = B.func_yd( B.s)
    B.phi_d1 = B.func_yd1(B.s) ## According to af, its strain rate direction
    B.phi_d2 = B.func_yd2(B.s)

    ## Strain rates referred in the groove axes
    B.d_grv = lib.rot_6d(B.phi_d1,psi)

    ## 2nd derivative in the groove
    B.F2xB      = np.zeros((6,6)) ## needs rotate to the groove
    # ------------------------------------------------------------------ #
    B.F2xB[0,0]=    B.phi_d2[0,0]*C2+B.phi_d2[0,1]*S2 +B.phi_d2[0,5]*SC
    B.F2xB[1,0]=    B.phi_d2[1,0]*C2+B.phi_d2[1,1]*S2 +B.phi_d2[1,5]*SC
    B.F2xB[5,0]=   (B.phi_d2[5,0]*C2+B.phi_d2[5,1]*S2 +B.phi_d2[5,5]*SC)/2.
    B.F2xB[0,1]=    B.phi_d2[0,0]*S2+B.phi_d2[0,1]*C2 -B.phi_d2[0,5]*SC
    B.F2xB[1,1]=    B.phi_d2[1,0]*S2+B.phi_d2[1,1]*C2 -B.phi_d2[1,5]*SC
    B.F2xB[5,1]=   (B.phi_d2[5,0]*S2+B.phi_d2[5,1]*C2 -B.phi_d2[5,5]*SC)/2.
    B.F2xB[0,5]= 2*(B.phi_d2[0,1]   -B.phi_d2[0,0])*SC+B.phi_d2[0,5]*(C2-S2)
    B.F2xB[1,5]= 2*(B.phi_d2[1,1]   -B.phi_d2[1,0])*SC+B.phi_d2[1,5]*(C2-S2)
    B.F2xB[5,5]=(2*(B.phi_d2[5,1]   -B.phi_d2[5,0])*SC+B.phi_d2[5,5]*(C2-S2))/2.

    # hardening parameters
    sigb  = B.func_hd(  B.x[3])
    dsigb = B.func_hd_d(B.x[3])
    siga  = A.func_hd(A.ebar)
    if np.isnan(sigb): raise IOError

    sr_ref = 1000.;  sr_m_a = 0.05;  sr_m_b = 0.05; dm = 0.

    ## objective functions (4)
    A.beta_grv = A.s_grv[5]/A.s_grv[0]
    ## ------------
    F    = np.zeros(4)
    F[0] = B.d_grv[1] ## d_tt
    F[1] = B.x[2] - B.x[1]*A.beta_grv
    F[2] = B.phi  - 1.0
    F[3] = -np.log(f*sigb/siga)+\
           (sr_m_a-sr_m_b)*np.log(sr_ref)-B.x[3]*B.d_grv[0]
    ## ------------

    if np.isnan(F).any(): raise IOError

    ## Jacobian
    J      = np.zeros((4,4))
    J[0,0] = S2  * B.F2xB[0,0] + C2 * B.F2xB[1,0] - 2*SC*B.F2xB[5,0]
    J[0,1] = S2  * B.F2xB[0,1] + C2 * B.F2xB[1,1] - 2*SC*B.F2xB[5,1]
    J[0,2] = S2  * B.F2xB[0,5] + C2 * B.F2xB[1,5] - 2*SC*B.F2xB[5,5]
    J[1,0] = -A.beta_grv
    J[1,2] = 1.
    J[2,0] = B.d_grv[0] ## ENNB
    J[2,1] = B.d_grv[1] ## ETTB
    J[2,2] = B.d_grv[5] *2 ## 2 ENTB
    J[3,3] = -dsigb / sigb - dm*np.log(sr_ref) - B.d_grv[0]

    return F, J




def main2(FLDA_pck_name='FLDA_MACRO_20160414_135512_9b435d',
         psi0=0,f0=0.990):
    """
    Reference:
    Crystallographic Texture, anisotropic yield
    surfaces and forming limits of sheet metals,
    F. Barlat, MSE, Vol 91, 1987

    Arguments
    ---------
    psi0 =0.
    f0   =0.990
    """
    t0 = time.time()
    ## Choice of material / loading condition
    mat = materials.iso_metal_vm()
    # mat = materials.iso_metal_hf8()
    uet(time.time() - t0,'Time for MAT');print

    """
    mat.func_hd as a function of ebar
    mat.func_sr as a function of ebar dot
    mat.func_yd as a function of cauchy stress
    """
    t0 = time.time()
    bnd = materials.prop_loading_long()
    # bnd = materials.bb() ## balanced biaxial
    uet(time.time() - t0,'Time for BND');print
    # bnd = materials.prop_loading_refine()

    if type(FLDA_pck_name).__name__=='NoneType':
        ## Save FLDa
        t0 = time.time()
        FLDA_pck_name = save_a(mat, bnd)
        uet(time.time() - t0,'Time for FLA');print

    hash_a = FLDA_pck_name[::-1][:6][::-1]
    ## Load from pickles.
    Paths = load_a_fc(hash_code=hash_a)

    print 'npath:', Paths.npath, '\n'

    ## Boundary condition
    beta    = bnd.beta
    sr_eq   = bnd.sr_eq
    dbar    = bnd.dbar
    ebar_mx = bnd.ebar_mx

    head= '%4s'%'step'+('%9s'+'%12s'+'%7s'*4)%(
        'Da','Db','Ea','Eb','Wa','Wb')+\
          ('%10s'*4)%('D1','D2','D3','D6')+\
          ('%10s'*4)%('E1','E2','E3','E6')+\
          ('%7s'*3)%('S1','S2','S6')+\
          '|'+\
          ('%10s'*4)%('D1','D2','D3','D6')+\
          ('%10s'*4)%('E1','E2','E3','E6')+\
          ('%7s'*3)%('S1','S2','S6')+\
          '%7s'%'f'+'%5s'%'psi'+'%8s'%'ratio'

    head=head+'\n'+'-'*276

    for ipath in xrange(Paths.npath):
        P   = Paths.paths[ipath]
        psi = psi0 * np.pi/180.
        f   = f0 * 1.0

        ## accumulative 6D strain for region B
        eb  = np.zeros(6)

        ## Equivalent quantifies for stress, strain and strain rate
        sb_ = 0.; eb_ = 0.; db_ = 0.

        ## Region B is saved as a class for path case.
        class matB: pass
        class A: pass
        matB.s=[];matB.e=[];matB.d=[]

        ## log file
        fn_log = os.path.join(lib.find_tmp(),'%s.log'%lib.gen_hash_code(6))
        fo_log = open(fn_log,'w')
        fo_log.write(head+'\n')

        t0 = time.time()
        for istp in xrange(P.nstp):
            t  = P.t[istp]
            dt = P.dt[istp]
            sa = P.s[:,istp]
            ea = P.e[:,istp]
            da = P.d[:,istp]

            wrate_a = np.dot(sa,da)
            da_ = wrate_a/P.sbar[istp]
            sa_ = P.sbar[istp]
            ea_ = P.ebar[istp]
            P.values_mk(istp,psi)

            ## Guess alpha for region B, -> obtaining also beta for region B
            ## This completes tilde stress state, which gives rho, gamma, and eta for region B.
            ## solve equation 18.
            wa = mat.func_hd(ea_);wb = mat.func_hd(eb_)

            ## function of alpha_b
            objf = Barlat_objf(P,psi,sa,da,f,mat.af,mat.func_yd,wa,wb)

            if istp==0:
                #x0 = 0. ## guessed alpha_b
                x0 = P.alpha
                y0 = None
            else:
                x0  = B.s[1]/B.s[0] ## guess alpha from previous result
                s33 = lib.s62c(B.s[::])
                nv,tv = lib.nt_psi(psi)
                stt = np.dot(np.dot(s33,tv),tv)
                y0  = stt

            ## Numerical conditions
            dx  = 5e-7
            tol = 1e-9
            nit = 1
            nmx = 50

            verbose=False
            # verbose=True
            head_iter = '   '+('%3s'+'%9s'*6)%(
                'stp','x0','Stt','F','s11','s22','s12')
            hist_iter=''

            while (True):
                if (nit>nmx):
                    if nit>1:
                        print head_iter
                        print hist_iter
                        #print cnt_iter

                    #raise IOError
                    break

                try:
                    rst = objf(x0,y0)
                except:
                    fo_log.close()
                    print fn_log
                    import plotter
                    _fn_ = plotter.plot_log2(fn=fn_log,yfunc=mat.func_yd)
                    os.system('open %s'%_fn_)
                    raise IOError


                F,B = rst
                if verbose:
                    if nit==1: print head_iter
                    print cnt_iter
                if (abs(F)<tol): break
                F0  = objf(x0-dx,y0)[0]
                F1  = objf(x0+dx,y0)[0]
                jac = (F1-F0)/(dx*2)
                if (jac==0): raise IOError
                x0 = x0 - F/jac

                s33=lib.s62c(B.s[::])
                nv,tv=lib.nt_psi(psi)
                stt = np.dot(np.dot(s33,tv),tv)
                y0 = stt
                cnt_iter= '-*-'+('%3i'+'%9.3f'*2+'%9.2e'+'%9.3f'*3)%(
                    nit,x0,y0,F,B.s[0],B.s[1],B.s[5])
                hist_iter=hist_iter+'\n'+cnt_iter
                nit = nit + 1

            # raise IOError
            B       = rst[1]
            sb, db  = B.s, B.d

            ## Normalize B.d such that B.dtt = A.dtt
            A.d=P.d[:,istp]
            B.d33 = lib.s62c(B.d[::])
            B.dtt = np.dot(np.dot(B.d33,tv),tv)
            A.d33 = lib.s62c(A.d[::])
            A.dtt = np.dot(np.dot(A.d33,tv),tv)
            scale_factor = B.dtt/A.dtt
            B.d = B.d/scale_factor
            wrate_b = np.dot(B.s,B.d)
            db_    = wrate_b / mat.func_yd(B.s)


            # ## eq 15
            # db_ = (B.eta*B.k1*P.k2)/(P.eta*P.k1*B.k2)*da_

            # # if db_<da_:
            # #     raise IOError

            # B.d=B.d*db_
            # A.d=P.d[:,istp]

            # ## check compatibility
            # B.d33 = lib.s62c(B.d[::])
            # B.dtt = np.dot(np.dot(B.d33,tv),tv)
            # A.d33 = lib.s62c(A.d[::])
            # A.dtt = np.dot(np.dot(A.d33,tv),tv)

            # print 'B.dtt', B.dtt
            # print 'A.dtt', A.dtt
            # print 'B.dtt-A.dtt',B.dtt-A.dtt
            # f = B.dtt/A.dtt
            # B.d33 = A.d33 / f

            # ## compatability



            # if db_1<da_:
            #     raise IOError

            ratio   = B.d[2]/P.d[2,istp]

            ## inhomogeneity
            df = eq20(da_,P.eta,P.k1,P.rho,B.rho,P.k2,B.k2,f)
            f = f + df * dt
            ## rotation
            psi = psi + update_psi(psi,P,da_) * dt

            if np.mod(istp,20)==0: print head

            stamp= '%4.4i'%istp+ '%9.1e'%da_+'%12.4e'%db_+\
                   ('%7.3f'*2+'%7.1f'*2)%(ea_,eb_,wa,wb)+\
                   ('%10.2e'*4)%(P.d[0,istp],P.d[1,istp],
                                 P.d[2,istp],P.d[5,istp])+\
                   ('%10.2e'*4)%(P.e[0,istp],P.e[1,istp],
                                 P.e[2,istp],P.e[5,istp])+\
                   ('%7.1f'*3 )%(P.s[0,istp],P.s[1,istp],
                                 P.s[5,istp])+' '+\
                   ('%10.2e'*4)%(B.d[0],    B.d[1],
                                 B.d[2],     B.d[5])+\
                   ('%10.2e'*4)%( eb[0],     eb[1],
                                  eb[2],      eb[5])+\
                   ('%7.1f'*3 )%(B.s[0],    B.s[1],
                                 B.s[5])+\
                   '%7.3f'%f+'%5.1f'%(psi*180/np.pi)+'%8.3f'%ratio
            print stamp
            fo_log.write(stamp+'\n')

            eb      = eb + B.d * dt
            B.e     = eb
            eb_     = eb_ + db_ * dt
            B.eb_   = eb_

            matB.s.append(B.s)
            matB.e.append(eb)
            matB.d.append(B.d)

            if ratio>10:
                fo_log.close()
                print fn_log
                import plotter
                _fn_ = plotter.plot_log2(fn=fn_log,yfunc=mat.func_yd)
                os.system('open %s'%_fn_)

                raise IOError
                #break

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fn',   type=str,help='file name for FLDA')
    parser.add_argument(
        '--psi',   type=float,help='Initial angle of the groove')
    parser.add_argument(
        '--f',   type=float,help='Initial inhomogeneity factor')
    args        = parser.parse_args()
    main(FLDA_pck_name=args.fn,psi0=args.psi,f0=args.f)
