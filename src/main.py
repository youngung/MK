from fldb import *
import modules; reload(modules)
from scipy import optimize
import os, glob
from MP.mat import mech
from MP import progress_bar
import lib
uet=progress_bar.update_elapsed_time

class Path:
    """
    Collection of loading history

    s: 6D sigma (stress)
    e: 6D strain accumulative strain
    d: 6D strian rate
    t     time stamps
    """
    def __init__(self,s,e,d,t,dt,ebar,sbar):
        self.s=s.T
        self.e=e.T
        self.d=d.T
        self.t=t
        self.dt=dt
        self.ebar=ebar
        self.sbar=sbar
        self.nstp = len(t)
        pass
    def rot_istp(self,istp,psi):
        s = lib.rot_6d(self.s[:,istp],psi)
        e = lib.rot_6d(self.e[:,istp],psi)
        d = lib.rot_6d(self.d[:,istp],psi)
        return s,e,d

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
        e6,s6,ebar,sbar,d6,t,dt = data_collection_FLDA[i]
        p = Path(s6,e6,d6,t,dt,ebar,sbar)
        P.add_path(p)

    ## draw
    fig = P.draw();fn='regiona.pdf'
    fig.savefig(fn)
    print '%s has been saved'%fn

    return P

def main(FLDA_pck_name='FLDA_MACRO_20160414_135512_9b435d'):
    t0 = time.time()
    ## Choice of material / loading condition
    mat = materials.iso_metal_vm()
    #mat = materials.iso_metal_hf8()
    uet(time.time() - t0,'Time for MAT');print

    """
    mat.func_hd as a function of ebar
    mat.func_sr as a function of ebar dot
    mat.func_yd as a function of cauchy stress
    """
    t0 = time.time()
    bnd = materials.prop_loading_long()
    uet(time.time() - t0,'Time for BND');print
    # bnd = materials.prop_loading_refine()
    # bnd = materials.bb() ## balanced biaxial

    if type(FLDA_pck_name).__name__=='NoneType':
        ## Save FLDa
        t0 = time.time()
        FLDA_pck_name = save_a(mat, bnd)
        uet(time.time() - t0,'Time for FLA');print

    hash_a = FLDA_pck_name[::-1][:6][::-1]
    ## Load from pickles.
    Paths = load_a_fc(hash_code=hash_a)

    print 'npath:', Paths.npath

    ## Boundary condition
    beta    = bnd.beta
    sr_eq   = bnd.sr_eq
    dbar    = bnd.dbar
    ebar_mx = bnd.ebar_mx

    fmt_head='%5i'*1+'%8.4f'*2+'%5.1f'*1+'%8.4f'*2+'%8.1f'*2+\
        '%8.4f'*2+'%8.2f'*2+'%3s'*1+'%8.3f'*4+'%10.2e'*4+\
        '%8.1f'*4+'%4s'+'%8.3f'*4+'%10.2e'*4+'%8.1f'*4
    head= ('%5s'+'%8s'*2+'%5s'+'%8s'*8+'%3s'*1+'%8s'*4+\
           '%10s'*4+'%8s'*4+'%4s'+'%8s'*4+'%10s'*4+'%8s'*4)%(
               'stp','x0','f','df','Eeq_A','Eeq_B','Seq_A','Seq_B',
               'EdeqA','EdeqB','ratio','psi','|',
               'E11_a','E22_a','E12_a','E33_a',
               'D11_a','D22_a','D12_a','D33_a',
               'S11_a','S22_a','S12_a','S33_a','||',
               'E11_b','E22_b','E12_b','E33_b',
               'D11_b','D22_b','D12_b','D33_b',
               'S11_b','S22_b','S12_b','S33_b')

    for ipath in xrange(Paths.npath):
        P = Paths.paths[ipath]
        psi = 0.
        f0 = 0.990
        f = f0

        da_ = 0.
        eb=np.zeros(6)
        sb_=0.
        eb_=0.
        db_=0.
        for istp in xrange(P.nstp):
            t  = P.t[istp]
            dt = P.dt[istp]
            sa = P.s[:,istp]
            ea = P.e[:,istp]
            da = P.d[:,istp]
            sa_ = P.sbar[istp]
            ea_ = P.ebar[istp]

            sag,eag,dag = P.rot_istp(istp,psi)## groove

            # FIND_S22
            # objf = modules.find_s22(
            #     sag[0],sag[1],sag[5],dag[0],dag[1],dag[5],
            #     f=f,af=mat.af,func_yld=mat.func_yd)
            # FIND_S2B

            objf = modules.find_s2b(
                H=mat.func_hd(eb_), g=mat.func_sr, sb=sag/f,
                db=dag[::],func_yld=mat.func_yd,af=mat.af)

            x0 = sag[1]
            dx = 1e-6
            nit = 0; it_mx = 20; tol=1e-5
            # print ('%2s '+'%9s'*7)%('nt','f_0','s1','s2','s6','d1','d2','d6')
            if db_==0: y0 = 1e-3
            else: y0 = db_
            while (True):
                if (np.isnan(x0)): raise IOError

                f_0,deq,sb,db,dum = objf(x0,y0)#[0]
                # print ('%2.2i '+'%9.1f'*4+'%9.1e'*3)%(nit,f_0,sb[0],sb[1],sb[5],db[0],db[1],db[5])
                if (abs(f_0)<tol): break
                f_1 = objf(x0-dx,y0)[0]
                f_2 = objf(x0+dx,y0)[0]

                jac = (f_2-f_1)/(dx*2)
                x0 = x0 - f_0/jac
                nit = nit + 1


            sb_ = mat.func_yd(sb)
            sa_ = mat.func_yd(sa)

            raise IOError

            wrate_b = 0.;wrate_a = 0.
            for i in xrange(3):
                wrate_b = wrate_b + sb[i]*db[i]+ sb[i+3]*db[i+3]*2
                wrate_a = wrate_a + sa[i]*da[i]+ sa[i+3]*da[i+3]*2

            db_ = wrate_b/sb_
            da_ = wrate_a/sa_

            # raise IOError
            eb_ = eb_ + db_ * dt

            # eb_=0.
            # sb_=0.


            fmt_s= '%5s '+'%5s '+\
                   '%8s '*3+\
                   '%8s '*11+\
                   '%4s '+\
                   '%8s '*3+\
                   '%8s '*11
            fmt_v= '%5.5i '+'%5.3f '+\
                   '%8.1f '+'%8.1e '*2+\
                   '%8.1f '*3+'%8.1e '*8+\
                   '%4s '+\
                   '%8.1f '+'%8.1e '*2+\
                   '%8.1f '*3+'%8.1e '*8
            if np.mod(istp,40)==0: print fmt_s%(
                    'step','f',
                    's_','e_','d_',
                    's1','s2','s6',
                    'd1','d2','d3','d6',
                    'e1','e2','e3','e6','||',
                    's_','e_','d_',
                    's1','s2','s6',
                    'd1','d2','d3','d6',
                    'e1','e2','e3','e6')
            print fmt_v%(istp,f,
                         sa_,ea_,da_,
                         sa[0],sa[1],sa[5],
                         da[0],da[1],da[2],da[5],
                         ea[0],ea[1],ea[2],ea[5],'||',
                         sb_,eb_,db_,
                         sb[0],sb[1],sb[5],
                         db[0],db[1],db[2],db[5],
                         eb[0],eb[1],eb[2],eb[5])

            ## updates
            eb=eb+db*dt

            f = f0 * np.exp(eb[2]-ea[2])
            # raise IOError
            if istp>20: raise IOError

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fn',   type=str,help='file name')
    args        = parser.parse_args()
    main(FLDA_pck_name=args.fn)
