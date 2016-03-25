"""
Materials database
"""
from flda import *
## Material characteristics
class Mat_A:
    def __init__(self,name):
        self.name=name
        pass
    def assign_hd(self,ihd_opt,**kwargs):
        self.func_hd = c_G(ihd_opt,**kwargs)
        pass
    def assign_sr(self,isr_opt,**kwargs):
        self.func_sr = c_F(isr_opt,**kwargs)
        pass
    def assign_yd(self,iyd_opt,**kwargs):
        self.func_yd = return_c_yld(iyd_opt,**kwargs)

def iso_metal_vm():
    """
    Isotropic Von Mises metal
    """
    my_iso_metal = Mat_A(name='iso_metal')
    my_iso_metal.assign_hd(ihd_opt=0,k=6.00e2,eps_0=4.23e-4,n=2.55e-1)
    my_iso_metal.assign_sr(isr_opt=0,ed0=1e-3,m=0.02)
    my_iso_metal.assign_yd(iyd_opt=1) ## Von Mises yield surface
    return my_iso_metal


def tension_test():
    eps_eq = np.linspace(0,0.2,100)
    eps_sr = [1e-3,1e-2,1e-1]

    my_mat = iso_metal_vm()
    F=my_mat.func_sr
    G=my_mat.func_hd

    fig=plt.figure()
    ax=fig.add_subplot(111)
    for j in xrange(len(eps_sr)):
        S=[];E=[]
        for i in xrange(len(eps_eq)):
            eps = eps_eq[i]
            edt = eps_sr[j]
            sig = F(edt) * G(eps)
            E.append(eps)
            S.append(sig)
        ax.plot(E,S,label=edt)

    ax.set_ylim(0.,)
    ax.legend(loc='best')
    fig.savefig('tension_test_mat.pdf',bbox_inches='tight')

class Bnd_A:
    def __init__(self,name,alpha0,alpha1,alphan,
                 beta,sr_eq,dbar,ebar_mx):
        self.name    = name
        self.alphas  = np.linspace(alpha0,alpha1,alphan)
        self.nprob   = alphan
        self.beta    = beta
        self.sr_eq   = sr_eq
        self.dbar    = dbar
        self.ebar_mx = ebar_mx

def prop_loading_short():
    """
    short proportional loading
    """
    my_loading = Bnd_A(name='proportional_Loading1',
                       alpha0=0,alpha1=1,alphan=3,beta=0,
                       sr_eq=1e-3,dbar=1e-6,ebar_mx=0.1)
    return my_loading

def prop_loading_long():
    """
    long proportional loading
    """
    my_loading = Bnd_A(name='proportional_Loading1',
                       alpha0=0,alpha1=1,alphan=3,beta=0,
                       sr_eq=1e-3,dbar=1e-4,ebar_mx=1.5)
    return my_loading
