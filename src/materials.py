"""
Collection of materials
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
        self.func_yld = return_c_yld(iyd_opt,**kwargs)


class Bnd_A:
    def __init__(self,name,alpha0,alpha1,alphan,
                 beta,sr_eq,dbar,ebar_mx):
        self.name=name
        self.alphas = np.linspace(alpha0,alpha1,alphan)
        self.nprob = alphan
        self.n     = self.nprob
        self.beta=beta
        self.sr_eq=sr_eq
        self.dbar=dbar
        self.ebar_mx=ebar_mx

def iso_metal_vm():
    """
    Isotropic Von Mises metal
    """
    my_iso_metal = Mat_A(name='iso_metal')
    my_iso_metal.assign_hd(ihd_opt=0,k=6.00e2,eps_0=4.23e-4,n=2.55e-1)
    my_iso_metal.assign_sr(isr_opt=0,ed0=1e-3,m=1.2)
    my_iso_metal.assign_yd(iyd_opt=1) ## Von Mises yield surface
    return my_iso_metal

def prop_loading_short():
    """
    short proportional loading
    """
    my_loading = Bnd_A(name='proportional_Loading1',
                       alpha0=0,alpha1=1,alphan=3,beta=0,
                       sr_eq=1e-3,dbar=1e-4,ebar_mx=0.1)
    return my_loading

def prop_loading_long():
    """
    long proportional loading
    """
    my_loading = Bnd_A(name='proportional_Loading1',
                       alpha0=0,alpha1=1,alphan=3,beta=0,
                       sr_eq=1e-3,dbar=1e-4,ebar_mx=1.5)
    return my_loading
