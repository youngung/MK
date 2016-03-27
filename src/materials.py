"""
Materials database
"""

import lib, flda
reload(lib)
reload(flda)
from flda import *

## Material Object
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
    def calc_ys(self):
        """
        yield lous in pi-plane
        """
        return lib.y_locus_c(100,self.func_yd)
    def tension_test(self):
        """
        """
        eps_eq = np.linspace(0,0.2,100)
        eps_sr = [1e-3,1e-2,1e-1]
        F=self.func_sr; G=self.func_hd
        fig=plt.figure(figsize=(7,3));
        ax1=fig.add_subplot(121)
        ax2=fig.add_subplot(122)
        for j in xrange(len(eps_sr)):
            S=[];E=[]
            for i in xrange(len(eps_eq)):
                eps = eps_eq[i]
                edt = eps_sr[j]
                sig = F(edt) * G(eps)
                E.append(eps)
                S.append(sig)
                pass
            ax1.plot(E,S,label=r'$\mathrm{\dot{\bar{E}}^{eq}}=10^{%i}$'%np.log10(edt))
            pass

        ys_ps, ys_pi = self.calc_ys()
        for j in xrange(len(eps_sr)):
            edt = eps_sr[j]
            ax2.plot(ys_ps[0]*F(edt),ys_ps[1]*F(edt),
                     label=r'$\mathrm{\dot{\bar{E}}^{eq}}=10^{%i}$'%np.log10(edt))

        ax2.set_xlabel(r'$\mathrm{\bar{\Sigma}^{11}}$',dict(fontsize=20))
        ax2.set_ylabel(r'$\mathrm{\bar{\Sigma}^{22}}$',dict(fontsize=20))

        ax1.set_ylim(0.,)
        ax1.legend(loc='best')
        ax1.set_title(self.name)
        ax1.set_xlabel(r'$\mathrm{\bar{E}^{eq}}$',dict(fontsize=20))
        ax1.set_ylabel(r'$\mathrm{\bar{\Sigma}^{eq}}$',dict(fontsize=20))
        plt.tight_layout()
        fig.savefig('tension_test_%s.pdf'%self.name,bbox_inches='tight')
        plt.close(fig)

## Material stocks
def iso_metal_vm():
    """
    Isotropic Von Mises metal
    """
    my_iso_metal = Mat_A(name='iso_vm')
    my_iso_metal.assign_hd(ihd_opt=0,k=6.00e2,eps_0=4.23e-4,n=2.55e-1)
    my_iso_metal.assign_sr(isr_opt=0,ed0=1e-3,m=0.10)
    my_iso_metal.assign_yd(iyd_opt=1) ## Von Mises yield surface
    return my_iso_metal

def iso_metal_hf8():
    """
    Isotropic Hosford
    """
    my_iso_metal = Mat_A(name='iso_hf')
    my_iso_metal.assign_hd(ihd_opt=0,k=6.00e2,eps_0=4.23e-4,n=2.55e-1)
    my_iso_metal.assign_sr(isr_opt=0,ed0=1e-3,m=0.10)
    my_iso_metal.assign_yd(iyd_opt=2,a=8) ## Hosford exponent of 8
    return my_iso_metal

def aniso_metal_hill():
    """
    Anisotropic Hill
    """
    my_aniso_metal = Mat_A(name='aniso_hill')
    my_aniso_metal.assign_hd(ihd_opt=0,k=6.00e2,eps_0=4.23e-4,n=2.55e-1)
    my_aniso_metal.assign_sr(isr_opt=0,ed0=1e-3,m=0.10)
    ## Anisotropic Hill's quadratic tuned for the IF steel
    my_aniso_metal.assign_yd(iyd_opt=0,R0=2.1,R90=2.7)
    return my_aniso_metal

## Visualize the constitutive behavior of the materials
def tension_tests():
    iso_metal_vm().tension_test()
    iso_metal_hf8().tension_test()
    aniso_metal_hill().tension_test()

## Boundary conditions
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
    my_loading = Bnd_A(name='proportional_short',
                       alpha0=0,alpha1=1,alphan=3,beta=0,
                       sr_eq=1e-3,dbar=1e-3,ebar_mx=0.1)
    return my_loading

def prop_loading_long():
    """
    long proportional loading
    """
    my_loading = Bnd_A(name='proportional_long',
                       alpha0=0,alpha1=1,alphan=15,beta=0,
                       sr_eq=1e-3,dbar=5e-4,ebar_mx=1.5)
    return my_loading

def prop_loading_refine():
    """
    long proportional loading with a refined step size
    """
    my_loading = Bnd_A(name='proportional_refine',
                       alpha0=0,alpha1=1,alphan=5,beta=0,
                       sr_eq=1e-3,dbar=1e-4,ebar_mx=1.0)
    return my_loading

def prop_loading_debug():
    """
    short proportional loading
    """
    my_loading = Bnd_A(name='proportional_short',
                       alpha0=0.1,alpha1=1,alphan=1,beta=0,
                       sr_eq=1e-3,dbar=1e-3,ebar_mx=0.10)
    return my_loading
