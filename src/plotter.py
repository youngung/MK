import matplotlib as mpl
mpl.use('Agg') ## In case X-window is not available.
import matplotlib.pyplot as plt
import os,lib
import numpy as np
from vpscyld import calc_normal
cos=np.cos
sin=np.sin
arctan2=np.arctan2

def plot_log2(fn,yfunc):
    """
    Argument
    --------
    fn
    """
    dat=np.loadtxt(fn,skiprows=2).T
    step,Da,Db,Ea,Eb,Wa,Wb,D1,D2,D3,D6,E1,E2,E3,E6,S1,S2,S6,\
        bD1,bD2,bD3,bD6,bE1,bE2,bE3,bE6,bS1,bS2,bS6,\
        f,psi,ratio = dat

    fig=plt.figure(figsize=(9,9))
    for i in xrange(9):
        fig.add_subplot(3,3,i+1)
    ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9=fig.axes

    ax1.plot(Ea,Wa,'--',label='A',zorder=10); ax1.plot(Eb,Wb,label='B')
    ax1.set_xlabel(r'$\mathrm{E^{eq}}$');ax1.set_ylabel(r'$\mathrm{w(E^{eq})}$')

    ax2.plot(E2,E1,'--',label='A',zorder=10);ax2.plot(bE2,bE1,label='B')
    ax2.set_xlabel(r'$\mathrm{E_{22}}$');ax2.set_ylabel(r'$\mathrm{E_{11}}$')

    ax3.plot(S2,S1,'--',label='A',zorder=10);ax3.plot(bS2,bS1,label='B')
    ax3.set_xlabel(r'$\mathrm{\Sigma_{22}}$');ax3.set_ylabel(r'$\mathrm{\Sigma_{11}}$')

    ax4.plot(Ea,Da,'--',label='A',zorder=10);ax4.plot(Ea,Db,label='B')
    ax4.set_xlabel(r'$\mathrm{E^{eq}}$');ax4.set_ylabel(r'$\mathrm{\dot{E}^{eq}}$')
    ax4.set_yscale("log")

    ax5.plot(Ea,Wa,'--',label='A',zorder=10); ax5.plot(Ea,Wb,label='B')
    ax5.set_xlabel(r'$\mathrm{E^{A,eq}}$');ax5.set_ylabel(r'$\mathrm{w(E^{eq})}$')

    ax6.plot(Ea,f,'-')
    ax6.set_xlabel(r'$\mathrm{E^{A,eq}}$')
    ax6.set_ylabel(r'$f$')

    ax7.plot(Ea,psi,'-')
    ax7.set_xlabel(r'$\mathrm{E^{A,eq}}$')
    ax7.set_ylabel(r'$\psi$ [Degree]')

    ax8.plot(Ea,  D1+  D2,'--',label='A')
    ax8.plot(Ea,+bD1+ bD2,'-', label='B')
    ax8.set_xlabel(r'$\mathrm{E^{A,eq}}$')
    ax8.set_ylabel(r'|$\mathrm{\dot{E}_{33}}$|')
    ax8.set_yscale('log')

    ax1.legend(loc='best'); plt.tight_layout()

    if type(yfunc).__name__=='NoneType': pass
    else:
        nstp = len(S1)
        ind_mid = int(nstp/2.)

        locus_ps_bar, locus_pi_bar = lib.y_locus_c(1000,yfunc)
        sbar_a = yfunc([S1[-1],S2[-1],0,0,0,S6[-1]])
        sbar_a0= yfunc([S1[0],S2[0],0,0,0,S6[0]])
        sbar_b = yfunc([bS1[-1],bS2[-1],0,0,0,bS6[-1]])

        sbar_am = yfunc([S1[ind_mid],S2[ind_mid],0,0,0,S6[ind_mid]])
        sbar_bm = yfunc([bS1[ind_mid],bS2[ind_mid],0,0,0,bS6[ind_mid]])

        locus_ps_a = locus_ps_bar*sbar_a
        locus_ps_am = locus_ps_bar*sbar_am
        locus_ps_bm = locus_ps_bar*sbar_bm
        locus_ps_0 = locus_ps_bar*sbar_a0
        locus_ps_b = locus_ps_bar*sbar_b

        ax9.plot(locus_ps_0[0],locus_ps_0[1],'b-',label='Initial')
        ax9.plot(locus_ps_a[0],locus_ps_a[1],'g-',label='A-final')
        ax9.plot(locus_ps_am[0],locus_ps_am[1],'g--',label='A-mid')
        ax9.plot(locus_ps_bm[0],locus_ps_bm[1],'r--',label='B-mid')
        ax9.plot(locus_ps_b[0],locus_ps_b[1],'r-',label='B-final')


        ## arrow - A final
        th = np.arctan2(D2[-1],D1[-1])
        x,y=calc_normal.find_stress(locus_ps_a[0],locus_ps_a[1],th)
        x,y=x[0],y[0]
        radi = np.sqrt(x**2+y**2)*0.1
        dx,dy = radi*cos(th),radi*sin(th)
        ax9.arrow(S1[-1],S2[-1],dx,dy)
        ax9.plot(S1[-1],S2[-1],'go')

        ## arrow - A - mid
        th = np.arctan2(D2[ind_mid],D1[ind_mid])
        x,y=calc_normal.find_stress(locus_ps_am[0],locus_ps_am[1],th)
        x,y=x[0],y[0]
        dx,dy = radi*cos(th),radi*sin(th)
        ax9.arrow(S1[ind_mid],S2[ind_mid],dx,dy)
        ax9.plot(S1[ind_mid],S2[ind_mid],'go')

        ## arrow - initial
        th = np.arctan2(D2[0],D1[0])
        x,y=calc_normal.find_stress(locus_ps_0[0],locus_ps_0[1],th)
        x,y=x[0],y[0]
        dx,dy = radi*cos(th),radi*sin(th)
        ax9.arrow(S1[0],S2[0],dx,dy)
        ax9.plot(S1[0],S2[0],'o',mfc='None',mec='k')

        ## arrow - B - mid
        th = np.arctan2(bD2[ind_mid],bD1[ind_mid])
        x,y=calc_normal.find_stress(locus_ps_bm[0],locus_ps_bm[1],th)
        x,y=x[0],y[0]
        dx,dy = radi*cos(th),radi*sin(th)
        ax9.arrow(bS1[ind_mid],bS2[ind_mid],dx,dy)
        ax9.plot(bS1[ind_mid],bS2[ind_mid],'ro')

        ## arrow - B final
        th = np.arctan2(bD2[-1],bD1[-1])
        x,y=calc_normal.find_stress(locus_ps_b[0],locus_ps_b[1],th)
        x,y=x[0],y[0]
        dx,dy = radi*cos(th),radi*sin(th)
        ax9.arrow(bS1[-1],bS2[-1],dx,dy)
        ax9.plot(bS1[-1],bS2[-1],'ro')

        mxx=max(locus_ps_a[0])
        mxy=max(locus_ps_a[1])
        mx=max([mxx,mxy])
        l0=-mx*0.2
        ax9.set_xlim(l0,)
        ax9.set_ylim(l0,)

        raise IOError, 'debug'

    _fn_='log_%s.pdf'%os.path.split(fn)[1].split('.log')[0]
    print 'saved to %s'%_fn_
    fig.savefig(_fn_,bbox_inches='tight')
    return _fn_

def plot_log(fn):
    e_a = []; e_b = []; d_a= [];  d_b=[]
    s_a = []; s_b = [];
    edeq_a=[]; edeq_b=[]; eeq_a=[];eeq_b=[]

    with open(fn,'r') as fo:
        lines = fo.readlines()[1:]
        for l in lines:
            l=l.split('\n')[0]
            stp,x0,f,df,Eeq_A,Eeq_B,Seq_A,Seq_B,\
                EdeqA,EdeqB,ratio,psi=map(float,l.split()[:12])
            E11_a,E22_a,E12_a,E33_a,D11_a,D22_a,D12_a,\
                D33_a,S11_a,S22_a,S12_a,S33_a=map(float,l.split()[13:25])
            E11_b,E22_b,E12_b,E33_b,D11_b,D22_b,D12_b,D33_b,S11_b,\
                S22_b,S12_b,S33_b=map(float,l.split()[26:])

            e_a.append([E11_a,E22_a,E12_a]); e_b.append([E11_b,E22_b,E12_b])
            d_a.append([D11_a,D22_a,D12_a]); d_b.append([D11_b,D22_b,D12_b])

            s_a.append([S11_a,S22_a,S12_a])
            s_b.append([S11_b,S22_b,S12_b])

            eeq_a.append(Eeq_A); eeq_b.append(Eeq_B); edeq_a.append(EdeqA); edeq_b.append(EdeqB)

    e_a = np.array(e_a).T; e_b = np.array(e_b).T
    d_a = np.array(d_a).T; d_b = np.array(d_b).T

    edeq_a= np.array(edeq_a); edeq_b= np.array(edeq_b)
    eeq_a= np.array(eeq_a); eeq_b= np.array(eeq_b)

    from matplotlib.backends.backend_pdf import PdfPages
    pdf_master = PdfPages('log_plotter_%s.pdf'%os.path.split(fn)[-1])

    ft = dict(fontsize=11)

    ## strain path
    fig = plt.figure(); ax=fig.add_subplot(111)
    ax.plot(e_a[1],e_a[0],label='A region')
    ax.plot(e_b[1],e_b[0],label='B region')
    ax.set_xlabel(r'$E_{11}$',ft); ax.set_ylabel(r'$E_{22}$',ft)
    ax.locator_params(nbins=4)
    fig.tight_layout()
    ax.legend();plt.close(fig);pdf_master.savefig(fig)

    ## stress path
    fig = plt.figure(); ax=fig.add_subplot(111)
    ax.plot(s_a[1],s_a[0],label='A region')
    ax.plot(s_b[1],s_b[0],label='B region')
    ax.set_xlabel(r'$S_{11}$',ft); ax.set_ylabel(r'$S_{22}$',ft)
    ax.locator_params(nbins=4)
    fig.tight_layout()
    ax.legend();plt.close(fig);pdf_master.savefig(fig)

    ## Eq. strain rate
    fig = plt.figure(); ax=fig.add_subplot(111)
    ax.plot(eeq_a,edeq_a,label='A region')
    ax.plot(eeq_b,edeq_b,label='B region')
    ax.set_xlabel(r'$E^{eq}$',ft)
    ax.set_ylabel(r'$\dot{E}^{eq}$',ft)
    ax.locator_params(nbins=4)
    fig.tight_layout()
    ax.legend();plt.close(fig);pdf_master.savefig(fig)

    ## Thickness strain rate
    fig = plt.figure(); ax=fig.add_subplot(111)
    ax.plot(eeq_a,-d_a[0]-d_a[1],label='A region')
    ax.plot(eeq_b,-d_b[0]-d_b[1],label='B region')
    ax.set_xlabel(r'$E^{eq}$',ft)
    ax.set_ylabel(r'$\dot{E}_{33}$',ft)
    ax.locator_params(nbins=4)
    fig.tight_layout()
    ax.legend();plt.close(fig);pdf_master.savefig(fig)

    pdf_master.close()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fn',type=str,help='file name (full path)')
    args        = parser.parse_args()

    plot_log2(args.fn)
    os.system('open *%s.pdf'%os.path.split(args.fn)[1].split('.log')[0])
