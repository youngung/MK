import matplotlib as mpl
mpl.use('Agg') ## In case X-window is not available.
import matplotlib.pyplot as plt
import os
import numpy as np


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

    plot_log(args.fn)
