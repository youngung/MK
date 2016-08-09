"""
"""
import os, glob, subprocess, tarfile, dill, shutil,MP, MP.lib.temp
import numpy as np
import matplotlib.pyplot as plt

def syncArchiveFromPal():
    """
    Sync archive files generated from Pal
    """
    fnRsync=MP.lib.temp.gen_tempfile()
    # cwd=os.getcwd()
    cwd = '~/repo/mk/archive'
    runRsync="""#!/bin/bash
    cd %s
    bash rsyncFromPal.sh
    cd -"""%cwd
    print runRsync
    with open(fnRsync,'w') as fo:
        fo.write(runRsync)
        fo.write('\n')
    rst=subprocess.check_output(['bash',fnRsync])

def loadFLDmin(fn='minFLD.txt'):
    dat=np.loadtxt(fn,dtype='str').T
    if len(dat)==0:
        ## no results found in minFLD.txt
        a=np.zeros((9,1))
        a[::]=np.nan
    else:
        a=np.array(dat[:9],dtype='float')
    # exx,eyy,thi,thf,sx,sy,sigbar,epsbar,dt=dat[:9]
    return a

class member:
    def __init__(self,fn=None):
        """
        Argument
        --------
        fn
        """
        self.fn = fn
        ## check the rule.
        eps_eq, yflab, hflab, path, f0 =  self.fn.split('_')[1:6]

        self.eps_eq = float(eps_eq)
        self.yflab  = yflab
        self.hflab  = hflab
        self.hardpath = path
        self.f0       = float(f0)

        ##
        cmd = ['tar','-tf',self.fn]
        self.members=subprocess.check_output(cmd).split('\n')[:-1]
        self.members.sort()
        self.fldfn=self.members[1]

        ## find out min file
        cmd=['tar','-xf',self.fn, self.fldfn]
        if subprocess.check_call(cmd)!=0:
            raise IOError
        exx,eyy,thi,thf,sx,sy,sigbar,epsbar,dt = loadFLDmin(self.fldfn)

        ths = np.arctan2(eyy,exx)##TD,RD
        rad = np.sqrt(exx**2+eyy**2)

        os.remove(self.fldfn)
        self.flc = np.array([exx,eyy])

class EachMK: ## for each MKVPSCHARD file (fixed f0 value)
    def __init__(self,fn):
        """
        Argument
        --------
        fn
        """
        cmd=['tar','-tf',fn]
        self.fns_member = subprocess.check_output(cmd).split()
        self.members=[]
        self.nmem = len(self.fns_member)
        self.eps_eqs = []
        self.yflabs  = []
        self.hflabs  = []
        self.hardpaths = []
        self.f0s       = []
        inds = np.zeros((len(self.fns_member),4))
        for i in xrange(len(self.fns_member)): ## 1d to (w,x,y,z)
            subprocess.check_call(['tar','-xf',fn,self.fns_member[i]])
            self.members.append(member(fn = self.fns_member[i]))
            if not(self.members[i].eps_eq in self.eps_eqs):
                self.eps_eqs.append(self.members[i].eps_eq)
            if not(self.members[i].yflab in self.yflabs):
                self.yflabs.append(self.members[i].yflab)
            if not(self.members[i].hflab in self.hflabs):
                self.hflabs.append(self.members[i].hflab)
            if not(self.members[i].hardpath in self.hardpaths):
                self.hardpaths.append(self.members[i].hardpath)
            if not(self.members[i].f0 in self.f0s):
                self.f0s.append(self.members[i].f0)
            os.remove(self.fns_member[i])

        w, x, y, z = self.yflabs, self.hflabs, self.hardpaths, self.eps_eqs
        self.find_dimension()

    def find_dimension(self):
        """
        Return the dimension of
        (YldFunc, Eps_eq, HrdFunc, HrdPath)
        """
        self.dimension = (len(self.yflabs),len(self.eps_eqs), len(self.hflabs), len(self.hardpaths))

    def find_mem(self,iyf,ihf,ipath,ieps,verbose):
        """
        Find the member of the tar files
        that correspondes to the metric given in terms of
        <iyf>, <ihf>, <ipath>, and <ieps>

        Arguments
        ---------
        iyf
        ihf
        ipath
        ieps
        verbose
        """
        yfl = self.yflabs[iyf]
        hfl = self.hflabs[ihf]
        hpt = self.hardpaths[ipath]
        eps = self.eps_eqs[ieps]

        if verbose:
            print 'yf lab:', yfl
            print 'hf lab:', hfl
            print 'hpt   :', hpt
            print 'eps   :', eps

        ind=-1
        for i in xrange(len(self.members)):
            eachMem=self.members[i]
            if eachMem.yflab==yfl and eachMem.hflab==hfl \
               and eachMem.hardpath==hpt and eachMem.eps_eq==eps:
                ind=i;break
        if ind==-1: raise IOError
        return ind, yfl, hfl, hpt, eps, eachMem

class master:
    def __init__(self,fns,fn_exp='/Users/yj/repo/vpsc/vpsc-dev-fld/ipynb/FLD/IFsteel_EXP/FLDexp.dill'):
        """
        Arguments
        ---------
        fns    = 'tar files for different values of f0.'
        fn_exp = 'Experimental data'
        """
        self.eachMKs=[]
        for i in xrange(len(fns)):
            self.eachMKs.append(EachMK(fns[i]))
        self.eachMKs=np.array(self.eachMKs)

        if len(self.eachMKs)==0:
            raise IOError, 'no MK data found.'

        self.get_exp(fn_exp)

        collectionMK = self.compareFit_ind(iyf=0,ihf=0,ipath=0,ieps=0,iplot=False)
        self.f0s=[]
        for i in xrange(len(collectionMK)):
            self.f0s.append(collectionMK[i].f0)
        self.find_dimension()

    def get_exp(self,fn):
        """
        Read experimental data.

        Argument
        --------
        fn
        """
        with open(fn,'r') as fo:
            self.exp_FLC=dill.load(fo)

    def find_dimension(self):
        """
        Find the dimension of the problems..
        """
        self.dimension = []
        self.dimension.append(len(self.f0s)) # nf0
        for i in xrange(len(self.eachMKs[0].dimension)):
            ## yf, eps, hflabs, hrdPath
            self.dimension.append(self.eachMKs[0].dimension[i])
        ## construct self.dimension
        ## f0, yf, eps, hflabs, hrdPath
        self.dimension=np.array(self.dimension)
        nf0, nyf, neps, nhf, npth = self.dimension
        self.membersContainer = np.empty((nf0,nhf,npth,nyf,neps),dtype='object')
        for ihf in xrange(nhf):
            for ipt in xrange(npth):
                for iyf in xrange(nyf):
                    for iep in xrange(neps):
                        cMK=self.compareFit_ind(iyf,ihf,ipt,iep,False)
                        for i in xrange(nf0):
                            self.membersContainer[i,ihf,ipt,iyf,iep]=cMK[i]

    def plot_all(self,label=None):
        """
        Plot all of the data...

        Argument
        --------
        label
        """
        from matplotlib.backends.backend_pdf import PdfPages
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        import MP.lib.mpl_lib, MP.lib.axes_label

        if type(label).__name__=='NoneType':
            fn='masterMKdata.pdf'
        elif type(label).__name__=='NoneType':
            fn='masterMKdata_%s.pdf'%label
        figPages=PdfPages(fn)

        ## Loop over hardening function
        ## Loop over path, based on which hardening parameters are calibrated
        nf0, nyf, neps, nhf, npth = self.dimension
        # membersContainer = np.empty((nf0,nhf,npth,nyf),dtype='object')
        for ihf in xrange(nhf): ## hardening function
            for ipt in xrange(npth): ## strain path
                fig=plt.figure(figsize=(3.5*nyf, 3.0*neps))
                grid=gridspec.GridSpec(neps,nyf,wspace=0,hspace=0) ## col/row
                axes_grid=[]

                for iyf in xrange(nyf):
                    axes_grid.append([])
                    for iep in xrange(neps):
                        ax=fig.add_subplot(grid[iep,iyf])
                        axes_grid[iyf].append(ax)
                        # cMK=self.compareFit_ind(iyf,ihf,ipt,iep,False)
                        for i in xrange(nf0):
                            cMK=self.membersContainer[i,ihf,ipt,iyf,iep]
                            ax.plot(cMK.flc[0],cMK.flc[1],label=cMK.f0)
                        y,x,yerr,xerr = self.exp_FLC[:,0],self.exp_FLC[:,1],\
                                        self.exp_FLC[:,2],self.exp_FLC[:,3]
                        ax.errorbar(x=x,y=y,xerr=xerr,yerr=yerr,marker='None',
                                    ls='None',color='k',label='Exp ISO')
                        ax.legend(loc='lower left')
                        MP.lib.axes_label.draw_guide(ax,r_line=[-0.5,0,1,2,2.5])
                        ax.set_xlim(-0.5,1); ax.set_ylim(-0.5,1)

                    MP.lib.mpl_lib.rm_all_lab(fig.axes)
                    # MP.lib.mpl_lib.tune_xy_lim(fig.axes)

                axes_grid = np.array(axes_grid)
                ax_deco   = axes_grid[0,neps-1]
                MP.lib.axes_label.deco_fld(ax=ax_deco,iopt=4,iasp=False)

                ## yield function labels
                for i in xrange(nyf): ## same strain offsets.
                    ax = axes_grid[i,0]
                    txt = self.membersContainer[0,ihf,ipt,i,0].yflab
                    ax.text(0.5,1.2,txt,transform=ax.transAxes,va='center',ha='center')
                for i in xrange(neps):
                    ax = axes_grid[0,i]
                    txt = self.membersContainer[0,ihf,ipt,0,i].eps_eq
                    txt = r'$\mathrm{\bar{E}^{eq}}$=%.3f'%txt
                    ax.text(-0.3,0.5,txt,transform=ax.transAxes,
                            rotation=90,
                            va='center',ha='center')
                # axes_grid[:,0] ## same yield functions

                ## figure text
                hpath = self.membersContainer[0,ihf,ipt,0,0].hardpath
                if hpath.lower()=='u':
                    hpath='Uniaxial Tension RD'
                elif hpath.lower()=='b':
                    hpath='Bulge test'
                txt=r'%s hardening function tuned by %s'%(
                    self.membersContainer[0,ihf,ipt,0,0].hflab,hpath)

                fig.text(0.5, 1.1, txt, transform=fig.transFigure,
                         va='center',ha='center',fontsize=15)

                figPages.savefig(fig,bbox_inches='tight')
                fig.clf()
                plt.close(fig)

        print '%s has been saved'%fn
        figPages.close()

    def compareFit_ind(self,iyf=0,ihf=0,ipath=0,ieps=0,iplot=False):
        """
        Arguments
        ---------
        iyf
        ihf
        ipath
        ieps
        iplot
        """
        try: self.dimension
        except:
            iDim = False
        else:
            iDim = True

        if iDim:
            sectionRequested = np.array([iyf,ieps,ihf,ipath])
            if (sectionRequested - self.dimension[1:]>0).any():
                print '-'*20
                print 'if0:',if0
                print sectionRequested, sectionRequested.shape
                print self.dimension[if0], self.dimension[if0].shape
                raise IOError, 'The given dimension exceeds what is available'

        inds=[]
        collectionMK=[]; collectionMK_fn=[]; collectionMems=[]
        for i in xrange(len(self.eachMKs)): ## each f0s
            each = self.eachMKs[i]
            verbose=False;#verbose=True
            ind, yfl, hfl, hpt, eps, mem = each.find_mem(
                iyf,ihf,ipath,ieps,verbose)
            self.eachMKs[i].members[ind].flc
            inds.append(ind)
            collectionMK.append(self.eachMKs[i].members[ind])
            collectionMK_fn.append(self.eachMKs[i].members[ind].fldfn)

        if iplot:
            fig=plt.figure(); ax=fig.add_subplot(111)
            y,x,yerr,xerr = self.exp_FLC[:,0],self.exp_FLC[:,1],\
                            self.exp_FLC[:,2],self.exp_FLC[:,3]
            ax.errorbar(x=x,y=y,xerr=xerr,yerr=yerr,marker='None',ls='None')
            ax.set_aspect('equal')
            ax.set_xlabel(r'$\mathrm{\bar{E}_{22}}$',fontsize=14)
            ax.set_ylabel(r'$\mathrm{\bar{E}_{11}}$',fontsize=14)
            ax.set_title('%s %s %s %.2f'%(yfl,hfl,hpt,eps))

            for i in xrange(len(collectionMK)):
                dat=collectionMK[i].flc
                f0 = collectionMK[i].f0
                ax.plot(dat[0],dat[1],label=f0)
                print 'fld file name:',collectionMK_fn[i]

            ax.grid('on')
            ax.legend()
        return collectionMK
