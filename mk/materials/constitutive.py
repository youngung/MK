## class to determine the constitutive description
## of subjected material.

## strain hardening description
## yielding description

import mk.library.lib as lib

"""
Snapshot may be more useful by allowing 'custom' printing
"""

class Snapshot:
    def __init__(self):
        self.logfn = lib.gen_tempfile(
            prefix='mk-constitutive',
            affix='log',ext='txt')
    def takeshot(self,**kwargs):
        """
        Take a snapshot of given key-worded arguments
        """
        with open(self.logfn,'a') as fo:
            for key, value in kwargs.iteritems():
                fo.write('%11.4e '%float(value))

    def linebreak(self):
        """
        Insert linebreak in the snapshot file
        """
        with open(self.logfn,'a') as fo:
            fo.write('\n')

class Constitutive:
    def __init__(self,f_yld=None,f_hrd=None,
                 hashcode=None,
                 params_yld=None, label_yld=None,
                 params_hrd=None, label_hrd=None):
        """
        Arguments
        ---------
        f_yld
        f_hrd
        hashcode
        params_yld
        label_yld
        params_hrd
        label_hrd
        """

        ## Create log-file to record
        ## the evolution of this material
        self.logfn = lib.gen_tempfile(
            prefix='mk-constitutive',
            affix='log',ext='txt')

        self.f_yld = f_yld
        self.f_hrd = f_hrd

        ## used when constructing <fortran> objects
        ## for hardening and yield function calculations
        ## -- see self.set_hrd and self.set_yld
        self.params_yld = params_yld
        self.params_hrd = params_hrd
        self.label_yld  = label_yld
        self.label_hrd  = label_hrd

    def set_hrd(self):
        """
        Manually set yield functions
        """
        label=self.label_hrd.lower()
        params=self.params_hrd
        m,qq=5e-2, 1e3
        if label[:4]=='voce':
            a,b0,c,b1 = self.params_hrd
            import mk.materials.func_hard_for
            self.f_hrd = mk.materials.func_hard_for.return_voce(
                a=a,b0=b0,c=c,b1=b1,m=m,qq=qq)
        elif label[:5]=='swift':
            raise IOError, 'swift is not supported yet.'
        else:
            raise IOError, 'Unexpected label found %s'%self.label_hrd

    def set_yld(self):
        """
        Manually set yield functions
        """
        label=self.label_yld.lower()
        params=self.params_yld
        import mk.yieldFunction.tuneYld2000, mk.yieldFunction.yf2
        if label=='h48r':
            yldFunc=mk.yieldFunction.yf2.wrapHill48R(params)
        elif label=='yld1':
            ys=params[1][::]
            rv=params[0][::]
            rb,yb=mk.yieldFunction.tuneYld2000.H48toYld_withYS(
                rv=rv,ys=ys,m=6,iopt=1)
            rv.append(rb)
            ys.append(yb)
            yldFunc = mk.yieldFunction.yf2.wrapYLD(
                r=rv,y=ys,m=6,k=2)
        elif label=='yld2':
            ys=params[1][::]
            rv=params[0][::]
            yldFunc = mk.yieldFunction.yf2.wrapYLD(
                r=rv,y=ys,m=6,k=2)
        else:
            raise IOError, 'Unexpected label for yield function %s %s'%(self.label_yld,label)

        self.f_yld = yldFunc

    def update_yld(self,stress):
        """
        Given the stress (and potentially
        many more state variables not included yet),
        update the results of yield function

        Argument
        --------
        stress
        """
        self.o_yld = self.f_yld(stress)
        self.stress, self.phi, self.dphi,\
            self.d2phi = self.o_yld

    def update_hrd(self,strain):
        """
        Given the strain (and potentially
        many more state variables not included yet),
        update the results of yield function

        Argument
        --------
        strain
        """
        self.eps = strain
        self.o_hrd = self.f_hrd(self.eps)

        self.m    = self.o_hrd[0]
        self.sig  = self.o_hrd[1]
        self.dsig = self.o_hrd[2]
        self.dm   = self.o_hrd[3]
        self.qq   = self.o_hrd[4]

    def recordCurrentStat(self):
        """
        Record the current status of
        this material to self.logfn
        """
        fmt1 = '%8.4f %8.1f' ## Effective strain-stress (hardening)
        fmt2 = '%8.3f' ## 6D stress
        fmt3 = '%8.3f' ## 6D strain

        with open(self.logfn,'a') as fo:
            ## hardening (effective strain/effective stress)
            fo.write(fmt1%(self.eps,self.sig))
            ## stress states on the yield locus phi of 1.
            for i in xrange(6):
                fo.write(fmt2%self.stress[i])
            ## derivatives of plastic potential
            for i in xrange(6):
                fo.write(fmt3%self.dphi[i])
            fo.write('\n') ## line breaker
