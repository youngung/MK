## class to determine the constitutive description
## of subjected material.

## strain hardening description
## yielding description

import lib

class Constitutive:
    def __init__(self,f_yld,f_hrd,hashcode=None):
        """
        Arguments
        ---------
        f_yld
        f_hrd
        """

        ## Create log-file to record
        ## the evolution of this material
        self.logfn = lib.gen_tempfile(
            prefix='mk-constitutive',
            affix='log',ext='txt')

        self.f_yld = f_yld
        self.f_hrd = f_hrd

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
        self.n, self.m, self.sig, self.dsig,\
            self.dm, self.qq = self.o_hrd

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
