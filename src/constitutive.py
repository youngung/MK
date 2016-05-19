## class to determine the constitutive description
## of subjected material.

## strain hardening description
## yielding description

class Constitutive:
    def __init__(self,f_yld,f_hrd):
        """
        Arguments
        ---------
        f_yld
        f_hrd
        """
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
        self.strain = strain
        self.o_hrd = self.f_hrd(self.strain)
        self.n, self.m, self.sig, self.dsig,\
            self.dm, self.qq = self.o_hrd
