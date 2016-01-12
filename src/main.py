"""
"""
def an_increment_a(rho=0.,er=1e-3,mat=None):
    """
    Deform the material an incremental step and
    return the stress tensor.
    """
    ## Material choice.
    ## if mat is None, assume isotropic mat.
    if type(mat)==None:
        mat='iso'
