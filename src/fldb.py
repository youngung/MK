"""
Marciniak Kuczynski model based on macro-mechanical descriptions
"""
from flda import *
import materials
# try:
#     import cPickle as pickle
# except:
import pickle
from MP import progress_bar
import lib; reload(lib)

def test_save_a():
    t0 = time.time()
    mat = materials.iso_metal_vm()
    #bnd = materials.prop_loading_short()
    bnd = materials.prop_loading_long()
    pickle_fn = save_a(mat, bnd)
    print 'Elapsed time:', progress_bar.convert_sec_to_string(
        time.time() - t0)

    hash_a = pickle_fn[::-1][:6][::-1]
    data = load_a(hash_code=hash_a)


def test_main():
    pass

def save_a(mat, bnd):
    """
    For initial, save FLA results.

    Arguments
    ---------
    mat (material class)
    bnd (boundary class)
    """
    data_collection_FLDA = []
    for i in xrange(len(bnd.alphas)):
        rst = FLDA_onepath(
            ## Boundary condition
            bnd.alphas[i],bnd.beta,bnd.sr_eq,
            bnd.dbar,bnd.ebar_mx,

            ## material
            mat.func_yd,mat.func_hd,mat.func_sr)

        data_collection_FLDA.append(rst)

    ## Archive the results.
    hash_a = lib.gen_hash_code(nchar=6)
    d_asc,t_asc=get_stime()
    pickle_fn = 'FLDA_MACRO_%s_%s_%s'%(d_asc,t_asc,hash_a)
    with open(pickle_fn,'wb') as fo:
        # pickle.dump(mat,fo,pickle.HIGHEST_PROTOCOL)
        # pickle.dump(bnd,fo,pickle.HIGHEST_PROTOCOL)
        ## I'd eventually need a method to save
        ## the above two commented-out lines
        ## for the sake of completion
        pickle.dump(data_collection_FLDA,fo,pickle.HIGHEST_PROTOCOL)

    print 'Results are pickled to %s'%pickle_fn
    return pickle_fn

def load_a(fn=None,hash_code=None):
    """
    Arguments
    ---------
    hash_code
    """
    import os, glob
    p=os.popen('ls FLDA_MACRO_*_*_%s'%hash_code)
    elems = p.readlines()
    if len(elems)!=1:
        raise IOError, 'Could not find the file'
    fn = elems[0].split('\n')[0]
    with open(fn,'rb') as fi:
        data_collection_FLDA = pickle.load(fi)
    print '%s was restored'%fn
    return data_collection_FLDA




def main():
    import modules
    from scipy import optimize


    t0 = time.time()
    ## Choice of material / loading condition
    mat = materials.iso_metal_vm()
    """
    mat.func_hd as a function of ebar
    mat.func_sr as a function of ebar dot
    mat.func_yd as a function of cauchy stress
    """
    bnd = materials.prop_loading_long()

    ## Save FLDa
    FLDA_pck_name = save_a(mat, bnd)
    hash_a = FLDA_pck_name[::-1][:6][::-1]

    ## Load from pickles.
    data_collection_FLDA = load_a(hash_code=hash_a)

    ## Boundary condition
    beta    = bnd.beta
    sr_eq   = bnd.sr_eq
    dbar    = bnd.dbar
    ebar_mx = bnd.ebar_mx


    ## Marciniak-Kuczynski parameters
    f  = 0.990
    psi= 0.
    n  = [1.,0.,0.]
    t  = [0.,1.,0.]


    for i in xrange(bnd.nprob):
        if i==0:
            eps_b_eq = 0.

        TT0 = time.time()
        ## Load from pickles.
        eps6,sig6,ebar,sbar,\
            sr6,times = data_collection_FLDA[i]

        t = 0.
        for j in xrange(len(eps6)):
            ## Find response of region B.
            eps6_a = eps6[j,:]
            sig6_a = sig6[j,:]
            ebar_a = ebar[j]
            sbar_a = sbar[j]
            sr6_a = sr6[j,:]
            dt=times[j]-t

            sig33_a = s62c(sig6_a)
            sig33_a_grv = rot_tensor(sig33_a,psi)
            sr33_a = s62c(sr6_a)
            sr33_a_grv = rot_tensor(sr33_a,psi)

            sig33_b_grv = np.zeros((3,3))
            sig33_b_grv = sig33_a_grv/f
            sig6_b_grv  = c2s6(sig33_b_grv)

            sr33_b_grv = np.zeros((3,3))
            sr33_b_grv[0,0] = 0. ## unknown
            sr33_b_grv[1,1] = sr33_a_grv[1,1]
            sr33_b_grv[2,2] = 0. ## unknown
            sr33_b_grv[0,1] = sr33_a_grv[0,1]
            sr33_b_grv[1,0] = sr33_a_grv[1,0]
            sr6_b_grv = c2s6(sr33_b_grv)

            ## Find e11_dot
            objf = modules.find_e11_dot(
                sig6_b_grv,mat.func_yd,sr6_b_grv,
                eps_b_eq,dt,mat.func_sr,mat.func_hd)

            rst = optimize.minimize(objf,x0=1.e-3)


            ## updates (Increase by time incremental step)
            print type(rst)
            print 'e11_dot:', rst.x
            print 'fun:    ', rst.fun
            return 


        print 'Elapsed time for #%i loading: %s'%(
            i+1,
            progress_bar.convert_sec_to_string(
                time.time() - TT0))

    print 'Total Elapsed time for main:',\
        progress_bar.convert_sec_to_string(
            time.time() - t0)
