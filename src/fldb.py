"""
Marciniak Kuczynski model based on macro-mechanical descriptions
"""
from flda import *
import materials;reload(materials)
import pickle
from MP import progress_bar
import lib; reload(lib)

def test_save_a():
    t0 = time.time()
    mat = materials.iso_metal_vm()
    #bnd = materials.prop_loading_short()
    #bnd = materials.prop_loading_long()
    bnd = materials.prop_loading_refine()
    pickle_fn = save_a(mat, bnd)
    print 'Elapsed time:', progress_bar.convert_sec_to_string(
        time.time() - t0)

    hash_a = pickle_fn[::-1][:6][::-1]
    data = load_a(hash_code=hash_a)

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
