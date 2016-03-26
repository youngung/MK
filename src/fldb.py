"""
Marciniak Kuczynski model based on macro-mechanical descriptions
"""
from flda import *
import materials;reload(materials)
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


def main(f0=0.999,psi0=0):
    """
    Arguments
    ---------
    f0
    psi0
    """
    # np.seterr(all='raise')

    import modules; reload(modules)
    from scipy import optimize

    t0 = time.time()
    ## Choice of material / loading condition
    # mat = materials.iso_metal_vm()
    mat = materials.iso_metal_hf8()
    """
    mat.func_hd as a function of ebar
    mat.func_sr as a function of ebar dot
    mat.func_yd as a function of cauchy stress
    """
    #bnd = materials.prop_loading_long()
    bnd = materials.prop_loading_refine()

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
    n  = [1.,0.,0.];   t  = [0.,1.,0.]
    f  = f0; psi = psi0

    rot_mat = rot(psi)
    n = lib.rot_vec(rot_mat,n)
    t = lib.rot_vec(rot_mat,t)

    fmt='%7i %7.4f %7.4f %7.1e %7.4f %7.4f %7.1f'+\
        ' %7.1f %7.4f %7.4f %7.2f'
    print ('%7s '*11)%(
        'stp','x0','f','df',
        'Eeq_A','Eeq_B','Seq_A','Seq_B',
        'EdeqA','EdeqB','ratio')

    mx_dp = 10000
    sig_A = np.zeros((bnd.nprob,6,mx_dp))*np.nan
    eps_A = np.zeros((bnd.nprob,6,mx_dp))*np.nan
    edt_A = np.zeros((bnd.nprob,  mx_dp))*np.nan
    Seq_A = np.zeros((bnd.nprob  ,mx_dp))*np.nan
    Eeq_A = np.zeros((bnd.nprob  ,mx_dp))*np.nan
    sig_B = np.zeros((bnd.nprob,6,mx_dp))*np.nan
    eps_B = np.zeros((bnd.nprob,6,mx_dp))*np.nan
    edt_B = np.zeros((bnd.nprob  ,mx_dp))*np.nan
    Seq_B = np.zeros((bnd.nprob  ,mx_dp))*np.nan
    Eeq_B = np.zeros((bnd.nprob  ,mx_dp))*np.nan

    times = np.zeros((bnd.nprob))*np.nan

    for i in xrange(bnd.nprob):
        eps_b_eq = 0.
        TT0 = time.time()
        ## Load from pickles.
        eps6,sig6,ebar,sbar,sr6,times \
            = data_collection_FLDA[i]

        t  =  0.
        x0 =  1e-5
        eps6_b=np.zeros((6,))
        # print 'nit, f1, f2, jac, x0, objf(x0)'
        plt.ioff()
        for j in xrange(len(eps6)):
            ## Find response of region B.
            eps6_a   = eps6[j,:]
            sig6_a   = sig6[j,:]
            eps_a_eq = ebar[j]
            sbar_a   = sbar[j]
            sr6_a    = sr6[j,:]
            dt       = times[j]-t

            sig33_a  = s62c(sig6_a)
            sr33_a   = s62c(sr6_a)
            sig33_a_grv = rot_tensor(sig33_a,psi)
            sr33_a_grv  = rot_tensor(sr33_a ,psi)

            sig33_b_grv = sig33_a_grv/f
            sig6_b_grv  = c2s6(sig33_b_grv)

            sr33_b_grv      = np.zeros((3,3))
            sr33_b_grv[0,0] =-1e-3 ## unknown
            sr33_b_grv[1,1] = sr33_a_grv[1,1]
            sr33_b_grv[2,2] = 1e-3 ## unknown but linked with sr33_b_grv[0,0]
            sr33_b_grv[0,1] = sr33_a_grv[0,1]
            sr33_b_grv[1,0] = sr33_a_grv[1,0]
            sr6_b_grv       = c2s6(sr33_b_grv)

            # print 'sr12_b_grv:',sr33_b_grv[0,1],sr33_b_grv[1,0]
            # print 'sr12_a_grv:',sr33_a_grv[0,1],sr33_a_grv[1,0]
            # return

            ## Find e11_dot
            objf = modules.find_e11_dot(
                sig6_b_grv,mat.func_yd,sr6_b_grv,
                eps_b_eq,dt,mat.func_sr,mat.func_hd)

            ## N/R
            irepeat = True; nit = 0; nit_mx = 100
            err_tol = 1.5e-2
            dx      = 1.e-6 ## debug?
            verbose = True
            if verbose and np.mod(j,20)==0:
                print '%4s %4s %4s %8s %8s %8s %8s'%(
                    'nit','f1','f2','jac','x0','objf','eb_eq')

            while (irepeat):
                f   = objf(x0)

                try:
                    ## Jacobian as middle value
                    f1  = objf(x0-dx)
                    f2  = objf(x0+dx)
                    jac = (f2-f1)/(dx*2)

                except:
                    ## Jacobian from forward direction
                    f1  = objf(x0+dx)
                    jac = (f1-f)/dx

                x0  = x0 - f/jac

                if verbose:
                    print '%4i %4.1f %4.1f %8.1e %8.5f %8.5f %8.5f'%(
                        nit, f1, f2, jac, x0, f, eps_b_eq)

                if abs(f) < err_tol:
                    sr33_b_grv[0,0] = x0
                    sr33_b_grv[2,2] = - sr33_b_grv[0,0] - sr33_b_grv[1,1]
                    irepeat         = False
                if nit>=nit_mx:
                    raise IOError, 'could not find the solution..'
                else:
                    nit=nit+1

            ## updates (Increase by time incremental step)
            ## Complete strain increment for region B:
            sr6_b_grv = c2s6(sr33_b_grv)

            ## rotate the strain rate back
            sr33_b  = rot_tensor(sr33_b_grv, -psi)
            sig33_b = rot_tensor(sig33_b_grv,-psi)
            sr6_b   = c2s6(sr33_b)
            sig6_b  = c2s6(sig33_b)

            ## work rate for B:
            sigb_eq = mat.func_yd(sig6_b)
            siga_eq = mat.func_yd(sig6_a)
            wrate_b = np.dot(sr6_b,sig6_b)
            eps_b_eq_dot = wrate_b / sigb_eq
            wrate_a = np.dot(sr6_a,sig6_a)
            eps_a_eq_dot = wrate_a / sigb_eq

            ## equivalent strain update for region B
            eps_b_eq = eps_b_eq + eps_b_eq_dot * dt
            eps6_b   = eps6_b + sr6_b*dt

            sig_A[i,:,j] = sig6_a[::]
            eps_A[i,:,j] = eps6_a[::]
            Seq_A[i,j]   = siga_eq
            Eeq_A[i,j]   = eps_a_eq
            edt_A[i,j]   = eps_a_eq_dot

            sig_B[i,:,j] = sig6_b[::]
            eps_B[i,:,j] = eps6_b[::]
            Seq_B[i,j]   = sigb_eq
            Eeq_B[i,j]   = eps_b_eq
            edt_B[i,j]   = eps_b_eq_dot

            times[i] = t

            t = t + dt

            ## psi update
            ## n update, t update
            ## f update
            new_f = f0 * np.exp(eps6_b[2] - eps6_a[2])
            df    = new_f - f
            f     = f + df
            ratio = sr6_b[2] / sr6_a[2]

            ## check the thickness strain rate difference.
            ## write values.
            if ratio>10:
                print fmt%(j+1,x0,f,df,eps_a_eq,eps_b_eq,
                           siga_eq,sigb_eq,eps_a_eq_dot,eps_b_eq_dot,
                           ratio)
                print 'Congratulations, you just failed'
                x0 = 1e-4
                break

        print 'Elapsed time for #%i loading: %s'%(
            i+1,
            progress_bar.convert_sec_to_string(
                time.time() - TT0))

    # ## Data vis
    fig=plt.figure(figsize=(6,7))
    plt.ioff()
    ax1=fig.add_subplot(221); ax2=fig.add_subplot(222)
    ax3=fig.add_subplot(223); ax4=fig.add_subplot(224)

    for i in xrange(bnd.nprob):
        ax1.plot(eps_A[i,1,:],eps_A[i,0,:],'b-')
        ax1.plot(eps_B[i,1,:],eps_B[i,0,:],'r-')
        ax2.plot(sig_A[i,1,:],sig_A[i,0,:],'b-')
        ax2.plot(sig_B[i,1,:],sig_B[i,0,:],'r-')
        l, = ax3.plot(Eeq_A[i,:], Seq_A[i,:],ls='-')
        ax3.plot(Eeq_B[i,:], Seq_B[i,:],color=l.get_color(),ls='--')

    print '-'*90
    print 'Congratulations, your sheet SUCCESSFULY failed'
    print '-'*90

    print 'Total Elapsed time for main:',\
        progress_bar.convert_sec_to_string(
            time.time() - t0)

    draw_guide(ax1)
    draw_guide(ax2,r_line=[0,0.5,1],max_r=1000)

    ax1.set_aspect('equal')
    ax2.set_ylim(0.,);ax3.set_ylim(0.,)
    ax2.set_aspect('equal')

    ax3.set_xlabel(r'$\mathrm{\bar{E}^{eq}_{22}}$',dict(fontsize=20))
    ax3.set_ylabel(r'$\mathrm{\bar{\Sigma}^{eq}_{22}}$',dict(fontsize=20))

    ax1.set_xlabel(r'$\mathrm{\bar{E}_{22}}$',dict(fontsize=20))
    ax1.set_ylabel(r'$\mathrm{\bar{E}_{11}}$',dict(fontsize=20))
    ax2.set_xlabel(r'$\mathrm{\bar{\Sigma}_{22}}$',dict(fontsize=20))
    ax2.set_ylabel(r'$\mathrm{\bar{\Sigma}_{11}}$',dict(fontsize=20))

    plt.ion();
    plt.draw()
    plt.tight_layout()
    fig.savefig('FLD_results.pdf',bbox_inches='tight')
