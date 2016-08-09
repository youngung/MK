"""
Run mk_run using serialized Python functions
that represents the yield function and strain hardening function

The serialization of the Python function was enabled by the use of
dill package, available in "https://github.com/uqfoundation/dill".
"""


def yf_plot(fnpickle=None,fnhrd_vpsc=None,hashcode='aaaaa'):
    """
    Arguments
    ---------
    fnpickle
    fnhrd_vpsc
    hashcode='aaaaa'
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import mk.tests.mechtests

    with open(fnpickle, 'rb') as fo:
        eps_eq     = dill.load(fo)
        yfs        = dill.load(fo)
        yfs_labels = dill.load(fo)
        hfs        = dill.load(fo)
        hfs_labels = dill.load(fo)
        results    = dill.load(fo) ## not important yet...

    with open(fnhrd_vpsc,'rb') as fo:
        WorkContour_vpsc_rho = dill.load(fo)
        rhops = dill.load(fo)

    nfs = len(hfs)
    nys = len(yfs[0])
    fig=plt.figure(figsize=(3.6*2, 3.))
    ax1=fig.add_subplot(121)
    ax2=fig.add_subplot(122)

    for j in xrange(nys):
        ## plot yield surface
        func = yfs[0][j]
        X, Y = mk.tests.mechtests.locus(func=func, nth=100)
        ax1.plot(X,Y,label=yfs_labels[j])

    for i in xrange(nfs):
        xs=np.linspace(0,1.0)
        ax2.plot(xs,hfs[i](xs),label=hfs_labels[i])

    ax1.set_xlim(0.,); ax1.set_ylim(0.,)
    ax1.legend(); ax2.legend()
    fn='yf_hf_%s.pdf'%hashcode
    fig.savefig(fn)
    print fn, ' saved.'

## example:
"""
-- Palmetto
$ python mk_run_pickles.py --f0 0.995 --r0 0    --r1 1 --nr 4 --fnpickle /home/younguj/repo/mk/matDatabase/IFsteel/yf_collection.dll --dry
$ python mk_run_pickles.py --f0 0.995 --r0 -0.5 --r1 1 --nr 4 --fnpickle /home/younguj/repo/mk/matDatabase/IFsteel/yf_collection.dll  --fnpickle_vpsc_hard /home/younguj/repo/mk/matDatabase/IFsteel/vpsc_FC_rhop.dll --dry


$ python mk_run_pickles.py --f0 0.995 --r0 -0.5 --r1 1 --nr 4 --fnpickle ~/repo/mk/matDatabase/IFsteel/yf_collection.dll  --fnpickle_vpsc_hard ~/repo/mk/matDatabase/IFsteel/vpsc_FC_rhop.dll --dry
"""
if __name__=='__main__':
    import numpy as np
    from mk.library.lib import gen_tempfile,rhos2ths
    import os,multiprocessing,time,argparse,dill,shutil,pickle,subprocess
    from MP import progress_bar
    from mk_paths import findCorrectPsi
    from MP.lib import etc

    parser = argparse.ArgumentParser()
    parser.add_argument('--f0',type=float,help='f0 value')
    parser.add_argument('--r0',type=float,help='rho start')
    parser.add_argument('--r1',type=float,help='rho end')
    parser.add_argument('--nr',type=int,help='number of rhos')
    parser.add_argument('--fnpickle', type=str,default=None,help='Pickle file name')
    parser.add_argument(
        '--fnpickle_vpsc_hard',type=str,default=None,
        help='pickle file name on the collection of VPSC-based hardening functions')
    parser.add_argument('--dry',dest='dry',action='store_true', default=False)
    args        = parser.parse_args()

    ## using VPSC-tuned strain-hardening parameters for each strain path
    ivpsc_hard=False
    if type(args.fnpickle_vpsc_hard).__name__!='NoneType':
        ivpsc_hard=True

    with open(args.fnpickle, 'rb') as fo:
        eps_eq     = dill.load(fo)
        yfs        = dill.load(fo)
        yfs_labels = dill.load(fo)
        hfs        = dill.load(fo)
        hfs_labels = dill.load(fo)
        results    = dill.load(fo) ## not important yet...
        yfs_params = dill.load(fo)
        hfs_params = dill.load(fo)

        # print 'yfs_params'
        # print yfs_params
        # print 'hfs_labels'
        # print hfs_labels
        # os._exit(0)

    fnCollect=[]

    if ivpsc_hard:
        ## over-write hfs
        nfs = 1
        hfs_labels=['VpscHard_FIT']
    else:
        nfs = len(hfs)
        pass

    # yf.plot(args.fnpickle, fnhrd_vpsc=args.fnpickle_vpsc_hard,
    #      hashcode='aaaaa')
    # os._exit(0)

    neps_eq = len(eps_eq)
    nyfs    = len(yfs[0])

    for ihrd in xrange(nfs): ## type of hardening function
        if not(ivpsc_hard):
            hfnDill = gen_tempfile(prefix='hfs',ext='dll')
            with open(hfnDill, 'wb') as fo:
                dill.dump(hfs_labels[ihrd],fo)
                dill.dump(hfs_params[ihrd],fo)
                ## dump hardening parameters.

        for ieps in xrange(neps_eq):
            for iyld in xrange(nyfs): ## type of yield function

                ## Save using dill
                yfnDill = gen_tempfile(prefix='yfs-%s'%yfs_labels[iyld],ext='pck')
                hashcode = etc.gen_hash_code2(nchar=6)
                with open(yfnDill,'wb') as fo_:
                    if fo_.closed:
                        raise IOError, 'File was closed %s'%yfnDill
                    try:
                        print 'yfs[ieps][iyld]:',yfs[ieps][iyld]
                    except:
                        print len(yfs)
                        raise IOError,'Error!!!'

                    # dill.dump(yfs[ieps][iyld],fo_)
                    # dill.dump(yfs_labels[iyld],fo_)
                    # dill.dump(yfs_params[ieps][iyld],fo_)

                    # pickle.dump(yfs[ieps][iyld],fo_)
                    pickle.dump(yfs_labels[iyld],fo_)
                    pickle.dump(yfs_params[ieps][iyld],fo_)
                    pass

                cmd = 'python mk_run.py --f0 %f --r0 %f --r1 %f --nr %i'
                cmd = cmd + ' --hash %s --fnyld %s'
                cmd = cmd%(args.f0, args.r0, args.r1, args.nr, hashcode,yfnDill)

                if ivpsc_hard:
                    cmd = cmd+' --fnhrd_vpsc %s'%args.fnpickle_vpsc_hard
                else:
                    cmd = cmd+' --fnhrd %s'%hfnDill

                print '\n'*2,'-'*20
                print 'cmd:'
                print cmd
                print '-'*20,'\n'

                if args.dry: pass
                else:
                    os.system(cmd)
                    fn_tar = gen_tempfile(prefix='mk_%.4f_%s_%s_%.3f_'%(
                            eps_eq[ieps],yfs_labels[iyld],hfs_labels[ihrd],
                            args.f0),ext='tar')
                    subprocess.check_call(['tar','-cvf',fn_tar,
                                           'allFLD-%s.txt'%hashcode,
                                           'minFLD-%s.txt'%hashcode])
                    #'mk_fld_pp_%s.pdf'%hashcode

                    shutil.move(fn_tar, os.getcwd())
                    fnCollect.append(os.path.split(fn_tar)[-1])

    ## exit if dry run
    if args.dry: os._exit(0)
    else:        pass

    ## archive the results (and the pickled filed passed to the simulation as well)
    # - prepare..
    shutil.copy(args.fnpickle, os.getcwd())

    if ivpsc_hard:
        prefix_on_tarfile='MKVPSCHARD_%.3f_%s'%(args.f0,time.strftime('%Y%m%d-%H%M'))
    else:
        prefix_on_tarfile='MKVARHARD_%.3f_%s'%(args.f0,time.strftime('%Y%m%d-%H%M'))

    fn_result_tarfile = gen_tempfile(prefix=prefix_on_tarfile,ext='tar')
    fn_result_tarfile = os.path.split(fn_result_tarfile)[-1]

    # - archive
    cmd = ['tar','-cvf',fn_result_tarfile]
    for i in xrange(len(fnCollect)):
        cmd.append(fnCollect[i])
    subprocess.check_call(cmd)

    # - remove archived files
    for i in xrange(len(fnCollect)):
        os.remove(fnCollect[i])
    os.remove(os.path.split(args.fnpickle)[-1])

    ## move fn_result_tarfile to archive folder
    path_parent = os.path.split(os.getcwd())[0]
    path_archive = os.path.join(path_parent,'archive')
    print 'fn_result_tarfile:', fn_result_tarfile
    print 'is moved to:', path_archive
    shutil.move(fn_result_tarfile, path_archive)
