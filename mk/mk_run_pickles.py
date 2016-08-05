"""
Run mk_run using serialized Python functions
that represents the yield function and strain hardening function

The serialization of the Python function was enabled by the use of
dill package, available in "https://github.com/uqfoundation/dill".
"""


## example:
"""
$ python mk_run_pickles.py --f0 0.995 --r0 0 --r1 1 --nr 4 --fnpickle /home/younguj/repo/vpsc-fld/ipynb/FLD/IFsteel_EXP/yf_collection.dll
"""
if __name__=='__main__':
    import numpy as np
    # from mk import main as mk_main
    from mk.library.lib import gen_tempfile,rhos2ths
    import os,multiprocessing,time,argparse,dill,shutil
    from MP import progress_bar
    from mk_paths import findCorrectPsi
    import subprocess
    from MP.lib import etc

    parser = argparse.ArgumentParser()
    parser.add_argument('--f0',type=float,help='f0 value')
    parser.add_argument('--r0',type=float,help='rho start')
    parser.add_argument('--r1',type=float,help='rho end')
    parser.add_argument('--nr',type=int,help='number of rhos')
    parser.add_argument('--fnpickle', type=str,default=None,help='Pickle file name')
    args        = parser.parse_args()

    with open(args.fnpickle, 'rb') as fo:
        eps_eq     = dill.load(fo)
        yfs        = dill.load(fo)
        yfs_labels = dill.load(fo)
        hfs        = dill.load(fo)
        hfs_labels = dill.load(fo)
        results    = dill.load(fo) ## not important yet...

    fnCollect=[]
    for ihrd in xrange(len(hfs)): ## type of hardening function
        hfs[ihrd]
        hfnDill = gen_tempfile(prefix='yfs',ext='dll')
        with open(hfnDill, 'w') as fo:
            dill.dump(hfs[ihrd],fo)
        for ieps in xrange(len(eps_eq)):
            for iyld in xrange(len(yfs[ieps])): ## type of yield function
                yfnDill = gen_tempfile(prefix='hfs',ext='dll')
                hashcode = etc.gen_hash_code2(nchar=6)
                with open(yfnDill,'w') as fo:
                    dill.dump(yfs[ieps][iyld],fo)
                cmd = 'python mk_run.py --f0 %f --r0 %f --r1 %f --nr %i'
                cmd = cmd + ' --hash %s --mat 0 --fnyld %s --fnhrd %s'
                cmd = cmd%(args.f0, args.r0, args.r1, args.nr, hashcode,
                           yfnDill, hfnDill)
                print '\n'*2,'-'*20
                print 'cmd:'
                print cmd
                print '-'*20,'\n'*2
                os.system(cmd)
                fn_tar = gen_tempfile(prefix='mk_%.4f_%s_%s_%.3f_'%(
                        eps_eq[ieps],yfs_labels[iyld],hfs_labels[ihrd],
                        args.f0),ext='tar')
                subprocess.check_call(['tar','-cvf',fn_tar,
                                       'allFLD-%s.txt'%hashcode,
                                       'minFLD-%s.txt'%hashcode])
                shutil.move(fn_tar, os.getcwd())
                fnCollect.append(os.path.split(fn_tar)[-1])

    ## archive the results (and the pickled filed passed to the simulation as well)
    # - prepare..
    shutil.copy(args.fnpickle, os.getcwd())
    fn_tar = gen_tempfile(
        prefix='mk_%.3f_%s'%(
            args.f0,time.strftime('%Y%m%d-%H%M')),ext='tar')
    fn_tar=os.path.split(fn_tar)[-1]

    # - archive
    cmd = ['tar','-cvf',fn_tar]
    for i in xrange(len(fnCollect)):
        cmd.append(fnCollect[i])
    subprocess.check_call(cmd)

    # - remove archived files
    for i in xrange(len(fnCollect)):
        os.remove(fnCollect[i])
    os.remove(os.path.split(args.fnpickle)[-1])

    ## move fn_tar to archive folder
    path_parent = os.path.split(os.getcwd())[0]
    path_archive = os.path.join(path_parent,'archive')
    print 'fn_tar:', fn_tar
    print 'is moved to:', path_archive
    shutil.move(fn_tar, path_archive)
