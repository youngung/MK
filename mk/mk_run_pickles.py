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

    parser = argparse.ArgumentParser()
    parser.add_argument('--f0',type=float,help='f0 value')
    parser.add_argument('--r0',type=float,help='rho start')
    parser.add_argument('--r1',type=float,help='rho end')
    parser.add_argument('--nr',type=int,help='number of rhos')
    parser.add_argument('--fnpickle', type=str,default=None,help='Pickle file name')
    args        = parser.parse_args()

    with open(args.fnpickle, 'rb') as fo:
        yfs        = dill.load(fo)
        yfs_labels = dill.load(fo)
        hfs        = dill.load(fo)
        hfs_labels = dill.load(fo)
        results    = dill.load(fo)

    fnCollect=[]
    for ihrd in xrange(len(hfs)):
        hfs[ihrd]
        hfnDill = gen_tempfile(prefix='yfs',ext='dll')
        with open(hfnDill, 'w') as fo:
            dill.dump(hfs[ihrd],fo)

        for iyld in xrange(len(yfs)):
            yfs[iyld]
            yfnDill = gen_tempfile(prefix='hfs',ext='dll')
            with open(yfnDill,'w') as fo:
                dill.dump(yfs[iyld],fo)
            cmd = 'python mk_run.py --f0 %f --r0 %f --r1 %f'
            cmd = cmd + ' --mat 0 --nr %i --fnyld %s --fnhrd %s'
            cmd = cmd%(args.f0, args.r0, args.r1, args.nr,
                       yfnDill, hfnDill)
            print '\n'*2,'-'*20
            print 'cmd:'
            print cmd
            print '-'*20,'\n'*2
            os.system(cmd)
            fn_tar = gen_tempfile(prefix='mk_%s_%s_%.3f_'%(
                    yfs_labels[iyld],hfs_labels[ihrd],
                    args.f0),ext='tar')
            subprocess.check_call(['tar','-cvf',fn_tar,'allFLD.txt','minFLD.txt'])
            shutil.move(fn_tar, os.getcwd())
            fnCollect.append(os.path.split(fn_tar)[-1])

    fn_tar = gen_tempfile(
        prefix='mk_%.3f_%s'%(
            args.f0,time.strftime('%Y%m%d-%H%M')),ext='tar')
    fn_tar=os.path.split(fn_tar)[-1]
    cmd = 'tar -cvf %s '%fn_tar
    for i in xrange(len(fnCollect)):
        cmd = '%s %s'%(cmd,fnCollect[i])
    os.system(cmd)
    for i in xrange(len(fnCollect)):
        os.remove(fnCollect[i])
