## run mk_run from a pickled file that contains
## various yield functions and strain-hardening functions


if __name__=='__main__':
    import numpy as np
    # from mk import main as mk_main
    from mk.library.lib import gen_tempfile,rhos2ths
    import os,multiprocessing,time,argparse,dill
    from MP import progress_bar
    from mk_paths import findCorrectPsi

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--f0',type=float,help='f0 value')
    parser.add_argument(
        '--r0',type=float,help='rho start')
    parser.add_argument(
        '--r1',type=float,help='rho end')
    parser.add_argument(
        '--nr',type=int,help='number of rhos')
    parser.add_argument(
        '--fnpickle', type=str,default=None,
        help='Pickle file name')
    args        = parser.parse_args()

    with open(args.fnpickle, 'rb') as fo:
        yfs        = dill.load(fo)
        yfs_labels = dill.load(fo)
        hfs        = dill.load(fo)
        hfs_labels = dill.load(fo)
        results    = dill.load(fo)

    # for ihrd in xrange(len(hfs)):
    for ihrd in xrange(1):
        hfs[ihrd]
        hfnDill = gen_tempfile(prefix='yfs',ext='dll')
        with open(hfnDill, 'w') as fo:
            dill.dump(hfs[ihrd],fo)

        # for iyld in xrange(len(yfs)):
        for iyld in xrange(1):
            yfs[iyld]
            yfnDill = gen_tempfile(prefix='hfs',ext='dll')
            with open(yfnDill,'w') as fo:
                dill.dump(yfs[iyld],fo)

            cmd ='python mk_run.py --f0 %f --r0 %f --r1 %f --mat 0 --nr %i --fnyld %s --fnhrd %s'%(
                args.f0, args.r0, args.r1, args.nr, yfnDill, hfnDill)

            print '\n'*2,'-'*20
            print 'cmd:'
            print cmd
            print '-'*20,'\n'*2
            os.system(cmd)
