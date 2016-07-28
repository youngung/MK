"""
"""
import matplotlib as mpl
mpl.use('Agg') ## In case X-window is not available.

## running multithreaded mk runs
def read(fn):
    f    = float(fn.split('mk-f')[1].split('-psi')[0])*1e-3
    psi0 = float(fn.split('psi')[1].split('-th')[0])
    th   = float(fn.split('-th')[1][:6])
    with open(fn) as FO:
        data_line = FO.readlines()[1]
        elements = data_line.split()
        matA_FN=elements[9]
        matB_FN=elements[10]
        ss_FN  =elements[11]
    return map(float, elements[:9]),f,psi0,th,\
        data_line,matA_FN,matB_FN,ss_FN

def postAnalysis(masterFileName):
    """
    A master is generated by mk_run that lists
    name of files generated by each calculation

    Argument
    --------
    masterFileName
    """
    import numpy as np
    import mk.library.parser
    numFail=0
    fileFail=[]

    fileFLDall = open('allFLD.txt','w')
    fileFLDmin = open('minFLD.txt','w')

    with open(masterFileName) as FO:
        blocks = FO.read().split('--\n')[:-1:]
        dat_min_master=[]
        # print 'number of blocks',len(blocks)
        for i in xrange(len(blocks)): ## each block
            eachBlock = blocks[i]
            linesInBlock = eachBlock.split('\n')[0:-1:]
            # print linesInBlock

            ## find the minimum |(E1,E2)|
            min_rad       = 2.0
            dat_min       = None
            ind_min       = None
            data_min_line = None
            for j in xrange(len(linesInBlock)):
                line    = linesInBlock[j]
                ind, fn = line.split()
                try:
                    data, f, psi0, th, data_line,\
                        matA_FN, matB_FN, ss_FN = read(fn)
                except:
                    pass
                else:
                    fileFLDall.write('%s'%data_line)
                    epsRD, epsTD, psi0, psif, \
                        sigRD,sigTD,sigA,T,dt = data[:9]

                    if np.isnan(epsRD) or np.isnan(epsTD):
                        fileFail.append(fn)
                        numFail=numFail+1
                    else:
                        rad = np.sqrt(epsRD**2+epsTD**2)
                        if rad<min_rad:
                            dat_min = data[::]
                            min_rad = rad
                            ind_min = j
                            data_min_line = data_line

            dat_min_master.append(
                [dat_min,matA_FN,matB_FN,ss_FN])
            if type(data_min_line).__name__!='NoneType':
                fileFLDmin.write('%s'%data_min_line)

    fileFLDall.close(); fileFLDmin.close()

    ## iplot?
    import matplotlib.pyplot as plt
    from mk.library.lib import draw_guide
    fig = plt.figure(figsize=(7,6))
    ax1=fig.add_subplot(221);ax2=fig.add_subplot(222)
    ax3=fig.add_subplot(223);ax4=fig.add_subplot(224)
    dat=np.loadtxt(fileFLDmin.name,dtype='str').T
    dat=dat[:9]

    ax1.plot(dat[1],dat[0],'o')
    dat=np.loadtxt(fileFLDall.name,dtype='str').T
    dat=dat[:9]
    ax2.plot(dat[1],dat[0],'o')
    draw_guide(ax1,r_line=[-0.5,0,1,2,2.5],max_r=2)
    draw_guide(ax2,r_line=[-0.5,0,1,2,2.5],max_r=2)
    ax1.set_aspect('equal');ax2.set_aspect('equal')

    ##
    for i in xrange(len(dat_min_master)):
        dat_min, matA_FN, matB_FN, ss_FN = dat_min_master[i]
        mk.library.parser.plotMat(matA_FN,ax=ax3,
                       color='red',linestyle='-')
        mk.library.parser.plotMat(matB_FN,ax=ax3,
                       color='blue',linestyle='--')
        mk.library.parser.plotEtc(ss_FN,ax=ax4)
    fig.savefig('mk_fld_pp.pdf')

def test_pp(fn='/local_scratch/MK-6e59e6-results.txt'):
    postAnalysis(fn)

def makeCommands(f0,psi0,th,logFileName,mat,fnyld,fnhrd):
    """
    Arguments
    ---------
    f0
    psi0
    th
    logFileName
    mat
    fnyld
    fnhrd
    """
    from mk.library.lib      import gen_tempfile
    stdoutFileName = gen_tempfile(
        prefix='stdout-mkrun')
    # stdoutFileName ='/tmp/dump'
    cmd = 'python main.py --fn %s -f %5.4f -p %+6.1f -t %+7.2f --mat %i --fnyld %s --fnhrd %s > %s'%(
        logFileName,f0,psi0,th,mat,fnyld,fnhrd,stdoutFileName)
    print 'cmd:',cmd
    return cmd

def test_run():
    run(f0=0.995,psi0=0.,th=0.,logFileName='dum')

def run(*args):
    import os
    cmd = makeCommands(*args)
    os.system(cmd)
    return cmd

def prepRun(*args):
    import os
    cmd = makeCommands(*args)
    return cmd

"""
$ python mk_run.py --f0 0.995 --r0 -0.5 --r1 1
             --nr 4 --mat 0 --fnyld <> --fnhrd <>
"""
if __name__=='__main__':
    import numpy as np
    # from mk import main as mk_main
    from mk.library.lib import gen_tempfile, rhos2ths
    import os,multiprocessing,time
    from MP import progress_bar
    from mk_paths import findCorrectPsi
    import argparse

    ## Arguments parsing
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
        '--mat', type=int, default=0,
        help='Material card in materials.py - see <def library> in materials.py\n(0: IsoMat, 1:)')
    parser.add_argument(
        '--fnyld', type=str,default=None,
        help='yield function pickled file name')
    parser.add_argument(
        '--fnhrd', type=str,default=None,
        help='strain hardening function pickled file name')

    ## rho to theta? ( should be later passed as arguments)
    args        = parser.parse_args()
    rhos = np.linspace(args.r0,args.r1,args.nr)

    print 'rhos:', rhos

    ths  = rhos2ths(rhos)
    uet  = progress_bar.update_elapsed_time
    logFileNames=[]
    k=0
    print '%3s %6s %5s %5s %60s'%('k','rho','th','psi0','logFileName')
    p0s=[]
    print '---'
    for i in xrange(len(ths)): ## each rho
        _psi0s_=findCorrectPsi(ths[i]*180/np.pi)
        p0s.append(_psi0s_)
        logFileNames.append([])
        for j in xrange(len(_psi0s_)): ## each psi
            psi0 = _psi0s_[j]
            logFileName = gen_tempfile(
                prefix='mk-f0%3.3i-psi%+6.3f-th%+6.3f'%(
                    int(args.f0*1e3),psi0,ths[i]),
                affix='log',i=(i*10000+j))
            logFileNames[i].append(logFileName)
            k=k+1
            print '%3i %6.2f %5.1f %5.1f %60s'%(
                k, rhos[i], ths[i]*180/np.pi,
                psi0*180/np.pi,logFileName)
            pass
        print '-'*83
        pass

    ncpu  = multiprocessing.cpu_count()
    pool  = multiprocessing.Pool(processes=ncpu)

    results = []
    for i in xrange(len(ths)):
        results.append([])
        for j in xrange(len(p0s[i])):
            psi0 = p0s[i][j]
            r=pool.apply_async(
                func=run,
                args=(args.f0,
                      psi0*180/np.pi,
                      ths[i]*180/np.pi,
                      logFileNames[i][j],
                      args.mat,
                      args.fnyld,args.fnhrd
                ))
            results[i].append(r)
            pass
        pass

    t0 = time.time()
    pool.close(); pool.join(); pool.terminate()
    wallClockTime = time.time()-t0

    ## end of calculation
    rstFileName = gen_tempfile(prefix='MK',affix='results')
    with open(rstFileName,'w') as rstFile:
        for i in xrange(len(results)): ## paths
            for j in xrange(len(results[i])): ## psi0s
                cmd = results[i][j].get()
                logFN = logFileNames[i][j]
                rstFile.write('%i %s\n'%(j,logFN))
                logFileNames.append(logFN)
            rstFile.write('--\n')
            pass
        pass

    for i in xrange(len(logFileNames)):
        print '%4.4i   %s'%(i,logFileNames[i])

    uet(wallClockTime,'Total wallclocktime:');print
    print 'All results are saved to %s'%rstFileName
    postAnalysis(rstFileName)
    print 'Fin --'
