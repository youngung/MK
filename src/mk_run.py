"""
"""

## running multithreaded mk runs

def makeCommands(f0,psi0,th,logFileName):
    """
    Arguments
    ---------
    f0
    psi0
    th
    logFileName
    """
    from lib      import gen_tempfile
    stdoutFileName = gen_tempfile(
        prefix='stdout-mkrun')
    cmd = 'python mk.py --fn %s -f %5.5f -p %+5.5f -t %+5.5f > %s'%(
        logFileName,f0,psi0,th,stdoutFileName)
    print 'cmd:',cmd
    return cmd

def test_run():
    run(f0=0.995,psi0=0.,th=0.,logFileName='dum')

def run(*args):
    import os,subprocess, time
    t0=time.time()
    cmd = makeCommands(*args)
    dt  = time.time()-t0
    os.system(cmd)
    return dt

if __name__=='__main__':
    import numpy as np
    from mk import main as mk_main
    from lib import gen_tempfile, rhos2ths
    import os,multiprocessing,time
    from MP import progress_bar
    from mk_paths import findCorrectPsi

    uet=progress_bar.update_elapsed_time

    f0 = 0.996
    print '---'

    ## rho to theta? ( should be later passed as arguments)
    rhos = [-0.5, -0.25, 0, 0.5, 1]#, 1.5, 2, 2.25, 2.5]
    ths  = rhos2ths(rhos)
    psi0 = 0. ## should be determined in mk_paths

    logFileNames=[]
    for i in xrange(len(ths)):
        _psi0s_=findCorrectPsi(ths[i])
        logFileNames.append([])
        for j in xrange(len(_psi0s_)):
            psi0 = _psi0s_[j]
            logFileName = gen_tempfile(
                prefix='mk-f0%3.3i-psi%+6.3f-th%+6.3f'%(
                    int(f0*1e3),psi0,ths[i]),
                affix='log',i=(i*10000+j))
            logFileNames[i].append(logFileName)

    ncpu  = multiprocessing.cpu_count()
    pool  = multiprocessing.Pool(processes=ncpu)

    results = []
    for i in xrange(len(ths)):
        _psi0s_=findCorrectPsi(ths[i])
        results.append([])
        for j in xrange(len(_psi0s_)):
            psi0 = _psi0s_[j]
            # func = mk_main
            func = run
            r=pool.apply_async(
                func=func,
                args=(f0,psi0,ths[i],logFileNames[i][j]))
            results[i].append(r)
    print '----'

    t0 = time.time()
    pool.close(); pool.join(); pool.terminate()
    wallClockTime = time.time()-t0

    ## end of calculation
    tTimes       =[]
    rstFileName = gen_tempfile(prefix='MK',affix='results')
    with open(rstFileName,'w') as rstFile:
        for i in xrange(len(results)): ## paths
            for j in xrange(len(results[i])): ## psi0s
                tTime = results[i][j].get()
                tTimes.append(tTime)
                logFN = logFileNames[i][j]
                rstFile.write('%i %s\n'%(j,logFN))
                logFileNames.append(logFN)
            rstFile.write('--\n')

    cpuRunTime=np.array(tTimes).sum()

    for i in xrange(len(logFileNames)):
        print '%4.4i   %s'%(i,logFileNames[i])

    uet(wallClockTime,'Total wallclocktime:');print

    # uet(cpuRunTime,   'CPU running time   :');print
    # speedup = cpuRunTime / wallClockTime
    # print 'speedup: %5.2f'%(speedup)
    print 'All results are saved to %s'%rstFileName
    print 'Fin --'

def read(fn):
    f    = float(fn.split('mk-f')[1].split('-psi')[0])*1e-3
    psi0 = float(fn.split('psi')[1].split('-th')[0])
    th   = float(fn.split('-th')[1][:6])
    with open(fn) as FO:
        data_line = FO.readlines()[1]
        elements = data_line.split()
    return map(float, elements),f,psi0,th, data_line

def pp(masterFileName):
    """
    Argument
    --------
    masterFileName
    """
    import numpy as np
    numFail=0
    fileFail=[]

    fileFLDall = open('allFLD.txt','w')
    fileFLDmin = open('minFLD.txt','w')

    with open(masterFileName) as FO:
        blocks = FO.read().split('--\n')[:-1:]
        dat_min_master=[]
        print 'number of blocks',len(blocks)
        for i in xrange(len(blocks)): ## each block
            eachBlock = blocks[i]
            linesInBlock = eachBlock.split('\n')[1:-1:]
            print linesInBlock

            ## find the minimum |(E1,E2)|
            min_rad =1.5
            dat_min = None
            data_min_line = None
            # print linesInBlock

            for j in xrange(len(linesInBlock)):
                line = linesInBlock[j]
                ind, fn = line.split()
                data, f, psi0, th, data_line = read(fn)
                fileFLDall.write('%s'%data_line)
                epsRD, epsTD, psi0, psif, \
                    sigRD,sigTD,sigA,T,dt = data[:9]
                rad = np.sqrt(epsRD**2+epsTD**2)
                if rad<min_rad:
                    data_min = data[::]
                    rad = min_rad
                    data_min_line = data_line

            dat_min_master.append(dat_min)
            fileFLDmin.write('%s'%data_min_line)

    fileFLDall.close(); fileFLDmin.close()

    print 'numFail:',numFail
    print 'FileFail:'
    for i in xrange(numFail):
        print fileFail[i]

def test_pp(fn='/local_scratch/MK-6e59e6-results.txt'):
    pp(fn)
