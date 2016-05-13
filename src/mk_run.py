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
    print cmd
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

    ncpu  = multiprocessing.cpu_count()
    pool  = multiprocessing.Pool(processes=ncpu)

    f0 = 0.996
    print '---'
    results = []

    ## rho to theta? ( should be later passed as arguments)
    rhos = [-0.5, -0.25, 0, 0.5, 1]#, 1.5, 2, 2.25, 2.5]
    ths  = rhos2ths(rhos)
    psi0 = 0. ## should be determined in mk_paths

    logFileNames=[]

    for i in xrange(len(ths)):
        _psi0s_=findCorrectPsi(ths[i])
        for j in xrange(len(_psi0s_)):
            psi0 = _psi0s_[j]
            logFileName = gen_tempfile(
                prefix='mk-f0%3.3i-psi%+6.3f-th%+6.3f'%(
                    int(f0*1e3),psi0,ths[i]),
                affix='log',i=(i*10000+j))
            logFileNames.append(logFileName)

    for i in xrange(len(ths)):
        _psi0s_=findCorrectPsi(ths[i])
        for j in xrange(len(_psi0s_)):
            psi0 = _psi0s_[j]
            # func = mk_main
            func = run
            r=pool.apply_async(
                func=func,
                args=(f0,psi0,ths[i],logFileNames[i]))
            results.append(r)
    print '----'

    t0 = time.time()
    pool.close(); pool.join(); pool.terminate()
    wallClockTime = time.time()-t0

    ## end of calculation
    tTimes       =[]
    rstFileName = gen_tempfile(prefix='MK',affix='results')
    with open(rstFileName,'w') as rstFile:
        for i in xrange(len(results)):
            tTime = results[i].get()
            tTimes.append(tTime)
            logFN = logFileNames[i]
            rstFile.write('%i %s\n'%(
                    i,logFN))
            logFileNames.append(logFN)

    cpuRunTime=np.array(tTimes).sum()

    for i in xrange(len(logFileNames)):
        print '%4.4i   %s'%(i,logFileNames[i])

    uet(wallClockTime,'Total wallclocktime:');print

    # uet(cpuRunTime,   'CPU running time   :');print
    # speedup = cpuRunTime / wallClockTime
    # print 'speedup: %5.2f'%(speedup)
    print 'All results are saved to %s'%rstFileName
    print 'Fin --'



def pp(fn):
    pass
