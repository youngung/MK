## running multithreaded mk runs
if __name__=='__main__':
    from lib import gen_tempfile
    import mk_paths
    import multiprocessing

    ncpu  = multiprocessing.cpu_count()
    pool  = multiprocessing.Pool(processes=ncpu)
    paths = mk_paths.returnPaths()

    f0 = 0.996

    print '---'
    results = []
    for i in xrange(len(paths)):
        r=pool.apply_async(func=main, args=(f0,paths[i]))
        results.append(r)

    print '----'
    raw_input('ready?')

    t0 = time.time()
    pool.close()
    pool.join()
    pool.terminate()
    wallClockTime = time.time()-t0

    ## end of calculation
    logFileNames =[]
    tTimes       =[]
    rstFileName = gen_tempfile(prefix='MK',affix='results')
    with open(rstFileName,'r') as rstFile:
        for i in xrange(len(results)):
            pathFuncName = paths[i].__name__
            logFN, tTime=results[i].get()
            rstFile.write('%i %6s %s %i\n'%(
                i,pathFUncName,logFN,tTime))
            logFileNames.append(logFN)
            tTimes.append(tTime)

    print 'All results are saved to %s'%rstFileName
