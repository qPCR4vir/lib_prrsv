#!/usr/local/bin/python -u

# 	$Id: makeIdentityScreener.py,v 1.3 2011/05/04 16:32:19 kael Exp $	
__version__ = tuple([int(x) for x in
                     '$Revision: 1.3 $'.split()[1].split('.')])
__author__ = "Kael Fischer"


import os,sys,os.path
import stat
import utils
import optparse


# defaults
origFasta = None
faNameFile = "IdentityScreenerSliceNames.txt"
outFastaName = "screened.fasta"
blastTimeout = 600
sliceSize=None
sliceCount=None

TEMPLATE="""#!/usr/local/bin/python -u
import os,sys,os.path
import shutil
import commands
import time
from utils.timeprofile import TimeProfiler
from sequence import fasta,screen
from types import *
from ncbi.giInfo import GiInfo
import grid
from grid import rpyc

#sys.exit(1)

originalFasta = "%(origFasta)s"
faNameFile = "%(faNameFile)s"
outFastaName="%(outFastaName)s"

def biggestRemoteRecord(servers):
    "return the largest remote record"

    nMax = 0
    winner = None
    for i,s in enumerate(servers):
        n = s.execute("len(screener.longestRemaining)","eval")
        if n > nMax:
            nMax = n
            winner = i
    return rpyc.obtain(servers[winner].namespace.screener.longestRemaining)


def makeCheckpointDir(dirSuffix):
    #
    # Make a checkpoint directory
    #
    dirName = os.path.join(os.getcwd(),"checkpoint_cycle%%s"%%dirSuffix)
    print "making checkpoint directory"
    os.system("mkdir %%s" %% (dirName))
    outFile.flush()
    shutil.copy(faNameFile,dirName)
    shutil.copy(outFile.name,dirName)
    shutil.copy(__file__,dirName)
    requests = []
    for i,asyncOs in enumerate(osSys):
        requests.append(asyncOs("cp %%s %%s" %%
                    (servers[i].namespace.screener.localPath,
                     os.path.join(dirName,
                     os.path.split(servers[i].namespace.screener.inPath)[1]))))
        
    while False in [r.is_ready for r in requests]:
        time.sleep(5)

    return dirName


def printLongInfo(longest):
    try:
        # if the title is just a Gi print it nicely
        longGi = int(longest.title.strip())
        gii = GiInfo(longGi)
        print str(gii)
        print gii.fgsnText()
    except:
        print longest.title
    print len(longest),longest.DRNABaseCount(),max(longest.DRNABaseLengths())

if __name__ == '__main__':

    RESUME = os.access(faNameFile,os.F_OK)
    RUN = True
    BLAST_TIMEOUT = %(blastTimeout)s

    timer = TimeProfiler()
    timer.mark('main')

    seqSliceSize = %(sliceSize)s
    partCt = %(sliceCount)s

    workingDir=os.getcwd()
    chkPointName = None

    outFile = file(os.path.join(workingDir,outFastaName),'a')

    qsubOptions=''

    # complete restart
    timer.mark('fastaSplit')
    if not RESUME:
        startName = originalFasta
        inFileName = startName

        if partCt == None:
            seqSize = os.stat(originalFasta).st_size
            partCt = seqSize/seqSliceSize +1
            
        seqSliceNames = fasta.splitFasta(originalFasta,partCt,tmpDir=workingDir)
        tmpf=file(faNameFile,'w')
        tmpf.write('\\n'.join(seqSliceNames))
        tmpf.close()
        
    # reuse split files    
    if RESUME:
        seqSliceNames = []
        tmpf = file(faNameFile)
        for l in tmpf:
            try:
                print 
                seqSliceNames.append(
                    os.path.join(os.getcwd(),
                                 os.path.split(l.rstrip())[1]))
            except:
                pass
        tmpf.close()

    fSplitTime = timer.lastdiff()

    print 'Input Slice Names'
    print seqSliceNames, len(seqSliceNames)

    if RUN:
        lastCheckpointName="_Done"
        timer.mark('rpycStart')
        serverCt=len(seqSliceNames)
        servers = rpyc.GridServers(serverCt,verbose=True)
        screeners = []
        droppers = []
        cleaners = []
        setups = []
        requests = []
        returnValues = []
        osSys = []

        #set up servers and make Async objects for
        #setting up, screening and cleaning up
        for s in servers:
            s.execute("from sequence import fasta")
            s.execute("screener=fasta.IterativeScreen()")
            setups.append(rpyc.Async(s.namespace.screener.setup))
            screeners.append(rpyc.Async(s.namespace.screener.screen))
            droppers.append(rpyc.Async(s.namespace.screener.dropShorties))
            cleaners.append(rpyc.Async(s.namespace.screener.cleanUp))
            osSys.append(rpyc.Async(s.modules.os.system))
        rpycStartTime = timer.lastdiff()
        try:
            timer.mark('copyStart')
            for i,s in enumerate(setups):
                requests.append(s(seqSliceNames[i]))
                print seqSliceNames[i]

            print "nodes copying chunks"
            while False in [r.is_ready for r in requests]:
                print [r.is_ready for r in requests].count(True), " nodes ready"
                time.sleep(0.5)
            copyTime = timer.lastdiff()
            print 'Chunk copy time: %%.1f' %% copyTime
            print "Longest Base Counts:",[r.result for r in requests]
            #
            # Main Loop coming up
            #
            cycle = 0
            longest = biggestRemoteRecord(servers)
            while max(longest.DRNABaseLengths())>=50:
                print 'beep'
                printLongInfo(longest)
                file("longest.fasta",'w').write(str(longest))
                os.system("formatdb -p F -i longest.fasta")

                outFile.write(str(longest))
                outFile.write("\\n")


                requests=[]
                for i,s in enumerate(screeners):
                    #submit screening job
                    requests.append(s(os.path.join(os.getcwd(),
                                                   'longest.fasta')))
                print "blasts running"

                timer.mark('blastSubmit')
                while False in [r.is_ready for r in requests]:
                    #print [r.is_ready for r in requests].count(True), " chunks done"
                    time.sleep(0.5)
                    if timer.lastdiff() > BLAST_TIMEOUT:
                        # server died?
                        # start from last check point
                        print "rpyc server timed out on megablast"
                        del requests
                        del setups
                        outFile.close()
                        requests = [cleaner() for cleaner in cleaners ]
                        time.sleep(60)
                        del cleaners
                        del droppers
                        del screeners
                        servers.kill()
                        del servers
                        print "closed output files, killed servers. some remote tmpFiles may be undeleted."

                        pid=subprocess.Popen("cd %%s ; nohup ./identityScreen_Rpyc.py > is.log 2>&1" %% (chkPointName),
                                             shell=True).pid
                        print "spawned new job (pid:%%s), in working directory:\\n%%s" (pid, chkPointName) 
                        print "goodbye <sniff>"
                        sys.exit(1)

                        
                print 'Remote BLAST time: %%0.1f' %% timer.elapsed('blastSubmit')        
                timer.mark('blastFinished')
                timer.observe('blastSubmit','blastFinished','blast')

                print [r.is_ready for r in requests].count(True), " chunks done"
                print "Longest Base Counts:",[r.result for r in requests]

                cycle +=1
                print "cycle ", cycle, " done"
                longest = biggestRemoteRecord(servers)
                if max(longest.DRNABaseLengths()) <= 50:
                    requests = [d() for d in droppers ]
                    while False in [r.is_ready for r in requests]:
                        time.sleep(1)
                    print "Shorties dropped"
                    print "Longest Base Counts:",[r.result for r in requests]
                    longest = biggestRemoteRecord(servers)
                    
                    
                if cycle %% 250 == 0 or cycle == 1:
                    chkPointName=makeCheckpointDir(str(cycle))

                    

                    
##         # comment finaly and uncomment here for mystery crash information
##         except:
##             requests = [cleaner() for cleaner in cleaners ]
##             while False in [r.is_ready for r in requests]:
##                 time.sleep(1)
##             raise


        finally:
            
            # make terminal checkpoint
            makeCheckpointDir(lastCheckpointName)
            outFile.close()
            requests = [cleaner() for cleaner in cleaners ]
            while False in [r.is_ready for r in requests]:
                time.sleep(1)

            print "remote files removed\\n"
            try:
                print "Total BLAST time: %%s" %% timer.total('blast')
                print "Average BLAST time: %%s" %% timer.averageElapsed('blast')
            except:
                print "Total BLAST time: None" %% timer.total('blast')
                print "Average BLAST time: None" %% timer.averageElapsed('blast')
                
            print "Total RunTime: %%s" %% timer.timegap()
            
"""

def main(sysargs=[]):
    """
    """
    oneLineUsage = "Usage: %prog [options] <input fasta file>"

    # set a long desccription for the help message
    description="""This program makes a multiprocess (RPC) screener to remove
redundant sequence from a fastafile.

The output file is named identityScreen_Rpyc.py.
"""

    op = optparse.OptionParser(
        oneLineUsage,description=description,
        version="%prog " + '.'.join([str(x) for x in __version__]))

    # OPTION DEFINITIONS
    op.add_option('-o','--outputfile',dest="outFileName",default='identityScreen_Rpyc.py',
                  help="name of output file. '-' means: STDOUT")
    op.add_option('-c','--clobber',action="store_true", default=False,
                  dest="clobber", help="allow output file to be overwritten.")
    op.add_option('-n','--slice-names',dest='faNameFile',
                  help='Name of file that will hold the slice names')
    op.add_option('-S','--slice-size',dest='sliceSize',type='int',
                  default=sliceSize,help='maximum size, in bytes, of each slice file')
    op.add_option('-C','--slice-count',dest='sliceCount',type='int',
                  default=sliceSize,help='number of slices')


    # add defaults to help messages     
    for o in op.option_list:            
        if o.type == None:              
            continue                    
        if o.help == None:              
            o.help = "Default: %default"
        else:                           
            o.help+= " Default: %default"

    # OPTION PARSING                    
    (opts,args) = op.parse_args(sysargs)

    try:
        if opts.sliceSize is None and opts.sliceCount is None:
            raise Exception, "sliceCount (-C) or sliceSize (-S) must be set."
        if len(args) != 1:
            raise Exception, "Input fasta file must me specified."
        
        outFile = utils.safeOFW(opts.outFileName,
                                append=False,clobber=opts.clobber) # open the output file
    except Exception, eData:
        #raise
        print >> sys.stderr, ("\nUsage Error: %s\n"%eData.message)
        print >> sys.stderr, op.format_help()
        return 1


    # Do It
    print >> outFile, TEMPLATE % { 
        'origFasta': args[0],
        'faNameFile':faNameFile,
        'outFastaName':outFastaName,
        'blastTimeout':blastTimeout,
        'sliceSize':opts.sliceSize,
        'sliceCount':opts.sliceCount
        }
    outFile.close()
    if opts.outFileName != '-':
        os.chmod(opts.outFileName,
                 os.stat(opts.outFileName)[stat.ST_MODE] | 0100)
        
        


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

 
