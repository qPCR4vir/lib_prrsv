#!/usr/local/bin/python
import sys
import os
from sequence import fasta
from sequence import bowtie as B
import subprocess
from utils import multiFile, mystemp,safeOFW
import ncbi

hg18db = '/r1/home/fishnet/sequences/Homo_sapiens/hg18/hg18'

def main(args,maxMismatch=0,seedLen=36):
    for f in args:
        bowtieFile,bowtiePath = mystemp(suffix='.bowtie')
        fStub,fExt = os.path.splitext(f)
        oFile = safeOFW('.'.join((fStub,'Hg18rm',fExt[1:])),'w')
        print os.curdir
        try:
            subprocess.call(['bowtie', '-f','-l',str(seedLen),'-n',str(maxMismatch),
                             hg18db, f, bowtiePath])
            humanQrys =B.readTitleSet(bowtiePath)
            if len(humanQrys)==0:
                raise RuntimeError, "No Alignments found"

            for rec in fasta.generalIterator(f):

                if rec.title in humanQrys:
                    continue
                else:
                    print >> oFile, rec
            oFile.close()

        except BaseException:
            os.unlink(oFile.name)
            print '%s -failed' %f
            raise

        finally:
            os.unlink(bowtiePath)

    return 0

if __name__ == '__main__':
    sys.exit( main(sys.argv[1:]))
