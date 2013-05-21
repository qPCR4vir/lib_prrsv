#!/usr/local/bin/python
import sys
import os
from sequence import fasta
from sequence import blastNoSQL as B
import subprocess
from utils import multiFile, mystemp,safeOFW
import ncbi

hg18db = '/r1/home/fishnet/sequences/Homo_sapiens/hg18/all.fa'

def main(args):
    for f in args:
        mbrFile,mbrPath = mystemp(suffix='.mbr')
        fStub,fExt = os.path.splitext(f)
        oFile = safeOFW('.'.join((fStub,'Hg18rm',fExt[1:])),'w')
        print os.curdir
        try:
            subprocess.call(['megablast', '-m8', '-i',f,
                             '-d', hg18db,'-o', mbrPath])
            humanQrys =B.m8TitleSet(mbrPath,qrySet=True)
            for rec in fasta.generalIterator(f):
                if len(humanQrys)==0:
                    break
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
            os.unlink(mbrPath)

    return 0

if __name__ == '__main__':
    sys.exit( main(sys.argv[1:]))
