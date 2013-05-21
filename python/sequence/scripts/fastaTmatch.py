#!/usr/local/bin/python
import sys
import os
from sequence import fasta
from sequence import blastNoSQL as B
import subprocess
from utils import multiFile, mystemp
import ncbi

def main(args):
    fFile= args[-1]
    titles=set(args[:-1])

    for rec in fasta.generalIterator(fFile):
        if rec.title in titles:
            print rec
            
    return 0
    

if __name__ == '__main__':
    sys.exit( main(sys.argv[1:]))
