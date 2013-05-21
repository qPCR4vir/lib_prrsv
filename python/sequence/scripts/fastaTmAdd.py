#!/usr/local/bin/python
import sys
import os
from sequence import fasta

def main(args):
    for rec in fasta.generalIterator(args[0]):
        rec.title+='|%.1f' % rec.Tm()
        print rec
            
    return 0
    

if __name__ == '__main__':
    sys.exit( main(sys.argv[1:]))
