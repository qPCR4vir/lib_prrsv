#!/usr/local/bin/python

from sequence import blastNoSQL
from utils import readArgsOrFiles

import sys
import os


def extractHSPs(inFile,ignoreTies=False):
    rv={}
    for m8rec in blastNoSQL.m8generator(inFile):
        if m8rec['query'] not in rv:
            rv[m8rec['query']]=[m8rec]
        elif  m8rec['score'] > rv[m8rec['query']][0]['score']:
            rv[m8rec['query']]=[m8rec]
        elif m8rec['score'] == rv[m8rec['query']][0]['score']:
            rv[m8rec['query']].append(m8rec)

    kList = rv.keys()
    kList.sort()
    for k in kList:
        for m8 in rv[k]:
            print m8['_str_']
            if ignoreTies:
                break
    


if __name__ == "__main__":
    if len (sys.argv) < 2 or sys.argv[1]=='-h':
        print "%s: [-T] <inFile>" % os.path.split(sys.argv[0])[0]
        print "\t-T ignore ties\n"
        sys.exit(1)
    ignoreTies = False 
    args = sys.argv[1:]
    if args[0]== '-T':
        ignoreTies = True
        args.pop(0)    

    extractHSPs(args[0],ignoreTies)

    

    
