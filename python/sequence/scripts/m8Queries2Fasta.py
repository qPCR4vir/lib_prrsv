#!/usr/local/bin/python
from sequence import blastNoSQL as B
import sys
from stat import *
import os
import os.path

try:
    minBitScore=float(sys.argv[1])
    queryPath = sys.argv[2]
    resultPath = sys.argv[3]
except:
    print "usage: m8Queries2Fasta.py <bit score> <query path> <result path>"
    sys.exit(1)


if os.path.isdir(queryPath):

    rStub = os.path.split(os.path.splitext(resultPath)[0])[1]

    queryTry = os.path.join(queryPath,rStub+'.fasta')
    
    print queryTry
    if os.access(queryTry,os.R_OK):
        queryPath=queryTry
    else:
        queryTry = os.path.join(queryPath,
                                os.path.splitext(rStub)[0]+'.fasta')
    print queryTry
    if os.access(queryTry,os.R_OK):
        queryPath=queryTry
       
        
outPath=os.path.splitext(resultPath)[0] + '.bs%s_qry.fasta' % int(round(minBitScore))

B.m8Queries2Fasta(resultPath,queryPath,outPath,minBitScore=minBitScore)



