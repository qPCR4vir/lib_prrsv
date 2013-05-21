#!/usr/local/bin/python

import sys

from utils import MultiDictCluster

iFile=sys.stdin
oFile=sys.stdout
mdc=MultiDictCluster()
iCol = iFile.next().strip().split('\t')[3:]
for g in iCol:
    mdc[g]={}
iFile.next()
for r in iFile:
    rF=r.strip().split()
    a=rF[0]
    vals=rF[3:]
    for g,v in zip(iCol,rF):
        mdc[g][a]=v
        
mdc.clusterFmt(oFile)
