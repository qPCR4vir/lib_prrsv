#!/usr/local/bin/python
import sys
import utils
from glob import glob

from sequence import blastNoSQL as B
from ncbi import taxonomy

targetTaxID = 54290
viralAlignmnets=glob("")
crossHybeAlignments = (,)
noCrossHybe = ""
minAlignLen=14

B.uniqueAlignments(viralAlignmnets,crossHybeAlignments,file(noCrossHybe,'w'))
gis=utils.noneDict(giInfo.GiInfo._table.Gi.valueList(
    Family_Tax_ID=targetFamilyTaxID) )                
#print len(famGis)

subjects = {}

for l in file(noCrossHybe):
    try:
        gi = int(B.giFromM8name(l.split()[1]))
    except IndexError:
        continue

    #print 'beep'
    
    if gi not in gis:
        continue

    #print '\tbeep'
    if int(l.split()[3]) >= minAlignLen:
        subj = l.split()[1]
        if subj not in subjects:
            subjects[subj] = []
        subjects[subj].append(l)
        
        
for aSet in subjects.values():
    if len(aSet) < 2:
        continue
    aSet.sort(key=lambda x: min([x[8],[9]]))
    s = 0
    setReads=[]
    for a in aSet:
        sS,sE = [int(x) for x in a[8:10]]
        e=max((sS,sE))
        qS,qE = a[6:8]
        if qE < qS:
            qTmp = qE
            qE=qS
            qS=qTmp
        setReads.append(reads[a[0]].subrecord(qS,qE))
        if s==0:
            print '\t'.join([str(x) for x in a])
        else:
            print '\t'.join([str(x) for x in a])+'\t%s' %(e-s)
        s=min((sS,sE))
    for r in setReads:
        print r
    print
              
