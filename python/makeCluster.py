#!/usr/local/bin/python
import re
from glob import glob
from types import ListType
from utils import MultiDictCluster
from Axon import ATF
from prrsv import chipMap,gprFN2chipNumber

varianceFile=file('replicateVar.txt','w')


def mean532(recs):
    rv=0
    v=[]
    if type(recs) != ListType:
        rec=[res]
    for rec in recs:
        if rec.Flags>=0: # and re.match(r'\d+',rec.ID):
            v.append(rec['F532 Median - B532'])
            rv+=rec['F532 Median - B532']
    if len(v)>0:
        rv/=float(len(recs))
        print >> varianceFile, '\t'.join([str(x) for x in [rec.Name]+v])
        return rv
    else:
        return None
    


mdc=MultiDictCluster()

for fn in glob('Prrsv_*idHack.gpr'):
    gpr=ATF.GPR(fn)
    eName=fn.split('.')[0]
    eName=chipMap[gprFN2chipNumber(eName)]
    expt={}
    for i in set((r.ID for r in gpr.dataLines)):
        expt[i]=mean532(gpr[i])
    mdc.eAppend(expt,eName,sumNorm=False)
    mdc.gNames.update(dict((r.ID,r.Name) for r in gpr.dataLines))

varianceFile.close()
mdc.clusterFmt()
