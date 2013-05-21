#!/usr/local/bin/python

from getopt import getopt
import os
import os.path
import sys

wrdSize = 16


opts,args = getopt(sys.argv[1:],'W:')
for o,a in opts:
    if o=='-W':
        wrdSize=int(a)
    else:
        print '-W is the only option!'
        sys.exit(1)

dbPath = args[0]

bCmd = 'megablast -f -R -d %s -W %s -D3 -i %%s -o %%s ' % (dbPath,wrdSize)
clobber = False

gridItpairs = []

for q in args[1:]:
    #make mb output file name
    # in file like: s_2_AA.hsEland_NM.fasta
    # out file like s_2_AA..hsEland_NM.12046.mbr

    qName,qExt=os.path.splitext(q)

    if qExt in ('.fa','.fasta'):
        o=qName + '.%s.mbr' % wrdSize
    else:
        o=q + '.%s.mbr' % wrdSize

    if os.access(o,os.F_OK) and not clobber:
        continue
    
    o=os.path.split(o)[1]
    gridItpairs.append((q,o))


cmd = ("gridIt.py \'%s\' "% (bCmd) +
           "\'%s\'" % "\' \'".join([' '.join(pair) for pair in gridItpairs]))

print "running:\n" + cmd
os.system (cmd)



    
