#!/usr/bin/env python

import sys

sys.path.append('..')

from AOS4.Oligo import Oligo
from AOS4.Parent import Parent
from AOS4.Score import GCcontent
from AOS4.ScoreCard import ScoreCard

fastaFilePath = 'test_input'
oLength = 70

GCnoSC = GCcontent()

seq=''
pid = None
fFile = file(fastaFilePath)
for line in fFile:
    if line[0] == '>':
        if seq == '':
            pid = line[1:].rstrip() 
        else:
            break
    elif pid != None:
        seq += line.rstrip()


parent = Parent(pid=pid,seq=seq)
print "Parent made"
oligos = parent.makeOligos(oLength)
print "Oligos made"
print GCnoSC(oligos[0])
print "GC calculated without scorecard"
sc = ScoreCard(oligos)
print "ScoreCard made"
sc.addScore(GCcontent)
print "ScoreCard - added GC content"
print "scorecard value for same oligo"
print oligos[0].score((oligos[0].scIndexes.keys()[0].scores.keys()[0]))
