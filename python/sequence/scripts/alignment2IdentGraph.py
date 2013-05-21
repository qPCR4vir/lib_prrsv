#!/usr/local/bin/python

###############################################################
#
# This script takes in a .aln file, and computes a graph
# of average pair-wise identity scores using an N window-size across
# the length of the alignment.
#
# Dale Webster 9/02/05
#
###############################################################

from sequence.MSA import MSA
import sys
import os
import math


# Given an alignment and flanking coordinates, calculates the average percent identity across all pairwise comparisons.
def computeIdentity( alignMSA, start, stop, reference ):
    matrix = alignMSA.matrix()

    seqList = []
    for row in matrix:
        seqList.append( row )

    same = 0.0
    total = 0.0001
    for seq1 in seqList:
        if reference == "-" or reference == seq1.title:
            for seq2 in seqList:
                for i in range( start, stop ):
                    if seq1.title != seq2.title:
                        if seq1.sequence[i] != "-" and seq2.sequence[i] != "-":
                            total += 1
                            if seq1.sequence[i] == seq2.sequence[i]:
                                same += 1

    return same/total




referenceSeq = "-"

if len(sys.argv) >= 3:
    alignFile = sys.argv[1]
    windowSize = int(sys.argv[2])
    if len(sys.argv) == 4:
        referenceSeq = sys.argv[3]
else:
    print "Usage: ./alignment2IdentGraph.py alignmentFile window-size [reference sequence] > outputFile.csv"
    sys.exit(0)

# Load the alignment
alignMSA = MSA()
alignMSA.loadAlignment( alignFile )

alnMatrix = alignMSA.matrix()
length = len( alignMSA.alignments[0].sequence )

# Default window size is the full alignment. (Returns a single value for the APWI for the entire alignment)
if windowSize == -1:
    windowSize = length-2

flank = int(windowSize/2)
cursor = flank

identGraph = []

while cursor < (length - flank - 1):
    identity = computeIdentity( alignMSA, int(cursor-flank), int(cursor+flank), referenceSeq )
    identGraph.append( identity )
    cursor += 1

for i in range( len(identGraph)-1 ):
    print str(int(identGraph[i]*1000)/10.0) + "\n",

print str(int(identGraph[ len(identGraph)-1 ]*100))


