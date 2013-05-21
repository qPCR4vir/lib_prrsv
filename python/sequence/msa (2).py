#!/usr/local/bin/python

###############################################################
#
# MSA.py
# This file contains the MSA class, which represents a multiple
# sequence alignment. In general it takes in a list of sequences
# and can perform various forms of alignment upon them, returning
# the output in various forms as well.
#
# Author: Dale Webster
# Date: 4/14/05
#
###############################################################


from Bio.Clustalw import MultipleAlignCL
from Bio import Clustalw
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import Translate

import fasta
from __init__ import *

import sys
import os
import math
import random

MAX_STOP_COUNT = 2

def makeFasta( title, seq ):
    "Quick helper for making a new FASTA record"
    fastaObj = fasta.Record()
    fastaObj.title = title
    fastaObj.sequence = seq
    return fastaObj

# Simple helper class
class Hit (object):
    def __init__(self, _strain, _start, _end):
        "constructor initializes fields."
        self.strain = _strain
        self.start = _start
        self.end = _end

class MSA (object):
    "Represents a Multiple Sequence Alignment"

    def __init__( self, clustalPath="clustalw" ):
        "Initializes default values"
        self.sequences = []
        self.alignments = []
        self.alignmentLength = -1
        self.tmpFileName = "tmpMSA.fasta"
	self.clustalPath = clustalPath

    def matrix(self):
        "Returns the 'alignment matrix', which is just a list of alignments"
        return self.alignments

    def length(self):
        "Returns the length of the alignment"
        return self.alignmentLength

    def count(self):
        "Returns the number of sequences in this aligment"
        return len( self.alignments )

    def loadFasta( self, fastaFileName ):
        "Given a FASTA file, parses and stores the sequences"
        self.sequences = list(fasta.iterator( fastaFileName ))
        
    def loadSequences( self, _sequences ):
        "Stores the given list of sequences"
        self.sequences = _sequences

    def loadAlignment( self, alignmentFile ):
        "Populates this object with the given alignment data from a CLUSTAL .aln file."

	# ***NOTE*** the CLUSTAL parser does not handle windows line breaks well...
        alignment = Clustalw.parse_file(alignmentFile)
        alignments = alignment.get_all_seqs()
        self.alignmentLength = alignment.get_alignment_length()
        
        for seq in alignments:
            sequence = fasta.Record()
            align = fasta.Record()

            sequence.title = seq.description
            align.title = seq.description

            align.sequence = seq.seq.tostring()
            sequence.sequence = seq.seq.tostring().replace("-","")

            self.alignments.append( align )
            self.sequences.append( sequence )

    def align(self):
        "Aligns the sequences using CLUSTAL, storing the results"

	if len(self.sequences) == 0:
		return

        self.sequencesToFile( self.tmpFileName )
        commandLine = MultipleAlignCL(os.path.join(os.curdir, self.tmpFileName), self.clustalPath)
        alignment = Clustalw.do_alignment(commandLine) 
        allRecords = alignment.get_all_seqs()
        length = alignment.get_alignment_length()
        
        alignmentStrings = []
        for record in allRecords:
            f = fasta.Record()
            f.title = record.description.strip()
            f.sequence = record.seq.tostring()
            alignmentStrings.append( f )

        self.alignments = alignmentStrings
        self.alignmentLength = length
            
        os.remove(self.tmpFileName)

    def sequencesToFile( self, filename ):
        "Saves the internal sequences to file"

        try:
            file = open(filename, 'w')
        except IOError:
            print "Error: No write permission in current directory for " + filename
            sys.exit(0)

        for record in self.sequences:
            file.write( str(record) + "\n" )

        file.close()

    def alignmentToPhylipFile( self, filename ):
        "Saves the alignment sequences to a Phylip file"

        nameLen = 10

        try:
            file = open(filename, 'w')
        except IOError:
            print "Error: No write permission in current directory for " + filename
            sys.exit(0)

        file.write( str(len(self.alignments)) + " " + str(len(self.alignments[0].sequence)) + "\n" )
        for record in self.alignments:

            length = min( len(record.title), nameLen )
            padding = nameLen - length
            
            file.write(str(record.title)[0:length] + " " * padding )
            file.write(str(record.sequence) + "\n")

        file.close()


    def alignmentToFile( self, filename, columns=80 ):
        "Given the number of horizontal chars, saves the alignment"

        try:
            file = open(filename, 'w')
        except IOError:
            print "Error: No write permission in current directory for " + filename
            sys.exit(0)

	file.write( self.clustalFormat( columns ) )
        file.close()


    def __repr__( self ):
        "Prints in CLUSTAL format with 80 columns per line."
	return self.clustalFormat()

    def clustalFormat( self, columns=80, labelLength=10 ):
	"Returns the string representaiton of this alignment formatted in the CLUSTAL .aln format."

        alignOffset = 16

	out = "CLUSTAL X (1.82) multiple sequence alignment\n\n\n"

        numRows = math.ceil( float(self.alignmentLength) / float(columns) )

        for i in range( int(numRows) ):
            low = i*columns
            high = min(low + columns, self.alignmentLength)
            for seq in self.alignments:
                out = out + seq.title[0:labelLength]
                padding = alignOffset - len(seq.title[0:labelLength])
                out = out + " " * padding
                out = out + seq.sequence[low:high] + "\n"

            out = out + " " * (alignOffset+columns+2) + "\n\n"

	return out


    def append( self, msa ):
        "Given another MSA, returns the result of appending it onto this MSA. The resulting MSA contains all common sequences."

        newSequences = self.appendSequences( self.sequences, msa.sequences )
        newAlignments = self.appendSequences( self.alignments, msa.alignments )
        newLength = self.alignmentLength + msa.alignmentLength

        newAlignment = MSA()
        newAlignment.sequences = newSequences
        newAlignment.alignments = newAlignments
        newAlignment.alignmentLength = newLength
        
        return newAlignment
        
    def appendSequences( self, sequenceListA, sequenceListB ):
        "Helper function takes in two lists of records, and returns the result of appending the sequences from A onto their matches in B."

        newSequences = []
        for sequenceA in sequenceListA :
            for sequenceB in sequenceListB:
                if sequenceA.title == sequenceB.title:
                    newSequence = fasta.Record()
                    newSequence.title = sequenceA.title
                    newSequence.sequence = sequenceA.sequence + sequenceB.sequence
                    newSequences.append( newSequence )

        return newSequences

    def backTranslate( self, nSequences ):
        "Given a list of nucleotide sequences to match the existing protein sequences, returns a back-translated version of the alignment."

        nDictionary = {}
        for seq in nSequences:
            nDictionary[seq.title.strip()] = seq.sequence

        pDictionary = {}
        for seq in self.sequences:
            pDictionary[seq.title.strip()] = seq.sequence

        alignDictionary = {}
        for seq in self.alignments:
            alignDictionary[seq.title.strip()] = seq.sequence

        sequenceNames = alignDictionary.keys()

        btAlign = MSA()
        alignments = []
        btAlign.sequences = nSequences
        for name in sequenceNames:
            newAlign = fasta.Record()
            newAlign.sequence = self.backTranslateHelp( nDictionary[name], pDictionary[name], alignDictionary[name] )
            newAlign.title = name
            alignments.append( newAlign )

        btAlign.alignments = alignments
        btAlign.alignmentLength = len( alignments[0].sequence )

        return btAlign

    def backTranslateHelp( self, nSeq, pSeq, alignSeq ):
        "Given a nucleotide, protein, and protein alignment sequence, returns the back-translated nucleotide alignment."

        newSequence = ""
        nCounter = 0
        pCounter = 0

        pMax = len(pSeq)
        nMax = len(nSeq)

        for i in range( len( alignSeq ) ):
            char = alignSeq[i]
            if char == "-":
                newSequence += "---"
            else:
                if pSeq[pCounter] == "X":
                    pCounter += 1
                    newSequence += "NNN"
                    nCounter += 3
                else:
                    pCounter+=1
                    newSequence += nSeq[nCounter:nCounter+3]
                    nCounter+=3

        return newSequence

    def alignUsingBackTranslation( self ):
        "Uses the internal sequences for alignment, first translating them to amino acid sequence for a better alignment."
        
        pSequences = []
        nSequences = []
        for sequence in self.sequences:
            pSequence = fasta.Record()
            pSequence.title = sequence.title
            pSequence.sequence = translate( sequence.sequence )
            # pSequence.sequence = pSequence.sequence.replace("*","X")
            if pSequence.sequence.count("*") > MAX_STOP_COUNT:
                print sequence.title + " : Too many stop codons. Check frame."
                print pSequence.sequence
                print sequence.title
                print sequence.sequence
                print
            else:
                pSequences.append( pSequence )
                nSequences.append( sequence )
        
        pMSA = MSA()
        pMSA.loadSequences( pSequences )
        pMSA.align()
        newAlign = pMSA.backTranslate( nSequences )

        self.sequences = newAlign.sequences
        self.alignments = newAlign.alignments
        self.alignmentLength = newAlign.alignmentLength

    def getConsensus( self ):
        "Returns a consensus sequence generated from the most numerous base in each column. If the higest two bases are represented the same number of times, the consensus base is chosen randomly."

        consensus = fasta.Record()
        consensus.title = "Consensus"
        sequence = ""

        matrix = self.matrix()
        for i in range( self.length() ):
            counts = { "A":0.0, "G":0.0, "C":0.0, "T":0.0, "-":0.0 }
            for row in matrix:
                rowSeq = row.sequence.upper()
                base = rowSeq[i]
                counts[base] += 1

            maxBase = "N"
            maxCount = -1
            for base in counts.keys():
                if counts[base]+random.random() > (maxCount+random.random()):
                    maxBase = base
                    maxCount = counts[base]

            if maxBase == "N":
                print "Internal Error."
                sys.exit()
            sequence += maxBase

        consensus.sequence = sequence
        return consensus

    def getAlignmentByName( self, name ):
        "Given the name of a sequence, returns its alignment string"

        for alignment in self.alignments:
            if alignment.title == name:
                return alignment.sequence

        return ""


    def getSubAlignment( self, xStart, xEnd ):
        "Returns a new MSA using the portion of this alignment specified by the start and end coordinates given."

        if( xEnd < xStart ):
            print "Error: xEnd < xStart in getSubAlignment."
            sys.exit(0)

        if( xEnd > self.length() ):
            print "Error: xEnd (" + str(xEnd) + ") > Length (" + str( self.length() ) + ") in getSubAlignment."
            sys.exit(0)

        if( xStart <= 0 ):
            print "Error: xStart <= 0 in getSubAlignment. ( index by 1, not 0 )"
            sys.exit(0)

        newMSA = MSA()

        for alignSeq in self.alignments:
            newAlignment = alignSeq.sequence[xStart-1:xEnd]
            newSequence = alignSeq.sequence[xStart-1:xEnd].replace("-","")
            newMSA.alignments.append( makeFasta(alignSeq.title, newAlignment) )
            newMSA.sequences.append( makeFasta(alignSeq.title, newSequence) )

        newMSA.alignmentLength = xEnd - xStart + 1

        return newMSA

    def findSequence(self, sequence):
        "Searches for the given sequence and it's reverse-compliment in the alignment. (ignores gaps) Returns a list of hits, where a hit consists of a strain, a start coordinate, and an end coordinate for the hit."

        sequenceRC = reverseComplement( sequence )

        forwardHits = self.findSequenceHelp( sequence )
        reverseHits = self.findSequenceHelp( sequenceRC )

        for reverseHit in reverseHits:
            start = reverseHit.start
            end = reverseHit.end

            reverseHit.start = end
            reverseHit.end = start

        return forwardHits + reverseHits

    def findSequenceHelp( self, sequence ):
        "Given a sequence, searches for an exact match in the alignment and returns all hits."
        hits = []
        
        for alignSeq in self.sequences:
            strain = alignSeq.title
            index = alignSeq.sequence.lower().find(sequence.lower())
            # print "strain: " + strain + " -> " + str(index)
            while index != -1:
                indexEnd = index + len(sequence)
                hits.append( Hit(alignSeq.title,
                            self.seqToAlignCoord( strain, index ),
                            self.seqToAlignCoord( strain, indexEnd )))
                index = alignSeq.sequence[index+1:].lower().find(sequence.lower())

        return hits


    def alignToSeqCoord( self, seqName, coord ):
        "Given a sequence name and an alignment coordinate, return the index of that position in the sequence."

        mySeq = ""
        for sequence in self.sequences:
            if sequence.title == seqName:
                mySeq = sequence.sequence

        myAlign = ""
        for alignment in self.alignments:
            if alignment.title == seqName:
                myAlign = alignment.sequence

        if coord > len( myAlign ) or coord < 0:
            print "Error: alignToSeqCoord received bad coordinate: " + str(coord)
            sys.exit(0)

        seqCoord = 0
        alignCoord = 0

        while alignCoord < coord:
            if myAlign[alignCoord] == "-":
                alignCoord += 1
            else:
                seqCoord += 1
                alignCoord += 1

        # print "seqToAlignCoord for " + seqName + " returning " + str(alignCoord-1)
        return seqCoord



    def seqToAlignCoord( self, seqName, coord ):
        "Given a sequence name and a coordinate, return the index of that position in the alignment."

        mySeq = ""
        for sequence in self.sequences:
            if sequence.title == seqName:
                mySeq = sequence.sequence

        myAlign = ""
        for alignment in self.alignments:
            if alignment.title == seqName:
                myAlign = alignment.sequence

        if coord > len( mySeq ) or coord < 0:
            print "Error: seqToAlignCoord received bad coordinate: " + str(coord)
            sys.exit(0)

        seqCoord = 0
        alignCoord = 0

        while seqCoord <= coord:
            if myAlign[alignCoord] == "-":
                alignCoord += 1
            else:
                seqCoord += 1
                alignCoord += 1

        # print "seqToAlignCoord for " + seqName + " returning " + str(alignCoord-1)
        return alignCoord - 1


    def basicMatrix(self):
        "Returns the nuc/aa letters in a python matrix"

        matrix = []
        rows = self.count()
        columns = self.length()

        for i in range(rows):
            matrix.append([])
            sequence = self.alignments[i].sequence
            for j in range(columns):
                matrix[i].append( sequence[j] )

        return matrix


    def translate( nSequence ):
        "Returns the amino acid translation of the given nucleotide sequence."
  
        myAlpha = IUPAC.unambiguous_dna
        mySequence = Seq(nSequence, myAlpha)
        
        standardTranslator = Translate.unambiguous_dna_by_id[1] 
        
        myProteinSequence = standardTranslator.translate(mySequence)
        
        return myProteinSequence.tostring()
