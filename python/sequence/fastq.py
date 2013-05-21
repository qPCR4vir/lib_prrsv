"""fastq
sequence utilities

kael and dale

$Id: fastq.py,v 1.11 2011/05/03 21:07:11 julio Exp $

"""
__version__ ="$Revision: 1.11 $"

import re
import os

from types import *

from __init__ import *
from utils import flatten,multiFile,fileIterator,getIteratable

import fasta

#
# Fastq Record Class incorporating quality values.
# Built on top of the fasta.Record class.
#
class Record (fasta.Record):
    """The Fastq Record class.
    
    Members:
    title       Title line ('@' character not included).
    sequence    The sequence.
    quality     The quality values.

    Record acts like a string in these cases:
    len(rec)
    rec[x:y]

    str(rec) returns a valid qfastq formatted string

    """

    def __init__ (self,title='',sequence='',colwidth=60,quality=[]):
        """Create a new Record.  colwidth specifies the number of residues
        to put on each line when generating qfastq format. Quality is an array
        of integers representing the quality of each base.
        """
        fasta.Record.__init__(self,title=title,sequence=sequence,colwidth=colwidth)
        self.quality=quality

    def __str__(self):
        s = []
        s.append('@%s' % self.title)
        i = 0
        while i < len(self):
            s.append(self[i:i+self.colwidth][0])
            i = i + self.colwidth
        s.append('+%s' % self.title)
        i = 0
        while i < len(self.quality):
            s.append( "".join( intToQual( self.quality[i:i+self.colwidth] ) ) )
            i = i + self.colwidth
        return "\n".join(map(str,s))

    def __getitem__(self,item):
        return (self.sequence[item], self.quality[item])
   
    def split(self,sliceBases=100000):
        """Returns an iterator of slices of the record
        """
        n=0
        for start in range(0,len(self.sequence),sliceBases):
            rec = Record()
            rec.title = self.title 
            rec.sequence = self.sequence[start:start+sliceBases]
            rec.quality = self.quality[start:start+sliceBases]
            rec.slice=n
            yield rec
            n+=1

    def fasta( self ):
        """Returns a fasta.Record version of this fastq.Record.
        Basically removes the quality scores but allows you
        to use the fasta.Record specific functions."""

        return fasta.Record( title=self.title,sequence=self.sequence,colwidth=self.colwidth )

    def prune( self, quality, harsh=False ):
        """Prunes bases from the 3' end of this read until a
        base with 'quality' or higher quality score is found.
        Returns a pruned version of this record. Does not modify
        this Record.

        Setting harsh = True will start from the 5' end, and prune
        everything after the first base that has lower quality than
        the specified value."""

        if not harsh:
            seq = self.sequence
            qual = self.quality
            while len(qual) > 0 and qual[-1] < quality:
                qual = qual[:-1]
                seq = seq[:-1]
        else:
            seq = ""
            qual = []
            for i in range(len(self.sequence)):
                if self.quality[i] >= quality:
                    seq += self.sequence[i]
                    qual.append( self.quality[i] )
                else:
                    break

        return Record(title=self.title,
                      sequence=seq,
                      colwidth=self.colwidth,
                      quality=qual )

    def trimEnds(self,quality):
        """Prunes bases from the ends of this read until a
        base with 'quality' or higher quality score is found.
        Moodify this Record."""
        for i,q in enumerate(self.quality):
            if q>=quality:
                if i>0:
                    self.quality=self.quality[i:]
                    self.sequence=self.sequence[i:]
                break
        
        for i,q in enumerate(reversed(self.quality)):   
            if q>=quality:
                if i>0:
                    self.quality=self.quality[:(-1*i)]
                    self.sequence=self.sequence[:(-1*i)]
                break
        
    

    def lowQ2N(self,quality,screenChar='N'):
        """replace low quality bases with N (or other screenChar).
        returns number of screened bases.
        """
        ct=0
        sChars=list(self.sequence)
        for i,c in enumerate(sChars):
            if self.quality[i] < quality:
                ct+=1
                sChars[i]=screenChar[0]
        self.sequence=''.join(sChars)
        return ct

    def minQ (self):
        """return minimum quality value
        """
        return min(self.quality)

    def maxQ (self):
        """return maximum quality value
        """
        return max(self.quality)

    def meanQ (self):
        """return mean quality value
        """
        return float(sum(self.quality))/float(len(self.quality))

    def medianQ (self):
        """return median quality value
        """
        import pylab 
        return pylab.median(self.quality)

    def phred33transform(self):
        """Add 33 to raw parsed quality.
        Raise ValueError if transformation would
        yield a qulaity greater than 42.
        """
        if max(self.quality)>9:
            raise ValueError, "Quality scores dont look like phred33"
        self.quality=[q+33 for q in self.quality]
        

def qualToInt( quals ):
    """Given one or more quality characters, returns the corresponding
    list of integer values, as defined by Solexa.
    """
    return map( lambda x: ord(x)-64, quals )

def intToQual( ints ):
    """Given a list of integers, returns the corresponding
    quality string as defined by Solexa."""
    return "".join( map( lambda x: chr(x+64), ints ) )

#
# Identifier
# 
def looksLikeFastq( path ):
    """Returns true if the given file handle appears to be a Fastq file.
    DOES NOT VALIDATE THE FILE. JUST CHECKS THE FIRST RECORD TO SEE IF
    IT LOOKS LIKE A FASTQ RECORD."""
    try:
        for record in FastqIterator( fh ): # TODO - BROKEN? Wrong variable?
            return True
    except:
        return False
    
#
# Iterator
#
def FastqIterator(files,raw=False,titleSet=None):
    """return an iterator of Records found in file handle, fh.
    if records are not needed raw can be set to True, and then 
    you can get (titleStr, seqStr, qualityStr).  With raw output,
    the sequence and quality strings have the newlines still in them.
    """
    def readTotitle(fh, titleChar):
        """returns a tuple ([lines before the next title line], next tile line)
        """
        preLines = []
        while True:
            l = fh.readline()
            if l.startswith(titleChar):
                return (preLines,l)
            elif l == '':
                return preLines,None
            else:
                preLines.append(l)

    
    for fh in fileIterator(files):
        preLines,nextTitleLine =readTotitle(fh,'@')

        while nextTitleLine != None:
            seqTitle = nextTitleLine[1:].rstrip()
            preLines,nextTitleLine=readTotitle(fh,'+')
            qualTitle = nextTitleLine[1:].rstrip()
            if len(qualTitle.strip()) > 0 and seqTitle != qualTitle:
                raise FastqParseError, ("Error in parsing: @title sequence entry must be immediately "
                                        "followed by corresponding +title quality entry.")
            seqLines = preLines
            qualLines = []
            for i in range(len(seqLines)): # Quality characters should be the same length as the sequence
                qualLines.append( fh.readline() )

            preLines,nextTitleLine=readTotitle(fh,'@')

            if titleSet!= None and seqTitle not in titleSetSet:
                continue

            seqLines = map(lambda x: x.strip(), seqLines)
            qualLines = map(lambda x: x.strip(), qualLines)
            if raw:
                yield (seqTitle, ''.join(seqLines), ''.join(qualLines))
            else:
                rec=Record()
                rec.title=seqTitle
                rec.sequence=''.join(seqLines)
                rec.quality=flatten(map(lambda x: qualToInt(x),qualLines))
                yield rec

iterator=FastqIterator


def phdQualIterator(fastaFiles, qualFiles, raw=False):
    def readTotitle(fh):
        """returns a tuple ([lines before the next title line], next tile line)
        """
        preLines = []
        while True:
            l = fh.readline()
            if l.startswith('>'):
                return (preLines,l)
            elif l == '':
                return preLines,None
            else:
                preLines.append(l)

    def qualityIterator(filename):
        fh = file(filename)
        preLines, nextTitleLine = readTotitle(fh)

        while nextTitleLine != None:
            title = nextTitleLine[1:].rstrip()
            preLines, nextTitleLine = readTotitle(fh)

            yield (title, ' '.join(preLines))

    qualFiles = getIteratable(qualFiles)
    for idx, fastaFh in enumerate(fileIterator(fastaFiles)):
        qualIter = qualityIterator(qualFiles[idx])

        for seqTitle, sequence in fasta.FastaIterator(fastaFh, raw=True):
            qTitle, qualities = qualIter.next()
            sequence = sequence.replace('\n', '')
            qualities = qualities.replace('\n', '')
            if raw:
                yield (seqTitle, ''.join(sequence), ''.join(qualities))
            else:
                qualities = qualities.split()
                if len(sequence) != len(qualities):
                    raise Exception, 'Invalid number of qualities'
                rec = Record()
                rec.title = seqTitle
                rec.sequence = sequence
                rec.quality = qualities
                yield rec
            

