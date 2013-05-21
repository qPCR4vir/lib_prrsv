"""fasta
sequence utilities

kael and dale

$Id: fasta.py,v 1.92 2012/12/17 22:06:29 kael Exp $

"""
__version__ ="$Revision: 1.92 $"

import commands
import copy
import subprocess
import time
import re
import string
import os
import sys
import tempfile
import shutil
import bisect
import random
from types import *

from utils import *

from __init__ import *

#import kdbom.exceptions 


#
# FASTA regular expressions
#
titlePat = re.compile(r'^>(?P<titleStr>.*)$',re.MULTILINE)
fTitlePat = re.compile(r'^>',re.MULTILINE)
seqPat = re.compile(r'^>(?P<titleStr>.*)$(?P<rawSeqStr>.*)>|',
                    (re.MULTILINE|re.DOTALL))
ReType = type(seqPat)

OS_READ_MAX = 2147483642

#
# re helpers
#

def titleMatch(rec,exps):

    for pattern in reTransformer(exps):
        mObj = pattern.search(rec.title)
        if mObj != None:
            return True
    return False



#
# Fasta Record Class
# For the timebeing this is derived from BioPython
# For heaven's sake don't use the BioPython parser!
#
class Record (StringSequence):
    """The Fasta Record class.
    
    Members:
    title       Title line ('>' character not included).
    sequence    The sequence.

    Record acts like a string in these cases:
    len(rec)
    rec[x:y]

    str(rec) returns a valid fasta formatted string

    """

    def __init__ (self,title='',sequence='',colwidth=60):
        """Create a new Record.  colwidth specifies the number of residues
        to put on each line when generating FASTA format.
        """
        #Fasta.Record.__init__(self,colwidth=colwidth)
        self.colwidth=colwidth
        self.title=title
        self.sequence=sequence

    def __str__(self):
        s = []
        s.append('>%s' % self.title)
        i = 0
        while i < len(self):
            s.append(self[i:i+self.colwidth])
            i = i + self.colwidth
        return os.linesep.join(s)
   

    def gi (self):
       
        """attempt to deuce gi from the title.
        """
        ff = self.fastaFields()
        if 'gi' in ff:
            return int(ff['gi'])
        else:
            return None

    def appendField(self,name,value):
        """
        """
        self.title += '|%s|%s' %(name,value)
        
       
    def fastaFields(self):
        """
        """
        rv = {}
        f = self.title.split('|')
        start = len(f)%2
        if f[0] == 'gi':
            start = 0
        for i in range(start,len(f),2):
            try:
                rv[f[i]] = f[i+1]
            except:
                rv[f[i]]=''
        return rv

       
    def m8name(self):
        return self.title.split()[0]

    def m8annotation(self):
        fields = self.title.split(None,1)
        if len(fields) == 2:
            return fields[1]
        else:
            return ''

    def split(self,sliceBases=100000):
        """Returns an iterator of slices of the record
        """
        n=0
        for start in range(0,len(self.sequence),sliceBases):
            rec = Record()
            rec.title = self.title 
            rec.sequence = self.sequence[start:start+sliceBases]
            rec.slice=n
            yield rec
            n+=1
        

    def DRNAsubrecords (self,titleCallback=None):
        """Returns an iterator of records made from the DRNA(i.e. not XXX)
        stretches of this record.

        If titleCallback is provided it should take this, the parent, record
        and the nucleotide offset of the subrecord (sequence starts at 1
        (not 0) per biology.

        The default title callback is giTileTitle.
        """

        if titleCallback == None:
            titleCallback = giTileTitle
        
        for subMatch in self.reSearch(DRNAbasesPat):
            rec = Record(
                title=titleCallback(self,subMatch.start()+1),
                sequence=subMatch.group(0))
            yield rec

    def shortySubrecords (self, shortyLength=50):
        '''finds shorties <= shortyLength and yields them as fasta records
        '''
                   
        for start,end in self.locateShorties(shortyLength=shortyLength):
            rec = Record(
                title=generalTileTitle(self,start+1),
                sequence=self[start:end])
            yield rec

    def tile(self, tileLength, startOffset=None,
             minLength=None,allDRNA=False, titleCallback=None,
             returnEndTile=False,colwidth=60):
        """simple tiling of current record, returns a generator.

        tileLength is an integer (all tiles are the same size).

        startOffset can be an integer specifying the tiles'
        start position relative to the prior tile's start, or
        it can be None for end-to-end tiling.

        minLength specifies the smallest of the tail records that
        will be return, if None, tileLength is used.

        titleCallback requirements are the same as for DRNAsubrecords.

        If allDRNA is true no tiles with invalid or degenerate bases
        will be generated.  The tile ennumeration begins after stretches
        of such bases.

        chopper is a historical alias for tile.
        """

        if titleCallback == None:
            titleCallback = giTileTitle

        if startOffset == None:
            startOffset = tileLength

        if minLength == None:
            minLength = tileLength

        if len(self.sequence) >= minLength:

            startPos=range(0,len(self)-(minLength)+1,startOffset)

            if returnEndTile:
                endStartPos = len(self)-tileLength
                if endStartPos not in startPos:
                    startPos.append(endStartPos)


            for start in startPos:
                yield Record(
                    title=titleCallback(self,start+1),
                    sequence=self.sequence[start:start+tileLength],
                    colwidth=colwidth
                    )

    # alias in honor of dave wang's C program
    chopper = tile
    
    def titleMatch(self,patterns):
        """True if title matches any of the given patterns
        """ 
        return titleMatch(self,patterns)

    def fasta(self):
        """Return self
        """
        return self


            
def giTileTitle (rec,n):
    """Make titles for tile based on gi and start of tile
    (sequence starts at 1 (not 0) per biology."""
    # has the record already been mangled?
    mangleMatch=re.match('(?P<gi>[0-9]+)_nt(?P<pos>[0-9]+)',rec.title)
    if mangleMatch != None:
        recOffset = int( mangleMatch.group('pos'))-1
        nStrLength = str(len(mangleMatch.group('pos')))
        base = mangleMatch.group('gi')

    else:
        recOffset=0
        # the string rep of the length of the longest n for this record
        nStrLength = str(len(str(len(rec))))
        
        base = rec.gi()
        if base == None:
            base = rec.m8name()

    fmtStr = "%s_nt%0" + nStrLength +"d"
    return fmtStr % (base,n+recOffset)

def generalTileTitle (rec,n):
    """Make titles for tile based current title and
    start of tile.
    Sub record format is recognized adjust start position
    using subrecord offset.
    (sequence starts at 1 (not 0) per biology."""
    # has the record already been mangled?
    mangleMatch=re.match('^(?P<prefix>.+)_nt(?P<pos>[0-9]+).?$',rec.title)
    if mangleMatch != None:
        recOffset = int( mangleMatch.group('pos'))-1
        nStrLength = str(len(mangleMatch.group('pos')))
        base = mangleMatch.group('prefix')

    else:
        recOffset=0
        # the string rep of the length of the longest n for this record
        nStrLength = str(len(str(len(rec))))
        
        base = rec.gi()
        if base == None:
            base = rec.m8name()

    fmtStr = "%s_nt%0" + nStrLength +"d"
    return fmtStr % (base,n+recOffset)


#
# Iterator
#
def FastaIterator(files,raw=False,titleSet=None):
    """return an iterator of Records found in files.
    If records are not needed raw can be set to True, and then 
    you can get (titleStr, seqStr).  With raw output, the sequence 
    string has the newlines still in it.

    If titleSet is not None it should support the 'in' opperation,
    and the records returned will be limited to those with tiles
    that match 
    """
    tSet = copy.copy(titleSet)
    
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

    for fh in fileIterator(files):
        preLines,nextTitleLine =readTotitle(fh)

        while nextTitleLine != None:
            title = nextTitleLine[1:].rstrip()
            preLines,nextTitleLine=readTotitle(fh)

            if tSet != None and title not in tSet:
                continue

            if raw:
                yield (title,''.join(preLines))
            else:
                rec=Record()
                rec.title=title
                rec.sequence=''.join(map(lambda x: x.rstrip(),preLines))
                rec.sequence=rec.dropNonLetters()
                yield rec
            if tSet != None:
                tSet.remove(title)
                if len(tSet) == 0:
                    break

iterator=FastaIterator

def fastaDictionary(inFile, keyFcn = lambda fastaRec: fastaRec.m8name() ):
    """return a dictionary of fasta Records (values) taken from inFile.
    dictionary's keys are generated using keyFcn, which must take one
    argument, the Record.  By default the keys are the m8names of the
    records.
    """
    rd = {}
    for rec in FastaIterator(inFile):
        rd[keyFcn(rec)] = rec
    return rd

FastaDictionary = fastaDictionary

#
# File Size Calculation Functions
#

def generalIterator(files,**kwds):
    """returns a fasta or fastq iterator as appropriate for each
    file in files. 
    """

    for fh in fileIterator(files):
        l = fh.readline()
        while len(l) == 0:
            l=fh.readline()
        fh.seek(0)
        if l.startswith('>'):
            g=FastaIterator(fh,**kwds)
        elif l.startswith('@'):
            import fastq
            g=fastq.FastqIterator(fh,**kwds)
        else:
            raise RuntimeError, "Sequence iteration error: can't determine type of file"
        for rec in g:
            yield rec
            
def fastaCount(things):
    """Count the number of titles in a file-like things.

    file like things can be one of more (in a list of tuple)
    file objects or paths .
    """
    if type(things) in (ListType,TupleType):
        return map(fastaCount,things)
    if type(things) in StringTypes:
        things = open(things)
    c = 0
    chunkSize = 10**4  # lowered - no big change in speed from 10^9.[kf]
    chunk = things.read(chunkSize)
    while len(chunk) > 0:
        c += len(fTitlePat.findall(chunk))
        chunk = things.read(chunkSize)
    return c


def fastaTitles(things):
    """Return the titles in file like things.

    file like things can be one of more (in a list of tuple)
    file objects or paths .
    """
    if type(things) in (ListType,TupleType):
        return map(fastaTitles,things)
    return titlePat.findall(toString(things))


def countRecsAndBases(things):
    """return the total # of fasta records and the number of
    valid DNA and RNA bases in fasta file-like things.

    file like things can be one of more (in a list of tuple)
    file objects or paths .
    """
    if type(things) in (ListType,TupleType):
        manyRslts = map(countRecsAndBases,things)
        tRecs = sum(map(lambda x: x[0],manyRslts))
        tBases = sum(map(lambda x: x[1],manyRslts))
        return (tRecs,tBases)
        
    recCount  = 0 
    baseCount = 0
    for t,s in FastaIterator(things,raw=True):
        recCount += 1
        baseCount += DRNABaseCount(s)
    
    return (recCount,baseCount)

def baseComposition(things,caseSensitive=False):
    """return the total of occurrences of each 'base' over every sequence
    in fasta file-like things, a dictionary is returned for each file.
    The dict. keys are the bases/letters (don't have to be valid D/RNA)
    values are the # of times seen.

    file like things can be one of more (in a list of tuple)
    file objects or paths .
    """
    if type(things) in (ListType,TupleType):
        return map(fastaCount,things)

    else:
        counts={}
        for t,s in FastaIterator(things,raw=True):
            for letter in s:
                if letter not in counts:
                    counts[letter]=0
                counts[letter]+=1

        for k in counts.keys():
            if k not in string.letters:
                del counts[k]
            else:
                if not caseSensitive:
                    if k!=k.upper():
                        K=k.upper()
                        if K not in counts:
                            counts[K] = 0
                        counts[K] = counts[k]
                        del counts[k]
                        

        return counts
        
   
    
#
# Fasta File Filtering and Partitioning
#
def fastaPartition(inFile,trueFile,notTrueFile,callback,
                   recordMassage=None,clobber=False):
    """partition a fasta file based on callback(Record)'s return
    value.

    inFile must be a file object or a valid path.  trueFile and
    notTrueFile may be file objects paths or None if coresponding records
    are to be discarded.

    if recordMassage is set to a function that takes a recrod as an argument,
    that function is called before the partioning callback. 

    Overwriting of output files is controled by the clobber flag.

    returns (trueCount,notTrueCount)
    """

    tCount = 0
    ntCount = 0

    if type(inFile) in StringTypes:
        inFile = file(inFile)

    if type(trueFile) in StringTypes:
        if os.access(trueFile,os.F_OK) and not clobber:
            raise RuntimeError, "file: %s exists" %trueFile
        trueFile=file(trueFile,'w')

    if type(notTrueFile) in StringTypes:
        if os.access(notTrueFile,os.F_OK) and not clobber:
            raise RuntimeError, "file: %s exists" % notTrueFile
        notTrueFile=file(notTrueFile,'w')

    for rec in FastaIterator(inFile):
        if recordMassage != None:
            recordMassage(rec)
        if callback(rec):
            tCount +=1
            if trueFile !=None:
                trueFile.write(str(rec))
                trueFile.write('\n')
        else:
            ntCount+=1
            if notTrueFile != None:
                notTrueFile.write(str(rec))
                notTrueFile.write('\n')

    return tCount,ntCount

def m8partition(inFile,foundOut,notFoundOut,m8Files,
                qryMatch=False,subjMatch=False,
                minScore=0,maxE=100000,clobber=False):
    """partition a fasta file into 2 new files based on presence
    of fasta title matching either query or subject in blast (m8/D3)
    output.  Liberal minScore and maxE values may be overridden.
    """
    import blastNoSQL as B
    
    if int(qryMatch) + int(subjMatch) != 1:
        raise ValueError, "Either qrySet or subjSet must be True"

    if type(foundOut) in StringTypes:
        foundOut = safeOFW(foundOut,clobber)
    if type(notFoundOut) in StringTypes:
        notFoundOut = safeOFW(notFoundOut,clobber)

    tSet = B.m8TitleSet(m8Files,qrySet=qryMatch,subjSet=subjMatch,
                        minScore=minScore,maxE=maxE)

    fastaPartition(inFile,foundOut,notFoundOut,lambda rec: rec.title in tSet)
    

def bowtiePartition(inFile,foundOut, notFoundOut,bowtieFiles,
                    subjMatch=False,qryMatch=False,clobber=False):
    """partition a fasta file into 2 new files based on presence
    of fasta title matching either query or subject in bowtie
    output.
    """
    import bowtie as B
    
    if int(qryMatch) + int(subjMatch) != 1:
        raise ValueError, "Either qrySet or subjSet must be True"

    if type(foundOut) in StringTypes:
        foundOut = safeOFW(foundOut,clobber)
    if type(notFoundOut) in StringTypes:
        notFoundOut = safeOFW(notFoundOut,clobber)

    tSet = B.readTitleSet(bowtieFiles,qrySet=qryMatch,subjSet=subjMatch)

    fastaPartition(inFile,foundOut,notFoundOut,lambda rec: rec.title in tSet)
     

def qualityTrim(sequence, qualities, badQualityTrim=(None, None),
                badQualityLimit=None ):
    """Returns (`sequence`, `qualities`) after trimming them based on the
    `badQualityLimit`.
    
    Arguments:
    - `sequence`: 
    - `qualities`:
    - `badQualityTrim`: Defines how many bad adjacent qualities must
    be found to trim on. It's a tuple to allow trimming in one or
    both sides of the sequence. For example: if `badQualityTrim` is
    (None, x) then `sequence` will only be trimmed on the 3'.
    - `badQualityLimit`: Limit start considering a quality as bad.
    """
    if badQualityLimit is None or \
       badQualityTrim == (None, None):
        return (sequence, qualities)

    if min(qualities) >= badQualityLimit:
        return (sequence, qualities)

    m = len(qualities)/2
    totalCount = 0
    if badQualityTrim[0] is not None:
        totalCount += badQualityTrim[0]
    if badQualityTrim[1] is not None:
        totalCount += badQualityTrim[1]
    totalCount -= 1
    p5idx = m - totalCount
    p3idx = m + totalCount
    t,c = None, 0
    l = min(filter(lambda x: x is not None, badQualityTrim))
    for i in xrange(p5idx, p3idx+1):
        if qualities[i] >= badQualityLimit:
            t = None
            c = 0
        else:
            if t is None:
                t = i
            c += 1
            if c >= l:
                break
    if c >= l:
        qualities = qualities[:t]
        sequence = sequence[:t]
        m = len(qualities)/2
    
    trim5, trim3 = None, None
    half5 = qualities[:m]
    half5.reverse()
    half3 = qualities[m:]

    isHalf5 = True
    for half, badCount in ((half5, badQualityTrim[0]),
                           (half3, badQualityTrim[1])):
        if badCount is None:
            if isHalf5:
                isHalf5 = False
            continue

        trimIdx = None
        curCount = 0
        for idx, qual in enumerate(half):
            if qual >= badQualityLimit:
                trimIdx = None
                curCount = 0
            else:
                if trimIdx is None:
                    trimIdx = idx
                curCount += 1
        
                if curCount >= badCount:
                    if isHalf5:
                        trim5 = trimIdx
                    else:
                        trim3 = trimIdx
                    break
        isHalf5 = not isHalf5

    seqHalf5 = sequence[:m]
    seqHalf3 = sequence[m:]
    
    if trim5 is not None:
        seqHalf5 = seqHalf5[-trim5:]
        half5 = half5[:trim5]
    half5.reverse()

    if trim3 is not None:
        seqHalf3 = seqHalf3[:trim3]
        half3 = half3[:trim3]

    return (seqHalf5 + seqHalf3, half5 + half3)


def filterFasta(fileNames,removeRepeats=False,
                minLZWsize=0,minLZWratio=0.0,
                regexRequire=None, regexRequireName='regex',
                regexProhibit=None, regexProhibitName='not_regex',
                regexName=None, badQualityTrim=(None, None),
                badQualityLimit=None, minTrimLength=0,
                applyRegexSequence=True, applyRegexTitle=False,
                maxHomo=0,fastqOutput=False, forceQuality=None,
                rejectQuality=(None,None),fivePTrim=0,threePTrim=0,
                barcodes=[],output=None, DEBUG=False):
    
    """ Perform several filtering steps on one or more fasta files.
    Outputfiles are named by inserting the filtering options before
    the extention (of any) of the input file name.

    If a regex (compared to sequence) is required or prohibited, a
    descriptive regexRequireName and/or regexProhibitName, can and should
    be given as well.

    Possible opperations:
       -Filter on LZW parameters (minLZWratio, minLZWsize)
       -remove duplicate records (removeRepeats) *see below
       -require a regex match to the sequence, title, or both. (regexRequire) 
       -require a regex no match to the sequence, title, or both. (regexProhibit)
       -require that the longest homopolymeric run be no longer that maxHomo (maxHomo)

    *Notes about Remove Duplicates
       -Quality restrictions specified with remove duplicates will result in
        artifacts based on different fastq qualities for different reads of
        the same sequence.
       -Duplicate removal over more than one file:  If there is more than one
        input file and duplicates are removed they are removed over the entire
        set of files.
    """
    import fastq
    if DEBUG:
        print "*** filterFasta ***\n Library versions:"
        print __file__, __version__
        import aos
        for l in (aos,):
            print l.__file__, l.__version__

        print "libaos: ", aos.__aos

    if type(fileNames) in StringTypes:
        fileNames = [fileNames]

    if regexProhibitName =='not_regex'  and regexName != None:
        regexProhibitName=regexName

    if regexRequireName == 'regex' and regexName != None:
        regexRequireName=regexName

    for fn in fileNames:
        filesToClean=[]
        try:

            if threePTrim>0 or fivePTrim>0 or len(barcodes)>0 or removeRepeats:
                if threePTrim > 0:
                    threePTrimIndex =-threePTrim
                else:
                    threePTrimIndex = None

                # first make number of barcodes + 1 temp files with suffix barcode
                barcodes.append('')
                bcSeqHashes = dict(zip(barcodes,
                                       [{} for x in range(len(barcodes))]))
                #check for dups
                for rec in generalIterator(fn):
                    ct=rec.countFromTitle()
                    if ct == None:
                        ct =1
                    # trim off 5' and 3' bases
                    seq = rec.sequence[fivePTrim:threePTrimIndex]
                    if hasattr(rec,'quality'):
                        qual = rec.quality[fivePTrim:threePTrimIndex]
                    curBC = seq[:len(barcodes[0])]
                    if curBC in barcodes:
                        seq = seq[len(curBC):]
                        if hasattr(rec,'quality'):
                            qual = qual[len(curBC):]
                        else:
                            qual = [99] * len(seq)
                    else:
                        curBC = ''
                    seq, qual = qualityTrim(seq, qual, badQualityTrim,
                                            badQualityLimit)
                    seqHash = hash(seq)
                    if seqHash not in bcSeqHashes[curBC]:
                        bcSeqHashes[curBC][seqHash]= 0
                    bcSeqHashes[curBC][seqHash]+=ct 

                    

                bcFileDict = dict(zip(barcodes,
                                      [mystemp(suffix = '.'+bc)
                                       for bc in barcodes]))
               
                filesToClean.extend([x[1] for x in bcFileDict.values()])    
                #write em out
                for rec in generalIterator(fn):
                    # trim off 5' and 3' bases
                    rec.sequence = rec.sequence[fivePTrim:threePTrimIndex]
                    if hasattr(rec,'quality'):
                        rec.quality = rec.quality[fivePTrim:threePTrimIndex]
                    else:
                        rec.quality = [99]*len(rec.sequence)
                    curBC = rec.sequence[:len(barcodes[0])]
                    if curBC in barcodes and curBC != '':
                        rec.title += "_" + curBC
                        rec.sequence = rec.sequence[len(curBC):]
                        if hasattr(rec,'quality'):
                            rec.quality = rec.quality[len(curBC):]
                    else:
                        curBC = ''
                    seq, qual = qualityTrim(rec.sequence, rec.quality,
                                            badQualityTrim,
                                            badQualityLimit)
                    rec.sequence = seq
                    rec.quality = qual
                    seqHash = hash(rec.sequence)
                    if DEBUG:
                        print "title: %s" % rec.title
                        print "sequence: %s" % rec.sequence
                        if hasattr(rec,'quality'):
                            print "quality: %s" % rec.quality
                        print "seqHash: %s" % seqHash

                    if removeRepeats:
                        rec.title = '%s|%s' % (
                            rec.title.split('|')[0],
                            bcSeqHashes[curBC][seqHash])

                    if removeRepeats:
                        if curBC in barcodes:
                            if bcSeqHashes[curBC][seqHash]!=0:
                                print >>bcFileDict[curBC][0], rec
                                bcSeqHashes[curBC][seqHash]=0
                        else:
                            if bcSeqHashes[''][seqHash]!=0:
                                print >>bcFileDict[''][0], rec
                                bcSeqHashes[''][seqHash]=0
                    else:
                        print >>bcFileDict[curBC][0], rec
                        
                map( lambda x: x.close(), [ y[0] for y in bcFileDict.values() ] )
                trimFiles = [x[1] for x in bcFileDict.values()]
            else:
                trimFiles = [fn]

            #print trimFiles

            for compFn in trimFiles:   
            
                # calculate rest of output file name
                if not output:
                    fnParts = os.path.splitext(fn)
                    outFn = fnParts[0]
                else:
                    outFn = output

                if len(barcodes) > 1:
                    curBC=os.path.splitext(compFn)[1]
                    if curBC=='.':
                        curBC='.nobarcode'
                    outFn+=curBC

                if not output:
                    if maxHomo >0:
                        outFn+='.maxH%s' % maxHomo

                    if minLZWsize >0:
                        outFn+='.LZW%s' % minLZWsize

                    if minLZWratio >0:
                        outFn+='.LZWr%s' % minLZWratio

                    if removeRepeats:
                        outFn+='.noDups'

                    if regexRequire != None:
                        outFn+= '.'+regexRequireName

                    if regexProhibit != None:
                        outFn+= '.'+regexProhibitName

                    if rejectQuality != (None,None):
                        outFn += '.Q%s_%sPos' % rejectQuality

                    if fivePTrim > 0:
                        outFn+='.5PT%s'% fivePTrim

                    if threePTrim > 0:
                        outFn+='.3PT%s' % threePTrim

                    if badQualityLimit is not None:
                        outFn+='.BQL%s' % str(badQualityLimit)
                        if badQualityTrim[0] is not None:
                            outFn+='_5PQC%s' % str(badQualityTrim[0])
                        if badQualityTrim[1] is not None:
                            outFn+='_3PQC%s' % str(badQualityTrim[1])

                    if minTrimLength > 0:
                        outFn+='.minBQT%s' % str(minTrimLength)
                
                if fastqOutput:
                    outFn+= '.fastq'
                else:
                    outFn += '.fasta'

                # if there is nothing more to be done
                # and repeats have been removed copy
                # that and go to the next file
                if (minLZWsize == 0 and minLZWratio == 0 and maxHomo == 0 and
                    regexRequire == None and regexProhibit == None and
                    fivePTrim==0 and outFn==0 and rejectQuality == None ):
                    if outFn != fn:
                        if compFn != fn:
                            shutil.move(compFn,outFn)
                            filesToClean.remove(compFn)
                            continue

                # filter line by line
                outFile = file(outFn, 'w')

                for rec in generalIterator(compFn):

                    # if all the trimming results in a null sequence, skip it
                    if len(rec) == 0:
                        continue
                    
                    # filters
                    if maxHomo > 0:
                        if rec.lenHomo() > maxHomo:
                            continue
                    if regexRequire != None:
                        if applyRegexSequence and regexRequire.search(rec.sequence) == None:
                            continue
                        if applyRegexTitle and regexRequire.search(rec.title) == None:
                            continue
                    if regexProhibit != None:
                        if applyRegexSequence and regexProhibit.search(rec.sequence) != None:
                            continue
                        if applyRegexTitle and regexRequire.search(rec.title) != None:
                            continue
                    if minLZWsize>0:
                        if rec.LZWsize() < minLZWsize:
                            continue
                    if minLZWratio >0:
                        if rec.LZWratio() < minLZWratio:
                            continue
                    if rejectQuality != (None,None):
                        if forceQuality != None:
                            if hasattr(rec,'quality') and rejectQuality[1] >= forceQuality:
                                continue
                        badCount = [q <= rejectQuality[1] for q in rec.quality ].count(True)
                        if badCount >= rejectQuality[0]:
                            continue
                    if minTrimLength > 0:
                        if len(rec.sequence) < minTrimLength:
                            continue

                    if fastqOutput:
                        if isinstance(rec,fastq.Record):    
                            print >> outFile, rec

                        else:
                            print >> outFile, fastq.Record(title=rec.title,
                                                           sequence=rec.sequence,
                                                           quality=[forceQuality] * len(rec.sequence))
                    else:
                        if isinstance(rec,fastq.Record):
                            print >> outFile, rec.fasta()
                        else:
                            print >> outFile, rec
                        

                outFile.close()

        finally:
            for fn in filesToClean:
                try:
                    os.unlink(fn)
                except:
                    pass
            


def removeDuplicates(inFilePaths,**kwargs):
    """Remove duplicate sequences over the domain of one or more fasta files.
    Add observed count information to the end of the record title '|<ct>',
    where <ct> is the number of times the sequence was observed in all the files.
    Output file(s), same number as inFilePaths,have '.noDups' added before the
    path extension.  If count information is present in the input files
    it is retained.

    if debug=True is passed in as a kwarg, terse progress information is printed
    on stderr.
    
    """

    if 'debug' in kwargs and kwargs['debug']:
        DEBUG=True
    else:
        DEBUG=False
    
    inFilePaths = list(getIterable(inFilePaths))
    seqHashes={}
    inReadTotal=0
    # get counts
    for f in inFilePaths:
        for rec in FastaIterator(f):
            seqHash = hash(rec.sequence)
            ct=rec.countFromTitle()
            if ct == None:
                ct =1
            inReadTotal+=1
            try:
                seqHashes[seqHash] += ct
            except KeyError:
                seqHashes[seqHash] = ct
        if DEBUG:
            print >> sys.stderr,(
                "removeDuplicates: %s cataloged, len(seqHashes): %s"
                %(f,len(seqHashes)))

    outStotal=0
    outRtotal=0
    for f in inFilePaths:
        inRoot, inExt = os.path.splitext(f)
        outName = inRoot+'.noDups'+inExt
        outFile=file(outName,'w')
        rCt=0
        sCt=0
        for rec in FastaIterator(f):
            seqHash = hash(rec.sequence)
            ct=seqHashes[seqHash]
            if ct is 0:
                # this seq already output
                continue
            if rec.countFromTitle() is not None:
                rec.title = '|'.join(rec.title.split('|')[:-1])
            rec.title = '%s|%s' %(rec.title,ct)
            seqHashes[seqHash]=0
            sCt+=1
            rCt+=ct
            print >> outFile, rec
        outFile.close()
        outRtotal+=rCt
        outStotal+=sCt
        if DEBUG:
            print >> sys.stderr,(
                "removeDuplicates: %s finished - %s reads - %s sequences"
                %(outName,rCt,sCt,))
    
    if DEBUG:
        print >> sys.stderr,("removeDuplicates: input  - %s reads - %s sequences"
                %(inReadTotal,len(seqHashes)))
        print >> sys.stderr,("removeDuplicates: output - %s reads - %s sequences"
                %(outRtotal,outStotal))
        print >> sys.stderr,("removeDuplicates: finished")

def screenHSPs(fastaFile,outFileName,MBfile,shortyLength=0,
               cushion=0,debug=False,excludeTitles=[]):
    """replace sequences that are query HSPs in megablast -D3 (or blast -m8)
    output file.

    Stretchs of bases =< shorty length are changed to XX...XX 
    
    Records are not output if:
    1) they have fewer ACGTU bases than shortyLength or
    2) their titles are in excludeTitles.

    output file will be opened for writing and then closed, on return.
    
    """

    import blastNoSQL

    filesToClose=[]

    if type(fastaFile) in StringTypes:
        fastaFile = file(fastaFile)
        filesToClose.append(fastaFile)
    if type(MBfile) in StringTypes:
        MBfile = file(MBfile)
        filesToClose.append(MBfile)
    outFile = file(outFileName,'w')
    filesToClose.append(outFile)
    
    hsps={}
    # read in HSPs
    i=0
    for hsp in blastNoSQL.m8generator(MBfile):
        hsp['query']=str(hsp['query'])
        if i % 10000 ==0 and debug: 
            print i
        i+=1
        if hsp['query'] not in hsps:
            hsps[hsp['query']]={hsp['q_start']: hsp['q_end']}

        else:
            if hsp['q_start'] in hsps[hsp['query']]:
                if hsp['q_end'] > hsps[hsp['query']][hsp['q_start']]:
                    hsps[hsp['query']][hsp['q_start']]= hsp['q_end']
            else:
                starts = hsps[hsp['query']].keys()
                starts.sort()
                for start in starts:
                    if hsp['q_start'] < start and  hsp['q_end'] > start:
                        end = hsps[hsp['query']][start]
                        del hsps[hsp['query']][start]
                        hsps[hsp['query']][hsp['q_start']]=end
                        break
                    elif hsp['q_start'] > start and  hsp['q_end'] <= hsps[hsp['query']][start]:
                        # do nothing
                        break
                    elif hsp['q_start'] > start and  hsp['q_end'] > hsps[hsp['query']][start]:
                        hsps[hsp['query']][ hsp['q_start']] = hsp['q_end']
                        break

    for rec in FastaIterator(fastaFile):
        if rec.title in excludeTitles:
            continue
        if  rec.title.split()[0] in hsps:
            newSeq = list(rec.sequence)
            rLen=len(newSeq)

            for start,end in hsps[rec.title.split()[0]].items():
                start -= cushion
                if start < 1:
                    start = 1
                end += cushion
                if end > rLen:
                    end = rLen
                newSeq[start-1:end] = 'X'*(end-start+1)
                if debug:
                    print start, end

##             if cushion == 0:
##                 for start,end in hsps[rec.title.split()[0]].items():
##                     newSeq[start-1:end] = 'X'*(end-start+1)
##                     print start, end
##             else:
##                 bounds = hsps[rec.title.split()[0]].items()
##                 bounds.sort()
##                 nonOL = [list(bounds[0])]
##                 for start,end in bounds[1:]:
##                     if start == nonOL[-1][0]:
##                         nonOL[-1][1] = end
##                     elif start <= nonOL[-1][0]:
##                         if end > nonOL[-1][0]:
##                             nonOL[-1][1] = end
##                     elif start - nonOL[-1][1] < cushion:
##                         nonOL[-1][1] -= cushion/2
##                         nonOL.append([start + cushion/2 + cushion%2,end])
##                     else:
##                         nonOL.append([start,end])
##                 for start,end in nonOL:
##                     newSeq[start-1:end] = 'X'*(end-start+1)
##                     print start, end
            rec.sequence = ''.join(newSeq)
        rec.screenShorties(shortyLength)
        if shortyLength == 0 or max(rec.matchLengths('[AGCTU]+')) >= shortyLength:
            print >> outFile, str(rec)    
    for f in filesToClose:
        f.close()    
    

            
screenMegaBLASThits = screenHSPs

def screenXmlBlastHits(fastaFile,outFileName,XMLfile,shortyLength=0,
               cushion=0,screenChar='X',
               debug=False,excludeTitles=None):
    """replace query sequences in fastaFile with XXX...XXX (or other
    screenChar) that are identical in HSPs in XML BLAST output, -m7,
    or megablast (-D2 -m7).

    Stretchs of bases =< shorty length are also changed to XX...XX 

    Screened fasta records are output to outFileName.
    
    Records are not output at all if:
    1) they have fewer ACGTU bases than shortyLength or
    2) their titles are in excludeTitles (which can be a set or a list).

    The output file will be deleted,opened for writing and then closed, on return.
    """

    import blastNoSQL

    filesToClose=[]

    if excludeTitles is None:
        excludeTitles = set()

    if type(fastaFile) in StringTypes:
        fastaFile = file(fastaFile)
        filesToClose.append(fastaFile)
    if type(XMLfile) in StringTypes:
        XMLfile = file(XMLfile)
        filesToClose.append(XMLfile)
    outFile = file(outFileName,'w')
    filesToClose.append(outFile)
    
    identRegions={}
    
    # read in Hits
    i=0
    for hit in blastNoSQL.ncbiXmlHitGenerator(XMLfile):
        if i % 10000 ==0 and debug: 
            print i, 'hits'
        i+=1
        qryRegions=blastNoSQL.ncbiHitIdenticalRegions(hit)['query']
        if cushion > 0:
            for i,r in enumerate(qryRegions):
                ## if r[0]==0:
                ##     qr1=r[0]
                ## else:
                ##     qr1=r[0]+cushion
                ## if r[1]==0:
                ##     print 'bad things'
                ##     qr2=r[1]
                ## else:
                ##     qr2=r[1]-cushion
                ## qryRegions[i]=(qr1,qr2)
                qryRegions[i]=(r[0]+cushion,r[1]-cushion)
        elif cushion < 0:
            for i,r in enumerate(qryRegions):
                s,e = r
                if s+cushion <0:
                    s=0
                else:
                    s=s+cushion
                qryRegions[i]=[s,e]
        myQuery=str(hit['query-def'])
        try:
            identRegions[myQuery].extend(qryRegions)
        except KeyError:
            identRegions[myQuery]=qryRegions

    #collapse overlapping ident regions
    for q in identRegions.iterkeys():
        rSort = sorted(identRegions[q])
        identRegions[q]=[]
        for s,e in rSort:
            if s>=e:
                continue
            try:
                rNow=identRegions[q][-1]
            except IndexError:
                identRegions[q]=[[s,e]]
            else:
                if s <= rNow[-1]:
                    identRegions[q][-1][-1]=e
                else:
                    identRegions[q].append([s,e])
    
    #now let's go over the input-output loop
    for rec in FastaIterator(fastaFile):
        if rec.title in excludeTitles:
            continue
        # do the masking screening
        if  rec.title in identRegions:
            newSeq = list(rec.sequence)
            rLen=len(newSeq)

            for start,end in identRegions[rec.title]:
                #check for if indent regions reach ends, screen all if they do
                #check for if the indent region is near the end, if so protect the minimum 2c+1
                if (start-cushion)<(2*cushion+1):
                    if (start-cushion)==0:
                        start=0
                    else:
                        start=2*cushion+1
                if (rLen-(end+1+cushion))<cushion:
                    if (end+cushion+1)==rLen:
                        end=end+cushion
                    else:
                        end=rLen-(2*cushion+2)

                
                newSeq[start:end+1] = 'X'*((end+1)-start)
                if debug:
                    print rec.title,start, end
                    print newSeq

            rec.sequence = ''.join(newSeq)
        
        # remove trivally short sequences
        rec.screenShorties(shortyLength)
        if shortyLength == 0 or max(rec.matchLengths('[AGCTU]+')) >= shortyLength:
            print >> outFile, str(rec)    
    #cleanup
    for f in filesToClose:
        f.close()    
    

            
def titleReSplit(inFile,matchFile,nomatchFile,patterns):
    """Partition a fasta file based on matching title to pattern(s).
    See also fastaPartition.
    """
    def cb (rec):
        return rec.titleMatch(patterns)
    
    return fastaPartition(inFile,matchFile,nomatchFile,cb)

def splitShorties(inFile,outFile,minDRNAStretch=50,screenAmbigAndShort=False):
    """Partition a fasta file based on the length of the
    longest stretch of valid bases in each record.  Records
    with fewer than  minDRNAStretch bases in their longest
    tract of valid bases are dropped. Others are place in
    outFile.
    
    See also fastaPartition.
    """    
    def cb (rec):
        if screenAmbigAndShort:
            s = screenRe(rec.sequence,notDRNAX)
            s = screenRe(s,DRNAbasesPat,maxLen=minDRNAStretch-1)
            rec.sequence=s
            
        if max(rec.matchLengths(DRNAbasesPat)) < minDRNAStretch:
            return False
        else:
            return True
        
    return fastaPartition(inFile,outFile,None,cb)


def dropLowLZW(inFile,outFile,minLZW):
    """Drop records with LZW sizes(or raitos if minLZW is a
    float) below minLZW.
    """
    if type(minLZW) in (types.IntType,types.LongType):
        cb = lambda r: r.LZWsize() >= minLZW
    elif type(minLZW) == types.FloatType:
        cb = lambda r: r.LZWratio() >= minLZW

    return fastaPartition(inFile,outFile,None,cb)

def dropDuplicates(inFile,outFile):
    """Count duplicate sequences in fasta file. Output
    distinct sequences.  Record titles are modified with |# . 
    returns number of unique sequences.
    """
    seqs={}
    for rec in FastaIterator(inFile):
        if rec not in seqs:
            seqs[rec]=1
        else:
            seqs[rec]+=1

    for rec,count in seqs.items():
        rec.title+='|%s'%count
        outFile.write(str(rec))
        outFile.write('\n')

    return len(seqs)            

def dropShorties(fileNames,minDRNAStretch=50,screenAmbigAndShort=False):
    """Drop Records with fewer than  minDRNAStretch bases in
    their longest tract of valid bases are dropped, per splitShorties.
    But the operation is done 'inplace'.  File Names is one or more
    paths of file to operate on.  Those file should be closed when
    this function is called.
    """

    if type(fileNames) in StringTypes:
        fileNames = [fileNames]
    for fileName in fileNames:
        tmpFile,tmpPath = mystemp(
            dir=os.path.split(fileName)[0])
        inFile = file(fileName)
        splitShorties(inFile,tmpFile,
                      minDRNAStretch=minDRNAStretch,
                      screenAmbigAndShort=screenAmbigAndShort)
        os.unlink(fileName)
        os.link(tmpPath,fileName)
        os.unlink(tmpPath)

def storeShorties (iFile,shortyLength=50, oFile='shortyRecords.fasta'):
    '''stores shorties from iFile into oFile in fasta format
    '''
    oFile=file(oFile,'a')
    iFile=file(iFile)
    for rec in FastaIterator(iFile):
        for s in rec.shortySubrecords(shortyLength=shortyLength):
            print >> oFile, s
    oFile.close()

            
def getLongestRecord(inFile):
    """Return the longest record (total valid D/RNA
    bases) from a file object.   
    """
    baseCount = 0
    longestRec = Record()
    for rec in FastaIterator(inFile):
        #print rec.DRNABaseCount()
        if len(rec) < baseCount:
            continue
        else:
            recCount = rec.DRNABaseCount()
            if recCount > baseCount:
                baseCount=recCount
                longestRec = rec
    return longestRec

def extractLongestRecord(inFileNames,leaveInFile=False):
    """Given one or more fasta file paths, find the longest
    record in the set.  remove it from the file it is found in
    and return the record.
    """

    if type(inFileNames) not in (ListType,):
        if type(inFileNames) in StringTypes:
            inFileNames = [inFileNames]
        else:
            inFiles = list(inFileNames)
         
    baseCount = 0
    longestRec = Record()
    recFileName = None
     
    for inFileName in inFileNames:
        inFile = file(inFileName)
        for rec in FastaIterator(inFile):
            if len(rec) < baseCount:
                continue
            else:
                recCount = rec.DRNABaseCount()
                if recCount > baseCount:
                    baseCount=recCount
                    longestRec = rec
                    recFileName = inFileName

    if recFileName != None and not leaveInFile:
        tmpFile,tmpPath = mystemp(
            dir=os.path.split(recFileName)[0])
        recFile = file(recFileName)
        for rec in FastaIterator(file(recFileName)):
            if rec.title != longestRec.title:
                tmpFile.write(str(rec))
                tmpFile.write('\n')
        recFile.close()
        os.unlink(recFileName)
        os.link(tmpPath,recFileName)
        os.unlink(tmpPath)

    return longestRec
       

def splitFasta(filename,splitCt=None,splitSize=None,
               tmpDir=None,dropTitles=[],nameGenerator=None):
    """Split fasta file into splitCt temporary files.  Return
    a list of paths to the files.

    If splitSize is specifed and splitCt is not, the file is split in
    to files with number of records equal to or less than splitSize (any
    remaining records are put in seperate file).

    A generator object can be specified
    which allows the user to specity custom paths for the slices.
    The generator must provide a file open for writing, and the path
    to the file.

    This could be as simple as:

    def mySpecialNames():
        choices = ['my','little','pony']
        for n in choices:
            yield (file(n,'w'),n)

    Unless tmpDir is set, new files are created in the current working
    directory.
            
    """


    if splitCt == None and splitSize == None:
        raise ValueError, "Either splitCt or splitSize mut be an integer" 
    

    rv=[]

    if tmpDir == None:
        tmpDir = os.getcwd()

    if nameGenerator == None:
        def nameFcn():
            while True:
                yield mystemp(suffix='.fasta',dir=tmpDir)
        nameGenerator = nameFcn()

    ff = file(filename)

    if splitCt != None:
        #make files
        LengthsAndFiles = []
        for n in range(splitCt):
            LengthsAndFiles.append([0,nameGenerator.next()])



        for rec in generalIterator(ff):
            if rec.title in dropTitles:
                continue
            LengthsAndFiles.sort()
            LengthsAndFiles[0][1][0].write(str(rec))
            LengthsAndFiles[0][1][0].write('\n\n')
            LengthsAndFiles[0][0] += len(rec.sequence)


        for length,tft in LengthsAndFiles:
            fObject, fPath = tft
            if length == 0:
                os.unlink(fPath)
            else:
                fObject.close()
            rv.append(fPath)
            
    elif splitSize != None:
        i=0
        oFile=None
        for rec in generalIterator(ff):
            if oFile == None:
                oFile,fName = nameGenerator.next()
                rv.append( fName)
            oFile.write(str(rec))
            oFile.write('\n\n')
            i+=1

            if i >= splitSize:
                i=0
                oFile.close()
                oFile=None
            
    else:
        raise ValueError, "Either splitCt or splitSize mut be an integer" 
    
    return rv


def tmpFile( fastaRecords, tmpDir=None ):
    """Given a fasta.Record or list of fasta.Records, creates
    a fasta file in the specified tmp directory containing said
    sequence(s), and returns the filename."""

    if not isinstance( fastaRecords, ( tuple, list ) ):
        fastaRecords = list(fastaRecords)
        
    if tmpDir == None:
        tmpDir = os.getcwd()
        
    (of, name) = mystemp(suffix='.fasta', dir=tmpDir)
    for record in fastaRecords:
        of.write( str(record) + "\n" )
    of.close()
    return name


def removeRecordsFromFile(fileName,recTitles,tmpDir=None,outFileName=None):
    """Given a path to a fasta file and a title (or titles) of fasta record(s),
    all matching fasta records for each title will be removed from the file.
    Leading '>' should should be removed (automatically stripping it causes
    ambiguious cases).  Trailing whitespace is not significant.
    Return value is the number (int) of records removed from the file.

    If outFileName is None, the original file will be overwritten.
    """
    import shutil
    if not os.access(fileName,os.W_OK):
        raise IOError, "no write access for fasta file: %s" % fileName

    if outFileName==None:
        outFileName = fileName

    inFile = file(fileName)

    (of, ofName) = mystemp(suffix='.fasta', dir=tmpDir)

    if type(recTitles) in StringTypes:
        recTitles = [recTitles]

    # remove training whitspace
    recTitles = [t.rstrip() for t in recTitles]


    

    BADREC = False
    rmCount=0
    for l in inFile:
        if l.startswith('>'):
            if l[1:].rstrip() in recTitles:
                BADREC=True
                rmCount += 1
            else:
                BADREC=False

        if not BADREC:
            of.write(l)
        
    of.flush()
    inFile.close()
    if outFileName==fileName:
        os.unlink(fileName)
    shutil.copyfile(ofName,outFileName)
    os.unlink(ofName)
    return rmCount
    

    
class IterativeScreen (object):
    """framework for iteratively screening through a fasta file.
    """

    def __init__ (self,inPath=None,MBparams="-W24 -D2 -m7 -E10 -G24 -FF",
                  maxRetries=1,cushion=0, shortyLength=35, mbTimeout=1000,
                  debug=False):
        """
        """
        self.MBparams = MBparams

        self.maxRetries=maxRetries
        self.cushion=cushion
        self.shortyLength=shortyLength
        self.mbTimeout=mbTimeout
        self.debug=debug


        self.longestRemaining=Record()
        self.retryCt=None
        self.lastMBstatus=None
        self.lastMBout=None
        self.lastMBerror=None
        self.minDRNAStretch=shortyLength

        if inPath != None:
            self.setup(inPath)

    def setup(self,inPath):
        """inFile is a Path.
        returns self.longestRemaining.DRNABaseCount() on success.
        """
        self.inPath = inPath
        localFile,self.localPath=mystemp(suffix='.fasta',dir='/tmp/md/')
        s,o = commands.getstatusoutput("cp %s %s" %(inPath,self.localPath))
        if s==0:
            self.dropAndStoreShorties()
            self.longestRemaining=getLongestRecord(file(self.localPath))
            return self.longestRemaining.DRNABaseCount()
        else:
            raise RuntimeError, "setup copy failed -  cp said: %s" %o



    def screen(self,screenFastaPath):
        """screen with MegaBLAST
        The longest record (largest DRNABaseCount) it  in stored as
        self.longestRemaining.
        Returns longestRemaining.DRNABaseCount() on success
        """

        self.retryCt=None
        dbTitles = fastaTitles(screenFastaPath)

        mbArgs = ['megablast'] + self.MBparams.split()
        mbArgs += [ '-d', screenFastaPath,
                    '-f', '-R', '-i', self.localPath,
                    '-o',self.localPath+'.MBout']

        if self.debug:
            print "running: '%s'" % ( ' '.join(mbArgs),)

        retryCt = 0
        self.lastMBstatus=None
        self.lastMBout=None
        self.lastMBerror=None
        while retryCt<self.maxRetries:
            mbProc=subprocess.Popen(
                mbArgs,
                stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            t0 = time.time()
            while mbProc.poll() == None:
                if time.time()-t0 > self.mbTimeout:
                    mbProc.poll()
                    break
                time.sleep(1)
            if mbProc.poll() != None:
                break
            os.system('kill %s' % mbProc.pid)
            time.sleep(1)
            mbProc.poll()
            retryCt+=1

        self.retryCt=retryCt
        self.lastMBstatus=mbProc.returncode
        self.lastMBout=mbProc.stdout.read()
        self.lastMBerror=mbProc.stderr.read()

        if self.lastMBstatus != 0:
            raise RuntimeError, (
                "---MEGABLAST FAILED---\nMEGABLAST command:%s\n"
                "End of stdout: %s\nEnd of stderr: %s\nExit Status: %s"
                % (' '.join(mbArgs),self.lastMBout[:-200],
                   self.lastMBerror[:-200],self.lastMBstatus))
        
        screenOutFile,screenOutPath = mystemp()
        screenXmlBlastHits(self.localPath,screenOutPath,self.localPath+'.MBout',
                           shortyLength=self.shortyLength,cushion=self.cushion,
                           excludeTitles=dbTitles)
        s,o = commands.getstatusoutput("mv %s %s" %(screenOutPath,
                                                    self.localPath))
        if s==0:
            self.longestRemaining=getLongestRecord(file(self.localPath))
            return self.longestRemaining.DRNABaseCount()
        else:
            raise RuntimeError, (
                "file replacment after screening failed-  mv said: %s" %o)
     

    def dropShorties(self):
        """
        """
        dropShorties(self.localPath,minDRNAStretch=self.shortyLength)
        self.longestRemaining=getLongestRecord(file(self.localPath))
        
    def storeShorties (self):
        '''
        '''
        name=str(self.localPath).split('/')
        name[-1]='shorty%s'%name[-1]
        shortyFileName='/'.join(name)
        storeShorties(self.localPath,shortyLength=self.shortyLength,oFile=name[-1])#shortyFileName)

    def dropAndStoreShorties(self):
        '''
        '''
        self.storeShorties()
        self.dropShorties()
    
    def cleanUp(self):
        os.unlink(self.localPath)
        

class RandomFastaPool (object):
    """Returns an object that acts as an iterator
    that returns a new random fasta record.  The
    title of the records is a sequential integer,
    starting at 0.

    records are also accessible via indexing of
    the object.    
    """

    def __init__(self,minLen,maxLen=None,**kwds):
        """    
        Arguments:
        - `self`:
        - `minLen`:
        - `maxLen`: optional used to specify a range of lengths
        - `**kwds`: keyword arguments passed to randomSequence
        """
        self.minLen=minLen
        self.maxLen=maxLen
        self.sequenceList=[]
        self.sequenceSet=set()
        self.rsKwds=kwds

    
    def next(self):
        """
        """
        if self.maxLen != None:
            if self.maxLen>self.minLen:
                l=random.randint(self.minLen,self.maxLen)
            else:
                l=self.minLen
        else:
            l=self.minLen

        s=None
        while s in self.sequenceSet:
            s=randomSequence(l,**self.rsKwds)
        self.sequenceSet.add(s)
        self.sequenceList.append(s)
        return Record(sequence=s,title=str(len(self.sequenceList)))


    def __getitem__(self,i):
        """
        """
        if i <0 or i>=len(self.sequenceList):
            raise IndexError
        return Record(sequence=self.sequenceList[i],
                      title=str(i))

    def __iter__(self):
        return self

class CattedFasta(object):
    """A class that represents a collection of fasta records, and
    can be printed out as a single record with a linker separating
    the sequnces.
    """
    def __init__ (self,records,linker='GGGGGGGGGG'):
        """records are be provided (as a python sequence-like
        object) here.
        """
        self.linker=linker
        self.records=tuple(records)
        self.rLengths=[len(r) for r in self.records]


    def catP2oldP(self,catP):
        """return the record and position in the original
        un-catted space, given the catP (catted position) in
        origin = 1 coordinates.
        """

        x=catP

        for i, l in enumerate(self.rLengths):
            x-=(l+len(self.linker))
            if x<0:
                break

        ## print i
        ## print x
        ## print x+self.rLengths[i]
        ## print catP

        if i==len(self.records)-1:
            rec=self.records[i]
            pos=catP - sum((self.rLengths[:i]))-(len(self.linker)*(i))
        else:
            rec=self.records[i+1]
            pos=catP - sum((self.rLengths[:i]))-(len(self.linker)*i)
        if pos >= len(rec):
            raise IndexError, 'catP falls in linker region'
        return rec,pos


    

    def fasta(self,title='catted fasta records'):
        """
        """
        return Record(
            title=title,
            sequence=self.linker.join((r.sequence for r in self.records)))
                      
                
        
            
        
        
    




def randomRecordIterator(filename, numResults=None, raw=False, visitedFunction=None,averageRecSize=10**4):
    """Returns a generator of `numResults` records in `file`. Supports
    fasta and fastq formats.

    Arguments:
    - `filename`: fasta or fastq file
    - `numResults`: number of random items to return. If None will
    iterate over all the records in the file.
    - `raw`: if True will return the raw record as it appears in the
    file instead of an instance of fasta.Record or fastq.Record.
    - `visitedFunction`: Function that accepts a record and
    returns True if the record should be returned in the
    iteration.
    """

    
    
    visitedFunction = visitedFunction
    visitedRanges = {}
    starts = []
    fastqRe = re.compile(r'^@([^\n]*)\n([^\+]*)\n\+\1\n([^@]*)\s+@([^\+]*)\n.*\n\+\4\n',
                         re.M)
    fastaRe = re.compile(r'^>([^\n]+)\n([^>]+)$', re.M)
    reObj = None
    def nextRandomPosition(maxPos):
        minPos = 0
        if len(starts) > 0 and starts[0] == 0:
            minPos = visitedRanges[starts[0]]
        if len(starts) > 1 and visitedRanges[starts[-1]] >= maxPos:
            maxPos = starts[-1]
        while True:
            pos = random.randint(minPos, maxPos)
            if visitedFunction is not None:
                return pos
            idx = bisect.bisect_right(starts, pos)
            if idx > 0 and visitedRanges[starts[idx-1]] > pos:
                continue
            if idx > 0 and idx < len(starts):
                lo = starts[idx-1]
                hi = visitedRanges[lo]
                if lo <= pos <= hi:
                    continue

            return pos

    def recInPosition(f, pos, recMark, retRaw):
        raw = ''
        startPos, endPos, oldStartPos = (long(pos),) * 3
        os.lseek(f, pos, os.SEEK_SET)
        chunkHalfSize = averageRecSize/2
        while True:
            # search to the left
            startPos -= chunkHalfSize
            if startPos < 0:
                startPos = 0
            os.lseek(f, startPos, os.SEEK_SET)
            chunk = os.read(f, oldStartPos-startPos)
            tpos = chunk.rfind(recMark)
            if tpos != -1:
                if tpos == 0:
                    startPos = 0
                    break
                if chunk[tpos-1] == '\n':
                    startPos += tpos
                    break
            oldStartPos = startPos
        while True:
            # search to the right
            os.lseek(f, endPos, os.SEEK_SET)
            chunk = os.read(f, chunkHalfSize)
            tpos = chunk.find('\n' + recMark)
            if len(chunk) == 0:
                break
            if tpos != -1:
                endPos += tpos
                break
            endPos += len(chunk)

        os.lseek(f, startPos, os.SEEK_SET)
        readSize = endPos - startPos
        
        if readSize > OS_READ_MAX:
            recordStr = ''
            nextSize = readSize - OS_READ_MAX
            while nextSize > 0:
                recordStr += os.read(f, OS_READ_MAX)
                nextSize -= OS_READ_MAX
        else:
            recordStr = os.read(f, readSize)
        if startPos > 0:
            startPos -= 2
        endPos += 2
        if not retRaw:
            match = reObj.search(recordStr)
            if recMark == '>':
                rec = Record(title=match.group(1),
                             sequence=match.group(2))
                rec.sequence = rec.dropNonLetters()
            else:
                raise Exception('Not implemented yet')
            return (rec, startPos, endPos)
        else:
            return (recordStr, startPos, endPos)
    
    f = open(filename)
    line = f.readline()
    while len(line) == 0:
        line = f.readline()
    if line.startswith('@'):
        recordMark = '@'
        reObj = fastqRe
        import fastq
    elif line.startswith('>'):
        recordMark = '>'
        reObj = fastaRe
    else:
        raise RuntimeError, ("Sequence iteration error: can't "
                             "determine type of file")
    f.close()
    
    f = os.open(filename, os.O_RDONLY)
    fsize = os.fstat(f).st_size
    curResults = 0
    while curResults < numResults:
        if curResults >= numResults:
            raise StopIteration()
        if ((len(starts) == 1) and 
            (visitedRanges[starts[0]] - starts[0] >= fsize)):
            raise StopIteration()

        curPos = nextRandomPosition(fsize)
        rec, sPos, ePos = recInPosition(f, curPos, recordMark, raw)

        if visitedFunction is None or visitedFunction(rec):
            idx = bisect.bisect_right(starts, sPos)
            lo = None
            if idx > 0:
                lo = starts[idx-1]
                hi = visitedRanges[lo]
                if lo <= sPos <= hi:
                    idx -= 1
                    visitedRanges[lo] = ePos
                else:
                    lo = sPos
                    starts.insert(idx, lo)
                    visitedRanges[lo] = ePos
            else:
                lo = sPos
                starts.insert(idx, sPos)
                visitedRanges[sPos] = ePos
            if idx < len(starts)-1:
                lo2 = starts[idx + 1]
                hi2 = visitedRanges[lo2]
                if lo2 <= ePos <= hi2:
                    visitedRanges[lo] = hi2
                    del visitedRanges[lo2]
                    del starts[idx + 1]


        if visitedFunction is None or not visitedFunction(rec):
            curResults += 1
            if curResults > numResults:
                break
            yield rec

    os.close(f)


def randomSubsequenceIterator(filename, length, numResults=None, visitedFunction=None):
    """Returns a generator of `numResults` random subsequences of
    length `length` from random fasta/fastq records in `filename`.
    """
    visitedRecords = dict()
    completedRecords = set([])
    isVisited = lambda x: x.title in completedRecords

    if type(visitedFunction) in (StringType,UnicodeType):
        visitedFunction = re.compile(visitedFunction)

    if type(visitedFunction) == ReType:
        searcher=visitedFunction.search
        visitedFunction=lambda x:searcher(x)

    f = open(filename)
    numRecords = fastaCount(f)
    f.close()
    
    curResultsCount = 0
    recordGenerator = randomRecordIterator(filename,
                                     visitedFunction=isVisited)
    while True:
        if curResultsCount == numResults:
            raise StopIteration()
        if len(completedRecords) == numRecords:
            raise StopIteration()
        rec = recordGenerator.next()
        if rec.title not in visitedRecords:
            visitedRecords[rec.title] = [[], dict()]

        starts, ranges = visitedRecords[rec.title]
        pos = random.randint(0, len(rec.sequence))
        if len(rec.sequence)-pos < length:
            continue
        idx = bisect.bisect_right(starts, pos)
        b = pos-length-1
        if b < 0:
            b = 0
        if idx > 0:
            if starts[idx-1] >= pos >= ranges[starts[idx-1]]:
                continue
            else:
                starts.insert(idx, b)
                ranges[starts[idx]] = pos+length
        else:
            starts.insert(idx, b)
            ranges[starts[idx]] = pos+length
        if len(starts) >= len(rec.sequence)/length:
            completedRecords.add(rec.title)
        subsequence = rec.sequence[pos:pos+length]
        if visitedFunction is None or (visitedFunction is not None and \
                                       not visitedFunction(subsequence)):
            curResultsCount += 1
            yield Record(title='%s_nt%d.%d' % (rec.title, pos,
                                               length),
                         sequence=subsequence)

