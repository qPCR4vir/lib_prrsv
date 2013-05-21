# compact.py
# Compact Sequence Representation
#
# by Kael Fischer
#
# The useful stuff here is in the numpy unsigned byte array classes:
# ByteArraySeq and ByteArraySeqBatch  and in the allied functions. 
#
# In python the overhead is such that you do not want to init 10^5 or more of any
# kind of object.  Hence the batch class idea.  The batch can be large but the
# sequences have to be the same length.
#
##############################################################################
# This used to be an interesting exercise in 3-bit packed sequences but
# that turned out to be a memory hog and very slow.  I wrote the sequence init
# code to left shift longs while reading the input sequence, this is a very bad
# idea vis-a-vis memory allocation speed.
#
# <historic stuff follows - and will be removed>
#
# 3-bit sequence representation
# stored in python long Objects
#
# the most significant bit is a 1 (placeholder)
# after that the most significant bit group is the 5' end
# 
# bit group is 'ABC'
# A  = ambiguity bit (call is 'low quality'), 0=high quality. 1=low quality.
# BC = the base call
#      00 -> A
#      01 -> C
#      10 -> G
#      11 -> T 
#
#############################################################################
__version__ = "$Id: compact.py,v 1.21 2011/12/05 20:01:31 kael Exp $"

import copy
import math
import string
import sys
from types import *

import numpy as N

from __init__ import *
import fasta
import utils
from utils import timeprofile

timer = timeprofile.timeprofile()

numByte = N.dtype('B')
numInt = N.dtype('i')

BASES = {
    'A':0,
    'a':0,
    'C':1,
    'c':1,
    'G':2,
    'g':2,
    'T':3,
    't':3,
    '.':4,
    'N':4,
    'n':4,
    0:'A',
    1:'C',
    2:'G',
    3:'T',
    4:'N'
    }

BASEMASK=7         # 111
AMBIGUOUS = 4      # 100
NOTAMBIGUOUS = 0   # 0
COMPLEMENT = 3     # 011 xor a base to the the complement 

BITSTR =    {
    'A':'00',
    'a':'00',
    'C':'01',
    'c':'01',
    'G':'10',
    'g':'10',
    'T':'11',
    't':'11'
    }

B0mask = string.atoi('00000011',2)
B1mask = string.atoi('00001100',2)
B2mask = string.atoi('00110000',2)
B3mask = string.atoi('11000000',2)
BaseMasks = (B0mask,B1mask,B2mask,B3mask)


def extBase(inByte,pos):
    """position can be 0,1,2,3 
    """
    shiftN = pos*2
    mask = string.atoi('00'*pos+'11',2)
    return BASES[(inByte >> shiftN) & mask]
    
    
def reverseSeqByte(b):
    """
    """

    if hasattr(b,"dtype"):
        rv=N.zeros_like(b)
    else:
        rv = 0


    pShift = {0:6,1:2,2:-2,3:-6}
    for p,m in enumerate(BaseMasks):
        if pShift[p] > 0:
            rv = rv | ((b & m) << pShift[p] )
        else:
            rv = rv | ((b & m) >> abs(pShift[p]) )
    return rv


def bases():
    choices = 'ACGT'
    for c in choices:
        yield c

bBytes = {}
for b1 in bases():
    for b2 in bases():
        for b3 in bases():
            for b4 in bases():
                seqKey = ''.join((b4,b3,b2,b1))
                seqNumber = ''.join((b1,b2,b3,b4))
                bBytes[seqKey] = string.atoi(''.join([BITSTR[x] for x in seqNumber]),2)


    

def seq2twoBitStr(s):
    """takes one of more DNA base letters and converts them to a
    string of 0's and 1's translating per the BITSTR dictionary.
    """

    if hasattr(s,'sequence'):
        s=s.sequence
    
    try:
        return ''.join([BITSTR[x] for x in s])
    except KeyError , msg:
        print "Non-DNA base in sequence:", msg

def seq2bytes(seq):
    i=range(0,len(seq),4)
    j=range(4,len(seq),4)
    if len(seq) % 4 == 0:
        j.append(len(seq))
    idx = zip(i,j)

    nibbles = [seq[s:e] for s,e in idx]
    # padd to a full byte
    if idx[-1][1] <=  (len(seq)-1):
        nibbles.append(seq[idx[-1][1]:]+'A'*(4-len(seq[idx[-1][1]:])))

    return [bBytes[n] for n in nibbles] 

def diffCount(xorBytes):
    try:
        i=int(xorBytes)
    except:
        i=None
        pass
    if type(i) in( IntType,LongType):
        return int((i & B0mask) > 0) + int((i & B1mask) > 0) +int((i & B2mask) > 0) + int((i & B3mask) > 0)
    elif type(xorBytes) == N.ndarray:
        return ( ((xorBytes & B0mask) != 0).astype(numInt) + ((xorBytes & B1mask) != 0).astype(numInt) +
                 ((xorBytes & B2mask) != 0).astype(numInt) + ((xorBytes & B3mask) != 0).astype(numInt)  )
    else:
        i = iter(xorBytes)
        return utils.flatten([diffCount(b) for b in xorBytes])

def oneZeroDiffTF(xorBytes):
    
    return (xorBytes == 0) | ((xorBytes & B0mask != 0) ^
                              (xorBytes & B1mask != 0) ^
                              (xorBytes & B2mask != 0) ^
                              (xorBytes & B3mask != 0)
                              )
    




def diffPos(xorByte):
    '''returns an ordered list of integer indices in the xorByte that
    correspond to base differences 
    '''
    return [n for n,mask in enumerate((B0mask,B1mask,B2mask,B3mask)) if xorByte & mask !=0 ]

def byte2seq (byte):
    return ''.join([extBase(byte,n) for n in (0,1,2,3)])


def appendOnBits(thing,count) :
    """shift thing left putting 1's in the right most positions after shifting.
    """
    while count > 0:
        thing = (thing << 1) | 1
        count -= 1
    return thing


def bitMap(thing):
    """Show the bitmap of the thing
    """
    rl=[]
    while thing:
        rl.append(str(int(thing & 1)))
        thing = thing >> 1

    rl.reverse()
    return ''.join(rl)


class ByteArraySeq  (object):
    """A single sequence packed in a 2 bit representation.
    Sequence is held in a numpy array of unsigned 1-byte integers.

    .seqArray holds one sequence per row, 4 bases per element,
    5' end of sequence corresponds to the 0th element in the row

    base packing is set up with 2 least significant bits being the 5'
    most base in the byte quartet. 

    Important Attributes:
    .seqArray - sequence data
    .length - length of each sequence in bases
    .lastByteMask - a mask that can be and'ed to
                    comparisons of these positions
                    to mask terminal undefined bases
 

    """

    def __init__ (self,seqBytes,length,title='ByteArraySeq'):
        if len(seqBytes.shape) > 1:
            raise ArgumentError, "seqBytes must be 1-dimensional"
        self.title=title
        self.seqArray=seqBytes
        self.length=length
        self.lastByteMask=255

    def __str__ (self):
        return ''.join([byte2seq(b) for b in self.seqArray])[:self.length]

    def fasta(self):
        import fasta
        return fasta.Record(sequence=str(self),title=self.title)

    def __getitem__(self,x):
        return ''.join([byte2seq(b) for b in self.seqArray])[x]

    def __len__(self):
        return self.length

    def complement(self):
        """Complement the sequences (in place).
        """
        self.seqArray = ~self.seqArray

    def reverse(self):
        """reverse the order of the bases of the sequences.
        """
        #flip rows
        try:
            self.seqArray=N.fliplr(self.seqArray)
        except ValueError:
            self.seqArray=N.flipud(self.seqArray)
        
        #flip base packing
        self.seqArray=reverseSeqByte(self.seqArray)
        
    def rc(self):
        """reverse and complement the batch (in place)
        """
        self.reverse()
        self.complement()

    def shift3p (self,count=1):
        """Shift sequence in 3-prime direction by count bases.
        A's will be introduced at the 5-prime end, and 3-prime
        bases will vanish.
        """

        # handle 1d case
        oShape = self.seqArray.shape
        if len(oShape) == 1:
            self.seqArray=self.seqArray.reshape(1,oShape[0])
            
        while count > 0:
            count -=1
            # get the new bases for the 5' end of each byte
            tempBase = self.seqArray & B3mask  # get 3' base
            tempBase[:,1:] = (tempBase[:,:-1]) >> 6 # shift to 5' bits and
                                                    # store in the next byte over
            tempBase[:,0] = 0                       # 'A' goes in the first pos
            self.seqArray = (self.seqArray << 2) + tempBase # push 3' bases out of the airlock
                                                            # add the stored 5' bases
        if len(oShape) == 1:
            self.seqArray=self.seqArray.reshape(oShape[0])
       

    def shift5p (self,count=1):
        """Shift sequence in 5-prime direction by count bases.
        A's will be introduced at the 3-prime end, and 5-prime
        bases will vanish.
        """

        # handle 1d case
        oShape = self.seqArray.shape
        if len(oShape) == 1:
            self.seqArray=self.seqArray.reshape(1,oShape[0])
            
        while count > 0:
            count -=1
            # get the new bases for the 3' end of each byte
            tempBase = self.seqArray & B0mask  # get 5' base
            tempBase[:,:-1] = (tempBase[:,1:]) << 6 # shift to 3' bits and
                                                    # store in the next byte over
            tempBase[:,-1] = 0                      # 'A' goes in the last pos
            self.seqArray = (self.seqArray >> 2) + tempBase # push 5' bases out of the airlock
                                                            # add the stored 3' bases
        if len(oShape) == 1:
            self.seqArray=self.seqArray.reshape(oShape[0])

    def ends(self,length):
        """Returns a new ByteArraySeqBatch instance where the sequence are the 3'
        terminal 'length' bases.
        """

        rv = copy.deepcopy(self)

        # trim extraneous 5' bytes
        extra5p = rv.length-length
        rv.seqArray = N.delete(rv.seqArray,range(extra5p/4),axis=1)

        # shift 5' some # of bases
        # so first base lines up with byte 0
        rv.shift5p(extra5p%4)

        # trim extraneous 3' bytes
        extra3p = (rv.seqArray.shape[1]*4)-extra5p%4-extra5p/4-length
        blankBytes = extra3p/4
        rv.seqArray = N.delete(rv.seqArray,range(length/4 + (length%4 !=0),rv.seqArray.shape[1]),axis=1)

        rv.length=length
        return rv
    

    def starts(self,length):
        """Returns a new ByteArraySeqBatch instance where the sequence are the 5'
        terminal 'length' bases.
        """
        rv = copy.deepcopy(self)

        # trim extraneous 3' bytes
        extra3p = rv.length-length
        rv.seqArray = N.delete(rv.seqArray,range(length/4 + (length%4 !=0),rv.seqArray.shape[1]),axis=1)

        # set extraneous 3' bases to 'A'
        retainBases = (length%4)
        if retainBases != 0:
            passthroughMask=0
            for i,m in enumerate(BaseMasks):
                if i >= retainBases:
                    break
                passthroughMask += m
            rv.seqArray[:,-1] =  rv.seqArray[:,-1] & passthroughMask

        rv.length=length
        return rv

    def nearMatchIdxTest(self,seqBatch,totalMm=0,firstQuadMm=0):
        """Returns the indices of sequences in seqBatch
        which match with mm or fewer mismatches.
        """

        if firstQuadMm > totalMm:
            raise ArgumentError , "firstQuadMm can not me greater that totalMm"  
        
        i=N.r_[0:seqBatch.seqArray.shape[0]]
        xorR = seqBatch.seqArray[:,1] ^ self.seqArray[-1]
        if firstQuadMm == 0:
            condArray = xorR!=0
        else: 
            condArray = diffCount(xorR) > firstQuadMm
        looseRows = N.extract(condArray,i)
        candidates = N.delete(seqBatch.seqArray,looseRows,axis=0)
        i=N.delete(i,looseRows)
        #print candidates.shape
        xorR=candidates ^ self.seqArray
                
        return N.extract(N.sum(diffCount(xorR),axis=1)<=totalMm,i)
        
    
    def nearMatchIdx(self,seqBatch,mm=0):
        """Returns the indices of sequences in seqBatch
        which match with mm or fewer mismatches.
        """

        i=N.r_[0:seqBatch.seqArray.shape[0]]
        xorR = seqBatch.seqArray ^ self.seqArray
                
        return N.extract(N.sum(diffCount(xorR),axis=1)<=mm,i)
        

class CompactBatch (object):
    def insertFasta(self,fastaPath,ignoreDifferentLengths=True,keepTitles=False):
        """Read a fasta file and initilize the batch sequences from
        the fasta file sequences.

        If the sequence length for this instance has not been
        specified, the length of the first sequence is used as
        the required length for all sequences in this batch.

        If ignoreDifferentLengths is True, fasta records with different
        lengths are silently skipped.  Otherwise warnings are printed
        for each such sequence.
        """

        if self.length == None:
            self.length=-1

        baddies = 0
        # check input
        for n,rec in enumerate(fasta.FastaIterator(fastaPath)):
           if self.length < 0:
               self.length=len(rec)
           else:
               if self.length!=len(rec):
                   baddies+=1
                   if not ignoreDifferentLengths:
                       print "Skipping %s, length %s" % (rec.title,len(rec))
                       continue
           

        # init zero array (all A's !!!)
        self.seqArray = N.zeros((n-baddies+1,self.length/4+int(self.length%4>0)),numByte)

        # insert sequences 
        baddies = 0
        for n,rec in enumerate(fasta.FastaIterator(fastaPath)):
            if keepTitles:
                self.titles.append(rec.title)
            if self.length!=len(rec):
                baddies+=1
                continue
            for j,byte in enumerate(seq2bytes(rec)):
                self.seqArray[n-baddies,j]=byte

    def __init__ (self,fastaPath=None,length=None,keepTitles=False):
        """make a new batch.
        if fasta path is given sequences will be read in
        from that path.

       

        
        """

        self.lastByteMask = 255
        self.seqArray = N.array([],numByte)
        self.length=length
        self.titles=[]

        if fastaPath!= None:
            self.insertFasta(fastaPath,keepTitles=keepTitles)

    def __getitem__(self,n):
        """Returns a ByteArraySeq from the batch.
        """
        if len(self.titles)>0:
            return ByteArraySeq(self.seqArray[n],self.length,self.titles[n])
        else:
            return ByteArraySeq(self.seqArray[n],self.length,str(n))

    def __len__ (self):
        return self.seqArray.shape[0]


class ByteArraySeqBatch (CompactBatch,ByteArraySeq):
    """A set of sequences, all of the same length, packed in a 2 bit representation
    """
    

    def _set3pMask(self):
        """
        """
        if self.lastByteMask ==None:
            return self.lastByteMask
        
        myMasks = list(BaseMasks[:])
        myMasks.reverse()
        
        self.lastByteMask = 255
        for n,m in enumerate(myMasks):
            if n > self.length%4-1:
                break
            self.lastByteMask = self.lastByteMask ^ m
        return bitMap(self.lastByteMask)


    def __str__(self):
        return '\n'.join([ "Packed Sequence Batch: %s %smers" % (len(self),self.length),
                           "first seq:\t%s"% str(self[0])[:self.length],
                           ".\n.\n.\n"
                           "last seq:\t%s"% self[-1]])        

    def findMarchJoins(self,marchOffset=11,mmAllowed=0):
        """Return a list of tuples corresponding to allowed joins.
        """
        pass
    
class JoinFinder (object):
    """For comparing two batches of sequences.  Called starts and ends b/c
    initial application was assembly of short reads, where one batch was
    the ends of all the reads and the starts were the start of every read.

    The two batches must have sequences of the same length.

    The methods all find starts and ends which are compatible, e.g.
    exact matches of near matches.

    Here, the "ends" are the ends of the reads, and the "starts" of the reads.
    Unfortunately if you are thinking about a "join" the "ends" are the 5'
    part of the join...  So that can be confusing.


    If your application is not assembly of a single set of reads, and the batches
    you want to compare have significantly different numbers of sequences, you
    should use the smaller set as 'ends'  
    
    """
    
    def __init__(self,ends,starts):
        self.e=ends
        self.s=starts
        self.joins={}


    def starts4end(self,endReadIdx,mm=0):
        """Return indices that correspond to starts that are compatible with
        the end(s) specified in endReadIdx, which can be a sequence, numpy.Array,
        or just and integer.  The closeness of the match requirement is contorled
        with mm which is the allowed number of mismatches in a pair comparison.

        The return value is a dictionary like d[endIdx] = [startIdx1, ..., startIdxn]
        The attribute .joins is also set to the return value.

        This is kinda slow and you should probably use slicingEndMatch, instead.
        
        note: in rpyc, if endReadIdx is a sequence it should be a tuple.
        """
        self.joins = {}

        if type(endReadIdx) in (IntType, LongType):
            endReadIdx = (endReadIdx,)

        for i in endReadIdx:
            starts = self.e[i].nearMatchIdx(self.s,mm=mm)
            if len(starts) > 0:
                self.joins[i]=starts

        return self.joins

    def starts4endTest(self,endReadIdx,mm=0):
        """note: in rpyc if endReadIdx is a sequence it should be a tuple
        """
        self.joins = {}

        if type(endReadIdx) in (IntType, LongType):
            endReadIdx = (endReadIdx,)

        for i in endReadIdx:
            starts = self.e[i].nearMatchIdxTest(self.s,mm)
            if len(starts) > 0:
                self.joins[i]=starts

        return self.joins
        

    def slicingEndMatch(self,eIdx1,eIdx2,mm=0):
        """Find compatible starts for multiple ends in a continuous set of sequences
        in 'ends' specified by eIdx1,eIdx2.
        Returns a 2-tuple of numpy arrays like (  array([end1, end2, ... ,end_n]),
                                                  array([start1, start2, ..., start_n]) )
        
        Hence a list of (end,start) index pairs is obtained from the above tuple, x,
        like so: zip(x[0],x[1])

        While the problem with starts4ends is speed, the problem with this method is memory.
        You should use fairly small batches.  When compating 
        """
        
        sumDiff=N.sum(diffCount(self.s.seqArray ^ N.array(N.vsplit(self.e.seqArray[eIdx1:eIdx2,:],eIdx2-eIdx1))),
                      axis=2)
        eArray,sArray = (sumDiff<=mm).nonzero()
        eArray += eIdx1
        return (eArray,sArray)



    


class MarchAssembler(CompactBatch):

##     def __init__ (self, joinsPath=None,readsPath=None):
        
##         self.insertFasta = ByteArraySeqBatch.insertFasta
##         pass

    



    def readJoinFile(self,joinPath):
        self.march=11
        self.cStarts={}
        self.cEnds={}
        for l in file(joinPath):
            eIdx,sIdx = [int(x) for x in l.strip().split()]
            if eIdx not in self.cStarts:
                self.cStarts[eIdx]=[]
            if sIdx not in self.cEnds:
                self.cEnds[sIdx]=[]
            
            self.cStarts[eIdx].append(sIdx)
            self.cEnds[sIdx].append(eIdx)
        self.used={}
        #self.joinTree={}
            
    def joinTree(self,startAt):

        self.lastS=startAt
        rv={}
        #print self.endsSeen
        if startAt not in self.endsSeen:
            self.endsSeen[startAt]=None
            if startAt in self.compat:
                subTrees = [self.joinTree(x) for x in self.compat[startAt] ]
                rv[startAt]=[x for x in subTrees if len(x) > 0]
            else:
                rv[startAt]=[]
        return rv


    def find5pJoined(self,startAt):
        """find the left most reads joined to the read with idx startAt. 
        """

        rv = []
        reads=[startAt]
        clean = False
        while len(reads) >0:
            newReads =[]
            clean = True
            for n,read in enumerate(reads):
                if read not in self.cEnds:
                    if read not in rv:
                        rv.append(read)
                else:
                    newReads.extend(self.cEnds[read])
                
                
            reads=newReads
        return rv
        
        
    def printTree(self,startAt,depth=0,outFile=sys.stdout):
        """
        """
        if depth == 0:
            outFile.write('\n\n')
            #self.endsSeen ={}
            
        self.lastS=startAt
        if startAt not in self.used:
            outFile.write( "%-10d%s%s\n" % (startAt,' '*(depth*11), self[startAt]))
            self.used[startAt]=None
            if startAt in self.cStarts:
                for  x in self.cStarts[startAt]:
                    self.printTree(x,depth=depth+1,outFile=outFile)


    def path2Fasta(self,indexList):
        """
        """
        pieces=[[self[indexList[0]][:self.march]]]
        for li,ci in enumerate(indexList[1:-1]):
            pieces.append(self[indexList[li]][self.march:])
            pieces.append(self[ci][self.length-self.march:])

        return pieces
                          
        



class TwoBitSequence (object):
    """2-bit sequence representation
    the most significant bit is a 1 (placeholder)
    most significant bit group is the 5' end

    bit groups are 2 bits
      00 -> A
      01 -> C
      10 -> G
      11 -> T
    
  
    """
    def __init__(self,sequence):
        """bases can be a list of Base objects, or a string of [ACGT]s.
        """
        self._seqBits=string.atol('1'+seq2twoBitStr(sequence),2)

    def bitMap(self):
        return bitMap(self._seqBits)

        
class CompactSequence (object):
    """3-bit sequence representation
    
    the most significant bit is a 1 (placeholder)
    most significant bit group is the 5' end
     
    bits are 'ABC'
    A = ambiguity bit (call is 'low quality'). 0=high quality. 1=low quality.
    BC = the base call
      00 -> A
      01 -> C
      10 -> G
      11 -> T

      When A=1, BC are not significant. 
    """

    def __init__(self,bases=[],qualityStr='',qualityCutoff=0):
        """bases can be a list of Base objects, or a string of [ACGT]s.
        """
        self._seqBits=1L

        if type(bases) == StringType:
            self.appendSeqStr(bases)
        else:
            self.appendBases(bases)


    def ambigMask(self):
        """Returns a bit mask of repeated 100 for the full length of the
        sequence.
        """
        rv=1L
        for n in range(len(self)):
            rv =  rv << 3 | AMBIGUOUS
        return rv

    def nonAmbigBits(self):
        """Returns a mask of 1's in non-ambiguous positions.
        """
        rv=1L
        for b in self:
            if b.ambiguity() == True:
                rv=rv << 3
            else:
                rv=rv << 3 |  BASEMASK
        return rv


    def placeHolderMask(self):
        """Returns a bit mask of 100000...000.
        """
        rv=1L
        return rv << self._bitCount()-1
    

    def unambiquious(self):
        """True if CompactSequence has no N's.
        """
        return self.placeHolderMask()==self.ambigMask() & self._seqBits
    
 
    def appendSeqStr(self,seqStr):
        """Add bases to the 3' end given a string.
        """

        for base in seqStr:
            self._seqBits = (self._seqBits << 3) | BASES[base] 

    def appendBases(self,bases):
        """Add bases to the 3' end given a list of Base objects.
        """

        if type(bases) not in (TupleType,ListType):
            bases = (bases,)
        while len(bases) > 0:
            base=bases.pop(0)
            if not isinstance(base,Base):
                raise ArgumentError, "% is not a base" % base
            self._seqBits=self._seqBits << 3 | base._bits


    def _bitCount(self):
        """Including the leftmost placeholder return how many bits _seqBits is.
        """
        return int(math.ceil(math.log((self._seqBits),2)))
    
    def __len__(self):
        return (self._bitCount()-1)/3

    def __cmp__(self,other):
        return cmp(self._seqBits,other._seqBits)
    
    def  __str__(self):
        """Normal string representation of sequence.
        """
        return ''.join([x.letter() for x in self])

    def __getitem__(self, key):
        """Single index return a base.
        Slices return a CompactSequence.
        """
        l = len(self)

        # slice returns a sequence
        if type(key) == SliceType:
            if key.start < 0:
                raise IndexError
            start = key.start
            if key.step == None:
                step = 1
            else:
                step=key.step
            if key.stop > l:
                stop=l
            else:
                stop=key.stop
            return CompactSequence([self[i] for i in range(start,stop,step)])

        if key > l-1 or key < -1*l:
            raise IndexError
        baseOffset = l-key-1
        return Base(self._seqBits>>(baseOffset*3) & BASEMASK)
        
    def __setitem__ (self,key,base):
        """set a position to a particular base (Base or [ACGT])
        """
        if type(base) == StringType:
            base=Base(BASES[base])  # wow what unfortunate syntax!

        if not isinstance(base,Base):
            raise ArgumentError, "% is not a base" % base
        
        mask = 0
        mask = appendOnBits(mask,len(self)*3 +1) # all 1 mask
        submask = BASEMASK << (3*(key)) # 1's in target location
        mask = mask ^ submask
        self._seqBits = self._seqBits & mask # target bits = 0
        newBits = base._bits << (3*(key))

        self._seqBits = self._seqBits |  newBits
        
    def __iter__(self):
        """bases 5'->3.
        """
        for n in range(len(self)):
            yield(self[n])

    def __hash__(self):
        return hash(self._seqBits)

    def find5(self,query,start=0):
        """search for query sequence starting from start in the 5'->3' direction.
        query can be a CompactSequence or a string.
        the sequence index is returned or -1 if not found (python style)
        """

        if not isinstance(query,CompactSequence):
            query = CompactSequence(query)
        
        passThough = 0
        passThough = appendOnBits(passThough,query._bitCount()) << (self._bitCount() - query._bitCount()-1-start*3)
        queryBits = query._seqBits << (self._bitCount() - query._bitCount() - start*3) & passThough
        
        pos=start
        while pos < len(self):
            #print "%35s" % bitMap(self._seqBits)
            #print "%35s" %  bitMap(passThough)
            #print "%35s" %  bitMap(queryBits)
            if self._seqBits & passThough ^ queryBits == 0:
                return pos
            pos += 1
            passThough = passThough >> 3
            queryBits = queryBits >> 3
            #print
        return -1
                
    def find3(self,query):
        """exercise for the reader. (not written yet)
        """
        pass

    def reverse(self):
        """Reverse sequence in place (returns None).
        """
        
        bases=list(self)
 
        self._seqBits=1L
        while len(bases) > 0:
            self._seqBits=self._seqBits << 3 | bases.pop(-1)._bits
        
    
    def complement(self):
        """Complement in place (returns None).
        """
        mask = 0
        for l in range(len(self)):
            mask = mask <<3 | COMPLEMENT

        self._seqBits= self._seqBits ^ mask
        

    def compareWith(self,other,ignoreN=True):
        """return long integer representing XOR'ed
        seq bits.  The sequences must be the same length.
        
        0L==identical

        if not identical, the place holder bit is retained.


        If ignoreN is true (default) ambiguous positions
        are set to 000.
        """

        if len(other) != len(self):
            return ValueError , "sequences must be the same length"
        
        rslt = self._seqBits ^ other._seqBits

        if rslt == 0L:
            return rslt

        if not ignoreN or (self.unambiquious() and other.unambiquious()):
            return rslt | self.placeHolderMask()
        else:
            rslt = rslt & self.nonAmbigBits() & other.nonAmbigBits()
            if rslt == 0L:
                return rslt
            else:
                return rslt | self.placeHolderMask()

    def countDiff(self,other,ignoreN=True):
        """return the number of sequence differences
        with other.
        """
        count=0

        diffBits = self.compareWith(other,ignoreN=ignoreN)
        if diffBits == 0L:
            return count
        
        mask=BASEMASK
        for n in range(len(self)):
            if diffBits & mask > 0:
                count+=1
            diffBits = diffBits >> 3
        return count
            

    def count(self,base):
        """Return the number of times a base occurs.
        """

        count=0
        
        if not isinstance(base,Base):
            base = Base(BASES[base])

        bBits=base._bits
        mask = BASEMASK
        myBits = self._seqBits
        for n in range(len(self)):
            if (myBits & mask) ^ bBits == 0:
                count += 1
            myBits = myBits >> 3
        return count
                

            

class Base (object):
    """A base.  CompactSequence.__iter__() returns these,
    and the CompactSequence constructor takes them.
    """

    def __init__ (self,bits,showN=True):
        self._bits = bits
        self.showN = showN


    def ambiguity(self):
        return self._bits & AMBIGUOUS >= AMBIGUOUS

    def letter(self):
        if self.showN and self.ambiguity():
            return 'N'
        else:
            return BASES[self._bits & 3]

    def __str__ (self):
        return self.letter()

    def __cmp__ (self, other):
        return cmp(self._bits, other._bits)
    
    
    
    
A = Base(0)
C = Base(1)
G = Base(2)
T = Base(3)
#N = Base(AMBIGUOUS)
