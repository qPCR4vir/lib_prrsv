#
# AOS Calculations
#
# Access to the ArrayOligoSelector compiled libraries.
#
# Smith Waterman uses the PAM47 matrix:
#              gap_open = -7
#              gap_extend = 0
#              match = 5
#              mismatch = -4
#
#
# Energy Calcultaions - comments from the C code (data.c):
#
#/* All energy parameters are in kcal/mol unit 
#   All energy parameters are in 37 degree C standard condition
#*/
#
#/*N. Sugimoto et.al.  NAR  1996, Vol 24, No22 4501-4505*/ [wc base pair energies]
#
#/*N. Peyret et.al. Biochemistry 1999, 38, 3468-3477
#  Biochemistry Vol37, No 26, 1998
#  NAR 1998, Vol 26, No.11
#  Biochemistry, Vol 37, No. 8, 1998
#  Biochemistry, Vol 36, No. 34, 1997 [internal mismatches]
#*/
#
# [loop params seem to come for Zuker web pages that are gone]
#/* [also] 
#   Peritz et. al. Biochemistry 1991 30, 6428-6436
#*/
#

__version__ = "$Revision: 1.10 $"
__author__ = "Kael Fischer"

from types import *
import ctypes
import os.path
import string
import re
#import sequence      # copied required items to this module

myDir = os.path.split(__file__)[0]
if myDir == '':
    myDir='.'
__aos = _aos=ctypes.cdll.LoadLibrary(os.path.join(myDir,"libaos.so.1.0.2"))
__aos.energy.restype=ctypes.c_double

def SW(s1,s2):
    """Return (alignment string, sw score).
    sequences are taken as read - no complementing or reversing.
    returns ('',-1) if there is no alignment.

    Score is an integer using the PAM47 matrix.
    """

    if type(s1) != StringType or type(s2) != StringType:
        raise AttributeError, "s1  and s2 must be strings"
        
    sizeDiff = len(s1)-len(s2)

    if sizeDiff > 0:
        s2+='x'*sizeDiff
    elif sizeDiff < 0:
        s1+='x'*abs(sizeDiff)
    s2=reverse(s2)
    sLen = len(s1)
    outBuffer =  ctypes.create_string_buffer(sLen*5)

    __aos.SW(ctypes.c_char_p(s1),
             ctypes.c_char_p(s2),
             ctypes.c_int(1),
             ctypes.c_int(sLen),
             outBuffer)

    outLines = outBuffer.value.split('\n')
    #print outLines
    if len(outLines) > 3:
        alignment = '\n'.join(outLines[0:2])
        score = int(outLines[2])
        return (alignment,score)
    else:
        return ('',-1)

def selfSW(s):
    """Return the self-alignment, i.e. the sequence to the
    reverse complement sequence.
    The better the alignment the more secondary structure.
    """
    return SW(s,reverseComplement(s))



def LZWsize(s):
    """Return size of LZW compresses string.
    A smaller size means lower complexity if you are
    comparing 2 sequences of the same length.
    """

    sLen=len(s)
    return __aos.LZWsize(ctypes.c_char_p(s),ctypes.c_int(sLen))


                        
def energy(s1,s2):
    """Return melting energy for association of s1/s2.
    """

    if type(s1) != StringType or type(s2) != StringType:
        raise AttributeError, "s1  and s2 must be strings"

    if len(s2) != len(s1):
        raise SequenceLengthError, "s1 and s2 are different lengths."

    s2=reverse(s2)

    return __aos.energy(ctypes.c_char_p(s1),ctypes.c_char_p(s2),ctypes.c_int(0))


def selfEnergy(s):
    """return the melting energy of 's' with itself.
    """
    return energy(s,s)

class AOSError(Exception):
    pass

class SequenceLengthError (AOSError):
    pass

class SequenceBaseError(AOSError):
    pass



def blastEnergy(qrySeq,subjSeq):
    """Return the energy for the association of the query and subject
    sequences.  The gaps are found by SW algnment.  
    """

    swAlign,swScore = SW(qrySeq,subjSeq)

    print >> sys.stderr,swAlign
    print >> sys.stderr,swScore
    print

    qrySw,subjSw =swAlign.split('\n')

    qrySw=qrySw.rstrip('.')
    subjSw=subjSw.rstrip('.')
    

    # replace gaps with -
    qrySw=qrySw.replace('.','-')
    subjSw=subjSw.replace('.','-')

    # replace sw spaces with missing seq
    while True:
        sPos = qrySw.find(' ')
        if sPos == -1:
            break
        qrySw=qrySw.replace(' ',qrySeq[sPos],1)

    while True:
        sPos = subjSw.find(' ')
        if sPos == -1:
            break
        subjSw=subjSw.replace(' ',subjSeq[sPos],1)

    # add dashes to right end if necessary
    sizeDiff = len(qrySw)-len(subjSw)
    if sizeDiff > 0:
        subjSw+='x'*sizeDiff
    elif sizeDiff < 0:
        qrySw+='x'*abs(sizeDiff)
    

    # if there is a gap in both sequences at the same place
    # (AOS SW does that sometimes) remove it (or them).
    goodPos=[]
    for i in range(len(qrySw)):
        if not(qrySw[i] == '-' and subjSw[i] == '-'):
            goodPos.append(i)
    qrySw=''.join([qrySw[i] for i in goodPos])
    subjSw=''.join([subjSw[i] for i in goodPos])
        
    try:
        e= energy(qrySw,reverseComplement(subjSw))
    except:
        print ("""energy calculation failed:

        Query:\t%s\t%s
        Subj:\t%s\t%s
        """ % (qrySeq,qrySw,subjSeq,subjSw))

        raise
        
    return e
    

#
# requirements from Fischer Lab sequence package
# here for protability
#

def reverseComplement(seq):
    """return the reverse complement of a sequence"""
    #seq=seq.upper()
    #if not isAllDRNA(seq):
    #    raise SequenceError, "Non-DRNA character in seq: %s" % (seq)

    # complement
    compl = complement(seq)
    # reverse
    return compl[::-1]

def reverse(seq):
    """return a copy of the reversed seq
    """
    rv = list(seq)
    rv.reverse()
    return ''.join(rv)

def complement(seq,transl=None):
    """Return complement of seq.
    """
    transl=DRNAcomplement
    if isAllDNA(seq):
        transl=DNAcomplement
    elif isAllRNA(seq):
        transl=RNAcomplement

    if transl == None:
        raise ValueError , "sequence is not all DNA or RNA bases." 
    
    compl = seq.translate(transl)
    return compl


def isAllDNA(seq):
    """True if all of sequence is valid unambiquious DNA.
    """
    if notDNA.search(seq) == None:
        return True
    else:
        return False

def isAllRNA(seq):
    """True if all of sequence is valid unambiquious RNA.
    """
    if notRNA.search(seq) == None:
        return True
    else:
        return False

def isAllDRNA(seq):
    """True if all of sequence is valid unambiquious RNA or DNA.
    """
    if notDRNA.search(seq) == None:
        return True
    else:
        return False
isAllNA = isAllDRNA

def isAllDRNAorSpace(seq):
    """True if all of sequence is valid RNA or DNA base or whitespace
    no ambiquity characters except n,N,x or X are allowd.
    """
    if notDRNAorSpace.search(seq) == None:
        return True
    else:
        return False
    
#
# Sequencetranslations
#
DNAcomplement = string.maketrans('aAcCgGtTnNxX-\t\n ','tTgGcCaAnNxX-\t\n ')
RNAcomplement = string.maketrans('aAcCgGuUnNxX-\t\n ','uUgGcCaAnNxX-\t\n ')
DRNAcomplement = string.maketrans('aAcCgGuUtTnNxX-\t\n ','uUgGcCaAaAnNxX-\t\n ')
NAcomplement=DRNAcomplement

#
# Sequence regular expressions
#
DNAbasesPat = re.compile(r'[ACGT]+',re.IGNORECASE)
RNAbasesPat = re.compile(r'[ACGU]+',re.IGNORECASE)
DRNAbasesPat = re.compile(r'[ACGUT]+',re.IGNORECASE)


NAbases = ('A','T','U','G','C')
DRNAbases = NAbases
DNAbases = ('A','T','G','C')
RNAbases = ('A','U','G','C')
wc_pairs = (('A','T'),
            ('T','A'),
            ('A','U'),
            ('U','A'),
            ('G','C'),
            ('C','G'))
purines = ('A','G')
pyrimidines = ('T','C','U')
bitBasesDict = {'A':0,'C':1,'G':2,'T':3,'a':0,'c':1,'g':2,'t':3,0:'A',1:'C',2:'G',3:'T'}
notNA =  re.compile('[^%s]' % string.join(NAbases,''))
notDRNA = notNA
notDNA = re.compile('[^%s]' % string.join(DNAbases,''))
notRNA = re.compile('[^%s]' % string.join(RNAbases,''))

notNAX =  re.compile('[^%s]' % (string.join(NAbases,'')+'X'))
notDRNAX = notNAX
notDNAX = re.compile('[^%s]' % (string.join(DNAbases,'')+'X'))
notRNAX = re.compile('[^%s]' % (string.join(RNAbases,'')+'X'))


notDRNAorSpace =re.compile('[^%s\s]' % string.join(NAbases,''))
lettersRE = re.compile('[%s]' % string.join(string.letters,''))
notLetters = re.compile('[^%s]' % string.join(string.letters,''))

stdStops = re.compile('TGA|TA[AG]')
stdStart = re.compile('ATG')
