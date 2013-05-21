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

__version__ = "$Revision: 1.13 $"
__author__ = "Kael Fischer"

import sys
from types import *
import ctypes
import os.path
import sequence

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
    s2=sequence.reverse(s2)
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
    return SW(s,sequence.reverseComplement(s))



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

    s2=sequence.reverse(s2)

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



def blastEnergy(qrySeq,subjSeq,findGaps=True,debug=False):
    """Return the energy for the association of the query and subject
    sequences.  The gaps are found by SW algnment.  
    """
    if debug:
        dbOut=file('beDebug.log','w')

    dbOut=sys.stdout

    if findGaps:
        swAlign,swScore = SW(qrySeq,subjSeq)

        #debug=False

        qrySw,subjSw = [line.strip().replace(' ','x')
                        for line in swAlign.split('\n')]

        if 'x' in qrySw or 'x' in subjSw :
            debug=True

        if debug:
            print >> dbOut,  "\ninput:"
            print >> dbOut,  qrySeq
            print >> dbOut,  subjSeq

            print >> dbOut,  "SW:"
            print >> dbOut,  qrySw
            print >> dbOut,  subjSw
            print >> dbOut,  swScore


        # add x to right end if necessary
        # this should not happen
        qrySw=qrySw.rstrip('.')
        subjSw=subjSw.rstrip('.')
        sizeDiff = len(qrySw)-len(subjSw)
        if sizeDiff > 0:
            subjSw+=['x']*sizeDiff
        elif sizeDiff < 0:
            qrySw+=['x']*abs(sizeDiff)

        qrySw=list(qrySw)
        subjSw=list(subjSw)

        #put mismatches back in
        for i in range(len(qrySw)):
            if qrySw[i] == '.':
                #print >> dbOut,  qrySw[:i].count('x')]
                qrySw[i]=qrySeq[i-qrySw[:i].count('x')]
        for i in range(len(subjSw)):
            if subjSw[i] == '.':
                subjSw[i]=subjSeq[i-subjSw[:i].count('x')]

        # go back to strings
        qrySw=''.join(qrySw)
        subjSw=''.join(subjSw)

        if debug:
            print >> dbOut,  "e-mangle:"
            print >> dbOut,  qrySw
            print >> dbOut,  subjSw

    else:
        qrySw=qrySw.replace('U','T')
        subjSw=subjSw.replace('U','T')
        qrySw=qrySw
        subjSw=subjSeq

    
    try:
        e= energy(qrySw,sequence.reverseComplement(subjSw))
    except:
        print ("""energy calculation failed:

        Query:\t%s\t%s
        Subj:\t%s\t%s
        """ % (qrySeq,qrySw,subjSeq,subjSw))

        raise
        
    return e
    

        

    
