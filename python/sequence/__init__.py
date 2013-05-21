"""
sequence utilities

$Id: __init__.py,v 1.46 2012/12/17 22:02:43 kael Exp $

This module is where the StringSequence base class lives.

Also there are many functions that operate either on raw strings
represnting sequence data or sequence titles.  These functions operate
on such raw strings or on StringSequence instances and they are usable
as bound instance methods too.  

"""
__version__ ="$Revision: 1.46 $"
__author__ = "Kael Fischer and others"


import copy
import functools
import math
import new
import re
import random
import string
import os
import tempfile
from types import *

from utils import mystemp,stats

class SequenceError (Exception):
    pass

class FastqParseError (Exception):
    pass

#
# Genetic Codes
#
stdGeneCode = {
    'TTT':'Phe','TCT':'Ser','TAT':'Tyr','TGT':'Cys',
    'TTC':'Phe','TCC':'Ser','TAC':'Tyr','TGC':'Cys',
    'TTA':'Leu','TCA':'Ser','TAA':'STOP','TGA':'STOP',
    'TTG':'Leu','TCG':'Ser','TAG':'STOP','TGG':'Trp',
    'CTT':'Leu','CCT':'Pro','CAT':'His','CGT':'Arg',
    'CTC':'Leu','CCC':'Pro','CAC':'His','CGC':'Arg',
    'CTA':'Leu','CCA':'Pro','CAA':'Gln','CGA':'Arg',
    'CTG':'Leu','CCG':'Pro','CAG':'Gln','CGG':'Arg',
    'ATT':'Ile','ACT':'Thr','AAT':'Asn','AGT':'Ser',
    'ATC':'Ile','ACC':'Thr','AAC':'Asn','AGC':'Ser',
    'ATA':'Ile','ACA':'Thr','AAA':'Lys','AGA':'Arg',
    'ATG':'Met','ACG':'Thr','AAG':'Lys','AGG':'Arg',
    'GTT':'Val','GCT':'Ala','GAT':'Asp','GGT':'Gly',
    'GTC':'Val','GCC':'Ala','GAC':'Asp','GGC':'Gly',
    'GTA':'Val','GCA':'Ala','GAA':'Glu','GGA':'Gly',
    'GTG':'Val','GCG':'Ala','GAG':'Glu','GGG':'Gly',    }

aaCodes = {
    'Ala':'A',
    'Arg':'R',
    'Asn':'N',
    'Asp':'D',
    'Cys':'C',
    'Gln':'Q',
    'Glu':'E',
    'Gly':'G',
    'His':'H',
    'Ile':'I',
    'Leu':'L',
    'Lys':'K',
    'Met':'M',
    'Phe':'F',
    'Pro':'P',
    'STOP':'*',
    'Ser':'S',
    'Thr':'T',
    'Trp':'W',
    'Tyr':'Y',
    'Val':'V',
    'A':'Ala',
    'C':'Cys',
    'E':'Glu',
    'D':'Asn',
    'G':'Gly',
    'F':'Phe',
    'I':'Ile',
    'H':'His',
    'K':'Lys',
    '*':'STOP',
    'M':'Met',
    'L':'Leu',
    'N':'Asp',
    'Q':'Gln',
    'P':'Pro',
    'S':'Ser',
    'R':'Arg',
    'T':'Thr',
    'W':'Trp',
    'V':'Val', 
    'Y':'Tyr' } 

#
# Sequencetranslations
#
DNAcomplement = string.maketrans('aAcCgGtTnNxX-\t\n ','tTgGcCaAnNxX-\t\n ')
RNAcomplement = string.maketrans('aAcCgGuUnNxX-\t\n ','uUgGcCaAnNxX-\t\n ')
DRNAcomplement = string.maketrans('aAcCgGuUtTnNxX-\t\n ','tTgGcCtTaAnNxX-\t\n ')
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
notLetters = re.compile('[^%s]' % string.join(string.letters+'-',''))

stdStops = re.compile('TGA|TA[AG]')
stdStart = re.compile('ATG')

NAambig = {'M':	('A', 'C'),
           'R' :('A', 'G'),
           'W' :('A', 'T'),
           'S' :('C', 'G'),
           'Y' :('C', 'T'),
           'K' :('G', 'T'),
           'V' :('A', 'C', 'G'),
           'H' :('A', 'C', 'T'),
           'D' :('A', 'G', 'T'),
           'B' :('C', 'G', 'T'),
           'N' :('A', 'C', 'G', 'T')}

NAall = list(NAbases)+ NAambig.keys()

#
# Title regexp
#
tileParts=re.compile(r'(\d+)_([n]?[t]?)(\d+)\.?(\d*)')

def sequenceWrapper(fcn):
    """call fcn with .sequence attribute of argument 0,
    if it exists.
    """
    def wrapper(*args, **kwargs):
        if hasattr(args[0],'sequence'):
            args=[args[0].sequence]+list(args[1:])
        return fcn(*args,**kwargs)
    functools.update_wrapper(wrapper,fcn)
    return wrapper


def titleWrapper(fcn):
    """call fcn with .title attribute of argument 0,
    if it exists.
    """
    def wrapper(*args, **kwargs):
        if (hasattr(args[0],'title')
            and type(args[0].title) in StringTypes):
            args=[args[0].title]+list(args[1:])
        return fcn(*args,**kwargs)
    functools.update_wrapper(wrapper,fcn)
    return wrapper
    
    
#
# Functions that operate on raw sequence (i.e. strings)
# or objects that have .sequence atributes (that are strings)
#
@sequenceWrapper
def isAllDNA(seq):
    """True if all of sequence is valid unambiquious DNA.
    """
    if notDNA.search(seq) == None:
        return True
    else:
        return False

@sequenceWrapper
def isAllRNA(seq):
    """True if all of sequence is valid unambiquious RNA.
    """
    if notRNA.search(seq) == None:
        return True
    else:
        return False

@sequenceWrapper
def isAllDRNA(seq):
    """True if all of sequence is valid unambiquious RNA or DNA.
    """
    if notDRNA.search(seq) == None:
        return True
    else:
        return False
isAllNA = isAllDRNA

@sequenceWrapper
def isAllDRNAorSpace(seq):
    """True if all of sequence is valid RNA or DNA base or whitespace
    no ambiquity characters except n,N,x or X are allowd.
    """
    if notDRNAorSpace.search(seq) == None:
        return True
    else:
        return False
    
def makeBarcodes(numBases, parityBit=False):
    """
    """
    barcodes = []
    
    for barcode in permutations(DNAbases,numBases):
        barcodes.append(barcode)

    if parityBit==True:
        for i in range(len(barcodes)):
            parityBase = bitBasesDict[reduce(lambda x,y: x^y, [bitBasesDict[z] for z in barcodes[i]])]
            barcodes[i] += parityBase
    return barcodes
        
def permutations(items, n):
    if n == 0:
        yield ''
    else:
        for i in range(len(items)):
            for base in permutations(items, n - 1):
                yield str(items[i]) + str(base)

@sequenceWrapper
def reverseComplement(seq, transl=DRNAcomplement):
    """return the reverse complement of a sequence"""
    #seq=seq.upper()
    #if not isAllDRNA(seq):
    #    raise SequenceError, "Non-DRNA character in seq: %s" % (seq)

    # complement
    compl = complement(seq, transl)
    # reverse
    return compl[::-1]

@sequenceWrapper
def reverse(seq):
    """return a copy of the reversed seq
    """
    rv = list(seq)
    rv.reverse()
    return ''.join(rv)

@sequenceWrapper
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

#
# re sequence helper
#
RePatternType = type(re.compile(''))
def _reTransformer (patterns):
    rv = []
    if type(patterns) not in (ListType,TupleType):
        patterns = (patterns,)
    
    for pattern in patterns:
        if type(pattern) == RePatternType:
            rv.append(pattern)
        elif type(pattern) in StringTypes:
            rv.append(re.compile(pattern))
    return rv
reTransformer = _reTransformer

#
# File Size Calculation Functions
#
def toString(thing):
    """takes a file object or a file name/path.
    returns a string with the contents of the file.
    files are not rewound or seek'ed to their previous
    positions.
    """
    thingType = type(thing)

    if thingType == FileType:
        return thing.read()
    elif thingType in StringTypes:
        return file(thing).read()
    

    

#
# Sequence Generation
#
def randomBase (excludeList = []):
    """Return a random DNA/RNA base, but not if it is in the excludeList.
    """
    choices = []
    
    for b in (DRNAbases):
        if b not in excludeList:
            choices.append(b)
    return random.choice(choices)


def randomBaseGenerator(pDict={'A':1, 'T':1, 'G':1, 'C':1},limit=None):
    """return a generator for random bases not in the exclude list.
    Can specify relative ferquencies for each base (this overrides excludeList),
    of the form: {'A':int,'T':int,'G',int,'C',int}.
    """
    choices =[]
    freqs=[]
    sumf=0
    
    for c,f in pDict.items():
        choices.append(c)
        if len(freqs) == 0:
            freqs.append(f)
        else:
            freqs.append(f+freqs[-1])
        sumf+=f
        

    while limit == None or limit > 0:
        n = random.randint(1,sumf)
        for i in range(len(freqs)):
            if n <= freqs[i]:
                yield choices[i]
                break
        if limit != None:
            limit -= 1
    
    

def randomSequence (length, excludeList=['U'],excludeSites=[],
                    maxAT=None,minAT=None,ATwindow=None,maxHomo=None):
    """Return a string of random DNA/RNA bases, but no bases from the
    excludeList.  Multibase sites can be excluded by specifying them
    in excludeSites.
    By default 'U' is excluded.
    """
    siteFound = True
    
    while siteFound:
        seq = []
        siteFound = False
        
        while length > len(seq):
            seq.append(randomBase(excludeList = excludeList))
        sStr = string.join(seq,'')
        
        # check for excluded sites
        for baddie in excludeSites:
            if sStr.find(baddie) != -1:
                siteFound = True
                break
            
        # check for TA
        if  siteFound == False and (maxAT != None or minAT != None):
            if ATwindow==None:
                ATwindow = length
            for start in range(length-ATwindow+1):
                window = sStr[start:ATwindow+start]
                ATcount =  window.count("A")+window.count("T")
                if maxAT != None:
                    window = sStr[start:ATwindow+start]
                    ATcount =  window.count("A")+window.count("T")
                    while ATcount > maxAT:
                        ats = []

                        idx=0
                        while idx != -1:
                            idx = window.find('A',idx+1)
                            if idx != -1:
                                ats.append(idx)

                        idx=0
                        while idx != -1:
                            idx = window.find('T',idx+1)
                            if idx != -1:
                                ats.append(idx)

                        target = random.choice(ats)
                        target += start
                        replacement = random.choice(('C','G'))
                        sStr=sStr[:start+target]+replacement+sStr[start+target+1:]
                        window = sStr[start:ATwindow+start]
                        ATcount =  window.count("A")+window.count("T")




                if minAT != None:
                    window = sStr[start:ATwindow+start]
                    ATcount =  window.count("A")+window.count("T")
                    while ATcount < minAT:
                        gcs = []

                        idx=0
                        while idx != -1:
                            idx = window.find('C',idx+1)
                            if idx != -1:
                                gcs.append(idx)

                        idx=0
                        while idx != -1:
                            idx = window.find('G',idx+1)
                            if idx != -1:
                                gcs.append(idx)
                        
                        target = random.choice(gcs)
                        target += start
                        replacement = random.choice(('A','T'))
                        sStr=sStr[:start+target]+replacement+sStr[start+target+1:]
                        window = sStr[start:ATwindow+start]
                        ATcount =  window.count("A")+window.count("T")


                
        if siteFound == False and (maxHomo != None and lenHomo(sStr)>maxHomo):
            siteFound=True                         
    return sStr


@sequenceWrapper
def expandAmbiguousSequence(sequence):
    """Returns an iterator with all the posible replacements of
    ambiguous bases.
    """
    sequence = sequence.upper()
    queue = [(sequence, 0)]
    nambig = NAambig.keys()
    while len(queue) > 0:
        seq, pos = queue.pop(0)
        if seq[pos] in nambig:
            for na in NAambig[seq[pos]]:
                seqTemp = list(seq)
                seqTemp[pos] = na
                queue.append((''.join(seqTemp), pos))
        elif pos == len(seq) - 1:
            yield seq
        else:
            queue.append((seq, pos + 1))


@sequenceWrapper
def disambiguateSequenceRandom(sequence,maxReplacements=None):
    """Returns an a sequence with ambiguous bases, replaced with a random
    base from the ambiguity spec.
    if maxReplacements is defined and there are more than that
    number of ambiguities a ValueError raised
    """
    ct=0
    nambig = NAambig.keys()
    letters = list(sequence)
    for i,b in enumerate(letters):
        ct+=1
        if (maxReplacements is not None) & (ct>maxReplacements):
            raise ValueError, (
                "sequence hase more than %s ambiguities" % maxReplacements)
        if b in nambig:
            letters[i]=random.choice(NAambig[b])
    return ''.join(letters)




@sequenceWrapper
def bashSequence(seqString,mutablePositions,
                 excludeList=[],
                 allowSameBase=False,
                 allowSameRings=True,
                 allowNonBase=True):
    """Retuns a mutated DNA sequence (as a string).
    mutable postions are specified as a single index or a
    list or tupple of same.
    """

    if type(mutablePositions) == IntType or \
       type(mutablePositions) == LongType:
        mutablePositions = [mutablePositions]

    
    seq = []
    for c in seqString:
        if allowNonBase == False:
            if c not in bases:
                raise SequenceError, "%s at position %s is not a DNA base, try one of %s" % \
                      (c,seqString.find(c),bases)
        seq.append(c)
        

    for pos in mutablePositions:
        
        if allowSameRings == False:
            if seq[pos] in purines:
                seq[pos] = randomBase(excludeList=purines+excludeList)
            elif seq[pos] in pyrimidines:
                seq[pos] = randomBase(excludeList=purines+excludeList)
            elif seq[pos] not in bases :
                raise SequenceError, "%s at position %s is not a DNA base, try one of %s" % \
                      (seq[pos],pos,bases)
            else:
                print seq[pos]
                raise Exception, "this should never happen, evacuate the building!\n************NOW**********\n"
        elif allowSameBase == True:
            seq[pos] = randomBase(excludeList=excludeList)

        else:
            seq[pos] = randomBase(excludeList=[seq[pos]]+excludeList)

    return string.join(seq,'')

@sequenceWrapper
def longestStretch(seqString, bases=()):
    """Returns the longest stretch of the sequence that contains
    only the bases in the string,tupple or list bases.  Returns the
    first sequence if more that one strech are tied for longest.
    """
    
    longest = ''
    
    if len(bases) == 0:
        return longest
    
    reBases = re.compile(r'[%s]+' % (string.join(bases,'')))
    matches = reBases.findall(seqString)
    for m in matches:
        if len(m) > len(longest):
            longest = m

    return longest
        

def xOutMatch(matchObj, maxLen=0):
    """return a string of Xes as long as a match object's match.
    If maxLen > 0 matches longer than maxLen are returned unchanged.  
    """
    mLen = matchObj.end()-matchObj.start()
    if maxLen == 0 or (maxLen > 0 and  mLen <= maxLen):
        return 'X' * mLen
    else:
        return matchObj.group(0)


@sequenceWrapper
def wordCount(seqString,wordSize):
    rv={}
    for s in xrange(0,len(seqString)+1-wordSize):
        e=s+wordSize
        try:
            rv[seqString[s:e]]+=1
        except KeyError:
            rv[seqString[s:e]]=1
    return rv
    
@sequenceWrapper
def screenRe(seqString,pattern,maxLen=0):
    """Replace occurences of the regex pattern with X,
    return resulting string.  pattern can be a string or a compiled re.
    If maxLen > 0 matches longer than maxLen  are unchanged.
    """

    replFcn = lambda x: xOutMatch(x,maxLen=maxLen)
    return re.sub(pattern,replFcn,seqString)


@sequenceWrapper
def screenShorties(seqString,shortyLength=50):
    """Replace stretches of valid RDNA bases <= shortyLength
    with the same number of X's.
    """
    return screenRe(seqString,'[ACGTU]+',maxLen=shortyLength)
    
    
@sequenceWrapper
def dropNonLetters(seqString):
    """remove non-letters from seqString, return new string.
    """
    return re.sub(notLetters,'',seqString)


@sequenceWrapper
def rePositions(seqString,pattern):
    """returns a list of positions of non-overlapping
    occurances of pattern.  Positions are in bio numbering (start=1).
    Note: non-overlapping!!
    """
    return [ m.start()+1 for m in  re.finditer(pattern,seqString)] 
    

#
# Useful Sequences
#
primerB   = 'GTTTCCCAGTCACGATA'
primerBrc = reverseComplement(primerB)

SPIKE70 = 'ACCTCGCTAACCTCTGTATTGCTTGCCGGACGCGAGACAAACCTGAACATTGAGAGTCACCCTCGTTGTT'
PROBE70 = reverseComplement(SPIKE70)

M13Forward = 'GTAAAACGACGGCCAG'
M13Reverse = 'CAGGAAACAGCTATGAC'

T7Promoter = 'TAATACGACTCACTATAGGG'

#
# other useful functions
#

#
# Functions that operate on raw strings or SequenceRecord.sequence
# these are (should be) accessible as StringSequence object methods.
#
# -NEW Easy Wrapper Methods-
#
# To add a wrapper method to a new function in StringSequence
# and its subclasses add the name of the function to
# StringSequence._fcnToWrap
#
@sequenceWrapper
def Tm(s,dnac=50,saltc=50,rna=False,debug=False):
    """Returns DNA/DNA tm using nearest neighbor thermodynamics. dnac is
    DNA concentration [nM] and saltc is salt concentration [mM].
    rna=0 is for DNA/DNA (default), for RNA, rna should be True.
    Sebastian Bassi <sbassi@genesdigitales.com>"""
    
    #Credits: 
    #Main author: Sebastian Bassi <sbassi@genesdigitales.com>
    #Overcount function: Greg Singer <singerg@tcd.ie>
    #Based on the work of Nicolas Le Novere <lenov@ebi.ac.uk> Bioinformatics. 17:1226-1227(2001)

    #This function returns better results than EMBOSS DAN because it uses updated
    #thermodinamics values and take into account inicialization parameters from SantaLucia
    #works (1998).
    
    #Things to do:
    #+Add a function to detect complementary sequences. Change K according to result.
    #+Add support for heteroduplex (see Sugimoto et al. 1995).
    #+Correction for Mg2+. Now supports only monovalent ions.
    #+Put thermodinamics table in a external file for users to change at will
    #+Add support for danglings ends (see Le Novele. 2001) and mismatches.
    
    dh=0 #DeltaH. Enthalpy
    ds=0 #deltaS Entropy

    def tercorr(stri):
        deltah=0
        deltas=0
        if rna==0:
            #DNA/DNA
            #Allawi and SantaLucia (1997). Biochemistry 36 : 10581-10594
            if stri[0]=="G" or stri[0]=="C":
                deltah=deltah-0.1
                deltas=deltas+2.8
            elif stri[0]=="A" or stri[0]=="T":
                deltah=deltah-2.3
                deltas=deltas-4.1
            if stri[-1]=="G" or stri[-1]=="C":
                deltah=deltah-0.1
                deltas=deltas+2.8
            elif stri[-1]=="A" or stri[-1]=="T":
                deltah=deltah-2.3
                deltas=deltas-4.1
            dhL=dh+deltah
            dsL=ds+deltas
            return dsL,dhL
        elif rna==1:
            #RNA
            if stri[0]=="G" or stri[0]=="C":
                deltah=deltah-3.61
                deltas=deltas-1.5
            elif stri[0]=="A" or stri[0]=="T" or stri[0]=="U":
                deltah=deltah-3.72
                deltas=deltas+10.5
            if stri[-1]=="G" or stri[-1]=="C":
                deltah=deltah-3.61
                deltas=deltas-1.5
            elif stri[-1]=="A" or stri[-1]=="T" or stri[0]=="U":
                deltah=deltah-3.72
                deltas=deltas+10.5
            dhL=dh+deltah
            dsL=ds+deltas
            if debug:
                print "delta h=",dhL
            return dsL,dhL

    def overcount(st,p):
        """Returns how many p are on st, works even for overlapping"""
        ocu=0
        x=0
        while 1:
            try:
                i=st.index(p,x)
            except ValueError:
                break
            ocu=ocu+1
            x=i+1
        return ocu

    sup=string.upper(s)
    R=1.987 # universal gas constant in Cal/degrees C*Mol
    vsTC,vh=tercorr(sup)
    vs=vsTC
    
    k=(dnac/4.0)*1e-8
    #With complementary check on, the 4.0 should be changed to a variable.
    
    if rna==0:
        #DNA/DNA
        #Allawi and SantaLucia (1997). Biochemistry 36 : 10581-10594
        vh=(vh+(overcount(sup,"AA"))*7.9+(overcount(sup,"TT"))*7.9
            +(overcount(sup,"AT"))*7.2+(overcount(sup,"TA"))*7.2
            +(overcount(sup,"CA"))*8.5+(overcount(sup,"TG"))*8.5
            +(overcount(sup,"GT"))*8.4+(overcount(sup,"AC"))*8.4)
        vh=(vh+(overcount(sup,"CT"))*7.8+(overcount(sup,"AG"))*7.8
            +(overcount(sup,"GA"))*8.2+(overcount(sup,"TC"))*8.2)
        vh=(vh+(overcount(sup,"CG"))*10.6+(overcount(sup,"GC"))*10.6
            +(overcount(sup,"GG"))*8+(overcount(sup,"CC"))*8)
        vs=(vs+(overcount(sup,"AA"))*22.2+(overcount(sup,"TT"))*22.2
            +(overcount(sup,"AT"))*20.4+(overcount(sup,"TA"))*21.3)
        vs=(vs+(overcount(sup,"CA"))*22.7+(overcount(sup,"TG"))*22.7
            +(overcount(sup,"GT"))*22.4+(overcount(sup,"AC"))*22.4)
        vs=(vs+(overcount(sup,"CT"))*21.0+(overcount(sup,"AG"))*21.0
            +(overcount(sup,"GA"))*22.2+(overcount(sup,"TC"))*22.2)
        vs=(vs+(overcount(sup,"CG"))*27.2+(overcount(sup,"GC"))*27.2
            +(overcount(sup,"GG"))*19.9+(overcount(sup,"CC"))*19.9)
        ds=vs
        dh=vh
        
    else:
        #RNA/RNA hybridisation of Xia et al (1998)
        #Biochemistry 37: 14719-14735         
        vh=(vh+(overcount(sup,"AA"))*6.82+(overcount(sup,"TT"))*6.6
            +(overcount(sup,"AT"))*9.38+(overcount(sup,"TA"))*7.69
            +(overcount(sup,"CA"))*10.44+(overcount(sup,"TG"))*10.5
            +(overcount(sup,"GT"))*11.4+(overcount(sup,"AC"))*10.2)
        vh=(vh+(overcount(sup,"CT"))*10.48+(overcount(sup,"AG"))*7.6
            +(overcount(sup,"GA"))*12.44+(overcount(sup,"TC"))*13.3)
        vh=(vh+(overcount(sup,"CG"))*10.64+(overcount(sup,"GC"))*14.88
            +(overcount(sup,"GG"))*13.39+(overcount(sup,"CC"))*12.2)
        vs=(vs+(overcount(sup,"AA"))*19.0+(overcount(sup,"TT"))*18.4
            +(overcount(sup,"AT"))*26.7+(overcount(sup,"TA"))*20.5)
        vs=(vs+(overcount(sup,"CA"))*26.9+(overcount(sup,"TG"))*27.8
            +(overcount(sup,"GT"))*29.5+(overcount(sup,"AC"))*26.2)
        vs=(vs+(overcount(sup,"CT"))*27.1+(overcount(sup,"AG"))*19.2
            +(overcount(sup,"GA"))*32.5+(overcount(sup,"TC"))*35.5)
        vs=(vs+(overcount(sup,"CG"))*26.7+(overcount(sup,"GC"))*36.9
            +(overcount(sup,"GG"))*32.7+(overcount(sup,"CC"))*29.7)
        ds=vs
        dh=vh

    ds=ds-0.368*(len(s)-1)*math.log(saltc/1e3)
    tm=((1000* (-dh))/(-ds+(R * (math.log(k)))))-273.15
    if debug:
        print "ds="+str(ds)
        print "dh="+str(dh)
        print "Tm="+str(tm)
        
    return tm

@sequenceWrapper
def sequenceMatchLengths(sequence,patterns):
    """return a list of the lengths of matches to patterns
    in the input sequence (either a record or a sequence).
    """

    patterns = reTransformer(patterns)
    rv = []
    for pattern in patterns:
        mLengths = map(len,pattern.findall(sequence))
        if len(mLengths) == 0:
            rv.append(0)
        else:
            rv.extend(mLengths)
    if len(patterns) == 1 and False:
        return rv[0]
    else:
        return rv
    
@sequenceWrapper
def DNABaseLengths(s):
    """return a list of the lengths of DNA base streatches in
    the sequence. 
    """
    return sequenceMatchLengths(s,DNAbasesPat)

@sequenceWrapper
def DNABaseCount(s):
    """return the number of DNA bases in
    the sequence. 
    """
    return sum(sequenceMatchLengths(s,DNAbasesPat))

@sequenceWrapper
def RNABaseLengths(s):
    
    """return a list of the lengths of RNA base stretches in
    the sequence. 
    """
    return sequenceMatchLengths(s,RNAbasesPat)

@sequenceWrapper
def RNABaseCount(s):
    """return the number of RNA bases in
    the sequence. 
    """
    return sum(sequenceMatchLengths(s,RNAbasesPat))

@sequenceWrapper
def DRNABaseLengths(s):
    """return a list of the lengths of RNA or DNA base stretches in
    the sequence. 
    """
    return sequenceMatchLengths(s,DRNAbasesPat)

@sequenceWrapper
def DRNABaseCount(s):
    """return the number of RNA and DNA bases in
    the sequence. 
    """
    return sum(sequenceMatchLengths(s,DRNAbasesPat))

@sequenceWrapper
def LZWsize(s):
    """return the size of the LZW compressed sequence
    """
    import aos
    return aos.LZWsize(s)

@sequenceWrapper
def LZWratio(s):
    """Return the LZW compression ratio.
    """
    import aos
    try:
        if len(s) == 0:
            raise aos.SequenceLengthError
        return float(LZWsize(s))/float(len(s)) 
    except aos.SequenceLengthError:
        print "Divide by Zero Error: length of sequence is zero!"

@sequenceWrapper
def lenHomo(s):
    """returns the length of the longest homopolymeric stretch
    in s"""
    letDict = dict(zip(s,[None]*len(s)))

    maxHomo = 0
    for letter in letDict.keys():
        i = len(longestStretch(s,letter))
        if i > maxHomo:
            maxHomo = i
    return maxHomo
    
@sequenceWrapper
def pctCG(s):
    """returns a float with the %CG of the sequence.
    """
    n=len(s)
    return float(s.count('c')+s.count('C')+s.count('g')+s.count('G'))/float(n)
pctGC = pctCG
    
@sequenceWrapper
def baseComp(s,letters='ACGTN'):
    """returns counts of of residues / bases in letters.
    """
    rd = {}
    for l in letters:
        rd[l]=s.count(l)
    return rd

@sequenceWrapper
def residueCount(s):
    """return a dict with the count of all letters in s
    
    Arguments:
    - `s`:
    """
    rv={}
    for i in xrange(len(s)):
        if s[i] in rv:
            continue
        rv[s[i]]=s.count[s[i]]
    return rv
    
@sequenceWrapper
def find(seq,s):
    """return a list of indexes of start positions of subsequence.
    s occorances may overlap. (postions are in python indexing: 0 start)
    """
    rv=[]
    i=0
    while True:
        try:
            rv.append(seq.index(s,i))
        except ValueError:
            return rv
        i=rv[-1]+1

@sequenceWrapper
def findRC(seq,s):
    """return a list of indexes of start positions of subsequence
    reverse complement. s occorances may overlap. (postions are
    in python indexing: 0 start)
    """
    s=reverseComplement(s)
    return find(seq,s)

#
# Functions that operate on sequence titles
#
#
# -Easy Function Wrappers for StringSequence (et al)-
#
# To wrap functions that are title manipulating functions
# add the functions' name to StringSequence._titleFcnToWrap
#

@titleWrapper
def splitTileTitle(title):
    """splits a record's title into (gi, start, end).  Any of those
    cannot be determined are set to None.  Any non-None values
    are integers.

    Note to get the sequence using slice notation use start=start-1
    """
    if type(title) not in StringTypes:
        raise TypeError, 'title must be StringType not: %s' % title
        
    rvL=[None,None,None]

    m=tileParts.match(title)
    if m!= None:
        rvL[0] = int(m.group(1))
        if m.group(2)=='':
            # acient 'segment' number format
            pass
        else:
            try:
                rvL[1] = int(m.group(3))
            except:
                pass
            try:
                rvL[2] = rvL[1]+int(m.group(4))-1
            except:
                pass

    return tuple(rvL)

@titleWrapper
def isTileTitle(title):
    """True if title is in the format of a tile
    """
    try:
        m=tileParts.match(title)
        return m!=None
    except TypeError:
        return False

@titleWrapper
def countFromTitle(title):
    """returns the terminal number from a tile that ends with '|####'.
    where #### is a integer.  If the split or numeric conversion fails,
    None is returned
    """
    try:
        return int(title.split('|')[-1])
    except (IndexError,ValueError):
        return None

occurrenceCount=countFromTitle

@titleWrapper
def hiSeqBarcode(title):
    """returns the bar code encoded in a HiSeq read title
    """
    try:
        return re.match('.*#([ACTG]+)',title).group(1)
    except AttributeError:
        return None

@titleWrapper
def fastaTitleFields(title):
    """Split NCBI fasta title into fields.
    return dictionary with keys like 'gi', 'ref', 'gb', etc
    """
    rv = {}
    f = title.split('|')
    start = len(f)%2
    if f[0] == 'gi':
        start = 0
    for i in range(start,len(f),2):
        try:
            rv[f[i]] = f[i+1]
        except:
            rv[f[i]]=''
    return rv


#
# Functions to be wrapped by StringSequence methods
# should be defined above.
#

class StringSequence (object):
    """Base class for sequences objects that represent sequence information as
    strings in a .sequence property.
    Provides stringy stuff like slicing.
    """

    def __init__(self,sequence='',title='',colwidth=60):
        self.title=title
        self.sequence=sequence
        self.colwidth=colwidth
        
    def __len__(self):
        """Return the length of sequence.
        """
        return len(self.sequence)

    length = __len__

    def __getitem__(self,item):
        return self.sequence[item]
    
    def wrappedSequence(self):
        """returns the wordwrapped sequence.
        """
        s=[]
        i = 0
        while i < len(self):
            s.append(self[i:i+self.colwidth])
            i = i + self.colwidth
        return os.linesep.join(s)


    def scramble(self):
        """randomly permute sequence.
        """
        self.sequence = ''.join(random.sample(self.sequence,len(self.sequence))) 
    

    def screenedLen(self,screenChrs=None):
        """Return the length of the sequence, less any
        of the characters in screenChrs.
        """
        if screenChrs == None:
            return len(self)

        count = 0
        for c in screenChrs:
            count += self.sequence.count(c)
        return len(self)-count


    def matchLengths(self,patterns):
        """return a list of lengths of the sequences matches to
        each pattern.  If patterns is a sequence of patterns a
        list of lists of matching lengths is returned otherwise
        an unnested single list of match lengths is returned.
        """
        return sequenceMatchLengths(self,patterns)

    def reSearch(self,patterns):
        """Return a list of Match Objects to the patterns given.
        If patterns is a sequence of patterns a
        list of lists of Match Objects is returned, otherwise
        an unnested single list of Match Objects is returned.
        """
        patterns = reTransformer(patterns)
        rv = []
        for pattern in patterns:
            matches = list(pattern.finditer(self.sequence))
            if len(matches) == 0:
                rv.append([])
            else:
                rv.append(matches)

        if len(patterns) == 1:
            return rv[0]
        else:
            return rv

    def screenShorties(self,shortyLength=50):
        """replace short stretches (<= shortyLength) of
        valid bases with X's
        """
        if shortyLength > 0:
            self.sequence=screenShorties(self.sequence,
                                         shortyLength=shortyLength)

    def tabFormat(self):
        """return string: <title>\t<sequence>
        """
        return "%s\t%s"%(self.title,self.sequence)

    ################################################################
    # module function wrappers
    #
    #
    _seqFcnToWrap = ['DNABaseCount','DNABaseLengths',
                     'DRNABaseCount','DRNABaseLengths',
                     'LZWratio','LZWsize','RNABaseCount',
                     'RNABaseLengths','Tm','baseComp','complement',
                     'expandAmbiguousSequence','disambiguateSequenceRandom',
                     'find','findRC','isAllDNA','isAllDRNA',
                     'isAllDRNAorSpace','isAllNA','isAllRNA',
                     'lenHomo','longestStretch',
                     'pctCG','pctGC','residueCount','reverse',
                     'reverseComplement','screenShorties',
                     'sequenceMatchLengths','dropNonLetters',
                     'splitTileTitle','isTileTitle',
                     'countFromTitle','occurrenceCount',
                     'tabFormat','wordCount','hiSeqBarcode']
    
    
    def __getattr__(self,name):
        if name in self._seqFcnToWrap:
            return functools.update_wrapper(functools.partial(
                eval(name),self),eval(name))
        else:
            raise AttributeError
    #
    #
    ################################################################


    def blastSequence(self):
        import blast
        import kdbom2.exceptions
        try:
            return blast.Sequence(Name=self.m8name())
        except kdbom2.exceptions.KdbomLookupError:
            return None

    def __getitem__(self,y):
        """
        """
        return self.sequence[y]

    def isTile(self):
        """True if instance's title is in the format of a tile
        """
        m=tileParts.match(self.title)
        return m!=None

    def splitTileTitle(self):
        """splits a record's title into (gi, start, end).  Any of those
        cannot be determined are set to None.
        """
        return splitTileTitle(self.title)

    def subrecord(self,start,end):
        """Return a subsequence (aka "tile") as specified by
        start and end. Base 1 = the first base, un like when
        you use slice notation.

        start must be less than end and > 0.

        title is systematicall generated like:
        <oldtitle>_nt<start>.<length>
        """
        if start < 0 or end <= start:
            raise ValueError, "start must be less than end and > 0."
        return self.__class__(title = "%s_nt%s.%s" %
                              (self.title, start, end-start+1),
                              sequence = self[start-1:end])


    def rcRecord(self,title=None):
        """Return a reverse complement Record.
        by default the title will be title+'_rc'
        """
        if title == None:
            title = '%s_rc' % self.title
        return self.__class__(title=title,sequence=self.reverseComplement())

    
class SequenceCollection (list):
    """A collection of sequences, which can report statistics about
    contained sequences 
    """


    def statTracker(self,seqAttr):
        """return a stats.Tracker instance for the given sequence attribute,
        e.g. length, pctGC, etc
        
        Arguments:
        - `self`:
        - `seqAttr`: string name of attribute to accumulate values/stats for 
        """
        tracker=stats.Tracker(lambda x: getattr(x,seqAttr)())
        tracker.append(self)
        return tracker
        

class Translator (object):
    """Takes an iterable NA sequence as a string or a StringSequence
    instance and provides an iteratior of the encoded amino acids. 
    """

    def __init__(self,seq,frame,gc=stdGeneCode,
                 dispDict=aaCodes):
        """make the Translator.
        
        Arguments:
        - `self`:
        - `seq`:
        - `frame`:
        """
        if abs(frame) >3:
            raise ValueError, "frame between -3 and 3"

        if frame > 0:
            self.seq = iter(seq)
        elif frame < 0:
            self.seq = reverseComplement(seq)
        else:
            raise ValueError, "frame must be > or < 0"

        self.frame = frame
        for i in range(abs(frame)-1):
            self.seq.next()


        self.gc=copy.copy(gc)
        if dispDict!= None:
            for codon,aa in self.gc.items():
                self.gc[codon] = dispDict[self.gc[codon]]
        


    def getCodon(self):
        """
        """
        return ''.join((self.seq.next(),
                        self.seq.next(),
                        self.seq.next()))

    def __iter__(self):
        """
        
        Arguments:
        - `self`:
        """
        while True:
            yield self.gc[self.getCodon()]
        

def stopPositions(seq):
    """Finds std genetic code stops: TAA, TGA and TAG
    in 3 forward frames.
    Return value is: ((frame1 stops), (frame2 stops), (frame3 stops)).
    Bio numbering is used (first base = 1)
    Positions are in accending order.
    """
    p = rePositions(seq,stdStops)
    return (
        tuple([x for x in p if x%3==0]),
        tuple([x for x in p if x%3==1]),
        tuple([x for x in p if x%3==2]),
        )
    
   
        
class ORFFinder(object):
    """Find ORFs in given frames (all 6 by default).
    Start and end positions (in biology type sequence #s python 0 = 1)
    are given in a dictionary (self.frames) after you make the object
    """

    def __init__(self,seq,frames=(1,2,3,-1,-2,-3),
                 minAA=50,maxAA=None ):
        """
        """
        self.seq=seq
        self.maxAA=maxAA
        self.minAA=minAA
        
        self.frames={}
        for f in frames:
            self.frames[f]=self.getFrameORFs(f,minAA=self.minAA,
                                             maxAA=self.maxAA)


    def getFrameORFs(self,frame,minAA=50,maxAA=None):
        """
        """
        rv=[]
        start=None
        end=None

        for pos,aa in enumerate(Translator(self.seq,frame,
                                           gc=stdGeneCode,
                                           dispDict=aaCodes)):
            print pos+1,aa
            if aa != '*':
                if start == None:
                    start = pos+1
                    print "start at %s" % start
            else:
                if start==None:
                    continue
                end = pos
                print "end at %s" % end
                length = end-start+1
                print "length %s" % length
                if (length >= minAA  and
                    (maxAA == None or length <= maxAA)):
                    rv.append((start,end))
                start=None
                end=None
        if start != None:
            length = pos+1-start
            print length
            if (length >= minAA  and
                (maxAA == None or length <= maxAA)):
                rv.append((start,end))
        
        return rv

    

        
        
        
      
        
    

    
