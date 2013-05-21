"""SAM/BAM Utilities that are free of DB connections.
"""
__version__ = tuple([int(ver) for ver in
                     "$Revision: 1.1 $".split()[1].split('.')])
__author__ = "Kael Fischer"

import sys
from __init__ import *
from utils import *


def samTupleGenerator(samFile,parseTags=False,skipBits=0):
    """Iterator over sam result alignments:
    
    (qName,flag,rName,pos,mapQ,cigar,rNext,pNext,tLen,seq,qual,tags)

    see: http://samtools.sourceforge.net/SAM1.pdf
    
    Skips comments and common SGE lines of output. 

    If parseTags is True, a dictionaty of optional tags is returned,
    otherwise the unparsed tag string is returned in the tags slot.
    
    """
    
    lineN=0
    for l in multiFile(samFile):
        lineN+=1
        f=l.strip().split('\t')
        try:
            (qName,flag,rName,pos,mapQ,cigar,
             rNext,pNext,tLen,seq,qual)   =  f[:11]
            flag=int(flag)
            if skipBits & flag != 0:
                continue
            if not parseTags:
                tags = '\t'.join(f[11:])
            else:
                tags={}
                for opt in f[11:]:
                    tag,typ,val = opt.split(':')
                    if typ in ('AZ'):
                        pass
                    elif typ=='i':
                        val = int(val)
                    elif typ=='f':
                        val = float(val)
                    elif val == 'H':
                        #hex byte array
                        pass
                    elif val == 'B':
                        # other numeric array
                        pass
                    tags[tag]=val
                
            
        except ValueError:
            if (l.startswith('@') or l.startswith('#') or
                l.startswith('Reported ') or
                l.startswith('Warning:') or
                l.startswith('Thus no job control in this shell.')):
                continue
            else:
                sys.stderr.write('parsing error on line %d:\n\t%s\n'%
                                 (lineN,l))
                raise 

        pos,mapQ,pNext,tLen = map(int,(pos,mapQ,pNext,tLen))

        yield (qName,flag,rName,pos,mapQ,cigar,
               rNext,pNext,tLen,seq,qual,tags) 


