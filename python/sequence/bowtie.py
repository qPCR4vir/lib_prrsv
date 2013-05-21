"""Bowtie Utilities.
"""
__version__ = tuple([int(ver) for ver in
                     "$Revision: 1.2 $".split()[1].split('.')])

__author__ = "Kael Fischer"

from types import StringTypes

from __init__ import *
from utils import *

from blastNoSQL import giFromM8name

def bowtieTupleGenerator(btrFile):
    """Iterator over default bowtie output; return tuples of this form: 

    (readName,direction,refName,position,
     readSeq,readQuals,iaCount,mmTxt)
     
    Skips comments and common SGE lines of output.
    """
    for l in multiFile(btrFile):
        if (l.startswith('#') or l.startswith('Mega BLAST run') or
            l.startswith('Warning: no access to tty (Bad file descriptor).') or
            l.startswith('Thus no job control in this shell.')):
            continue
        try:
            yield l2t(l)
        except ValueError:
            raise


def l2t (line):
    """parse a line to a low overhead  tuple, where
    the elements have the correct data type
    """
    f=line.strip().split('\t')
    if len(f) == 8:
        (readName,direction,refName,position,
         readSeq,readQuals,iaCount,mmTxt) =  line.strip().split('\t')

    elif len(f) == 7:
        (readName,direction,refName,position,
         readSeq,readQuals,iaCount) =  line.strip().split('\t')
        mmTxt=''
    else:
        raise ValueError, 'Fewewr than 7 fields: %s' % l

    mm = tuple(mmTxt.split(','))
    position=int(position)


    return (readName,direction,refName,position,
            readSeq,readQuals,iaCount,mm)

def readTitleSet(rsltPaths,qrySet=False,subjSet=False):

    """return a set of query titles present in
    one or more bowtie alignment files
    """
    if type(rsltPaths) in StringTypes:
        rsltPaths = [rsltPaths,]
    if not qrySet and not subjSet:
        raise ValueError, "Either qrySet or subjSet must be True"
   
    rv=set()
    idxs=[]
    if qrySet:
        idxs.append(0)
    if subjSet:
        idxs.append(2)
        
    for btFile in rsltPaths:
        for l in file(btFile):
            for i in idxs:
                rv.add(l.split('\t')[i])
    return rv

