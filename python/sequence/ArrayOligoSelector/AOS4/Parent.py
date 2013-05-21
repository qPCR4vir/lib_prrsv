#!/usr/bin/env python
#
# ArrayOligoSelector
# oligo parent sequences 
#
# A sequence that oligos are specifed relative to
#
#

#
# 	$Id: Parent.py,v 1.2 2006/03/02 06:41:43 kael Exp $	
#

from Oligo import Oligo


# package level defaults
import __init__ as AOS
DEBUG = AOS.DEBUG

class Parent:
    """
    """
    def __init__ (self,pid=None,seq=None, fastaFilePath=None):
        """
        """

        self.pid=pid
        self.seq=seq
        

    def makeOligos(self,oLength):
        """
        """

        def addSeqNumber(base, startNt):
            return "%s_nt%s" % (base, startNt)
        uids = map(addSeqNumber,[self.pid]*(len(self.seq)-oLength+1),
                   range(len(self.seq)-oLength+1))


        
        return map(Oligo,uids,[self]*len(uids),
                   range(len(self.seq)-oLength+1),[oLength]*len(uids))

        

class ORF (Parent):
    pass

class Genome (Parent):
    pass

class GinormousMetaGenomicSpace (Parent):
    pass


