#!/usr/bin/env python
#
# ArrayOligoSelector
# oligo candidate classes
#
#

#
# 	$Id: Oligo.py,v 1.3 2006/03/02 18:24:31 kael Exp $	
#

# package level defaults
import __init__ as AOS
DEBUG = AOS.DEBUG


class Oligo:
    """AOS Oligo -
    Tracks position of xx-mer in a parent sequence.
    Tracks status of flagging by ScoreCards.
    Can retrieve scores from ScoreCards.
    """

    def __init__ (self, uid, parent=None, start=None, length=0):
        """
        """

        # should we track strand here?
        
        self.uid=uid
        self.parent=parent
        self.start=start
        self.length=length
        
        # indexing information for score lookup.
        # keyed by scorecard object; value is my postion
        # that scorecard's oligo list.
        self.scIndexes={}

        # has a sc flagged me
        self.scFlagged = {}


    def seq (self):
        """Returns a string coresponding to the oligo sequence
        """

        # this is a method so we can do some error checking
        # like if parent refuses to return a sequence
        
        #cache this?
        
        if self.parent.seq == None:
            raise TypeError, \
                  ("Oligo (%s) can't find its parent's (%s) sequence." %
                   (self.uid,self.parent.id))


        try:
            return self.parent.seq[self.start:self.start+self.length+1]
        except IndexError:
            if DEBUG:
                print ("""Oligo's indices overflow parent's sequence.
uid:\t%s
parent:\t%s
start:\t%s
length:\t%s
end:\t%s""" %
                       (self.uid,self.parent.id,
                        self.start,self.length,self.start+self.length+1))
            raise

        except TypeError:
            if DEBUG:
                print ("""Oligo's parent's sequence is not a string.
uid:\t%s
parent:\t%s""" %
                       (self.uid,self.parent.id))
            raise

            
    def score(self,sc_score):
        """Return an oligo's score.
        """
        sc = sc_score.scorecard
        try:
            pos = self.scIndexes[sc]
        except KeyError:
            if DEBUG:
                print "Oligo and score card not associated"
                return None
        except:
            raise

        return sc.scores[sc_score][pos]
        
        

    def rank(self,sc_score):
        """Return the rank of an oligos score.
        """
        pass

    def percentile(self,sc_score):
        """Return the percentile of an oligo's score.
        """
        pass

    def __str__ (self):
        """FASTA representation
        """
        return ">%s\n%s\n" % (self.uid,self.seq())



        
    
        
        
        
