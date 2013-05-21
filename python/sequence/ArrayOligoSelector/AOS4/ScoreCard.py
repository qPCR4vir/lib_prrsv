#!/usr/bin/env python
#
# ScoreCard
# ArrayOligoSelector
#
# This is the workhose class for scoring of oligos,
# keeping track of the various scoring results.
#
# Decisons about optimal order of scoring functions may
# also be included here.[?]
#


#
# 	$Id: ScoreCard.py,v 1.3 2006/03/02 06:41:43 kael Exp $	
#

import Oligo
import Score

# package level defaults
import __init__ as AOS
DEBUG = AOS.DEBUG


class ScoreCard:
    """ScoreCard instances must know about oligos and background sequences.
    The scores themselves are stored in oligo-order in list-like arrays
    which are the values in the 'scores' dictionary.  This dictionary is keyed
    on the scoring instance.
    """

    def __init__ (self, oligos):
        """ Make a new score card for a list of oligos.
        """
        self.oligos = oligos
        self.scores={}
        self.background={}

        for i in range(len(self.oligos)):
            oligos[i].scIndexes[self] = i
            oligos[i].scFlagged[self] = False

    def addScore(self,scoreClass):
        """Add a new score vector.
        Overwrites old vector if any.
        """

        try:
            score = scoreClass(self)
            if not isinstance(score,Score.Score):
                raise TypeError
        except:
            raise TypeError, "scoreClass must be a Score or subclass"
        
        # these should probably be Numeric-like arrays 

        self.scores[score] = map(score,self.oligos)
        if DEBUG:
            print "score (%s) added to ScoreCard" % (score.__class__.__name__)
        
        

    def updatesScore(self,score,oligoIndexes):
        """Update score vector at particular oligo positions.
        oligoIndexes is a 1D bool. vector - 1 = update; 2 = ignore.
        """

        # the idea is that recently unflagged positions
        # can be added to the score vectors.
        pass

    def setBackground(self):
        pass
    

        
