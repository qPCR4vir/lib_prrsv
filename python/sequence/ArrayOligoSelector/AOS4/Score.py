#!/usr/bin/env python
#
# Score
# ArrayOligoSelector
#
# oligo scoring classes
#


#
# 	$Id: Score.py,v 1.2 2006/03/02 06:41:43 kael Exp $	
#

import Oligo
import ScoreCard


# package level defaults
import __init__ as AOS
DEBUG = AOS.DEBUG

class Score:
  """Scores are initilized with their score card, if any.
  They are called with an oligo as the input and return a
  number, a bool or None if the oligo is flagged it the context
  of the scores ScoreCard.
  """

  def __init__ (self,scorecard=None):
    """
    """
    self.scorecard=scorecard

  def calcScore(self,oligo):
    """
    """
    return None
  
  def __call__(self,oligo):
    """
    """
    if not isinstance(oligo,Oligo.Oligo):
      raise TypeError, "oligo must be an Oligo instance"

    # if there is a ScoreCard
    # check for filtering
    if isinstance(self.scorecard,ScoreCard.ScoreCard):
      try:
        if oligo.scFlagged[self.scorecard]:
          return None
      except KeyError:
        if DEBUG:
          print "oligo is not associated with scorecard."
        raise

    # otherwise calculate score
    return self.calcScore(oligo)



class GCcontent(Score):
  """
  """
  def calcScore(self,oligo):

    seq = oligo.seq().upper()
    return float(seq.count('C')+seq.count('G'))/float(len(seq))
    


class SelfHyb_SW(Score):
  """
  """
  pass


class CrossHyb(Score):
  """
  """
  pass

class Complexity_LZW(Score):
  """
  """
  pass


    

