#!/usr/local/bin/python
#
# Redundant and Human sequence screening
#
#
__version__ = tuple([int(ver) for ver in
                     "$Revision: 1.8 $".split()[1].split('.')])

__author__ = "Kael Fischer"

from grid import GridThreads
import grid
import re
import sys
import os
import os.path
import string
import tempfile
import threading
import time
import commands
from types import *

from copy import copy

import fasta 

from utils import timeprofile
timer = timeprofile()

FastaIterator = fasta.FastaIterator


HERVtitles = [re.compile('endogenous.*retro',re.IGNORECASE),
              re.compile('HERV',re.IGNORECASE)]

HIVtitles = [re.compile('HIV-[12][,) ]',re.IGNORECASE),
             re.compile('Human immunodeficiency',re.IGNORECASE)]




class MBGridScreen(GridThreads.GridThreadManager):
    """
    """
    def __init__ (self, dbName, qryFileNames,
                  MBparams="-W24 -D3 -E10 -G24 -FF"):


        GridThreads.GridThreadManager.__init__(self)

        self.dbName = dbName
        self.qryFileNames = copy(qryFileNames)
        self.MBparams = MBparams
        self.doneQryFiles = []
        self.errorQryFiles = []
        self.started=False
        self.status=[]
        self.output=[]


    def submitJobs(self):
        """
        """
        for f in self.qryFileNames:
            f = os.path.split(f)[1]
            mbCmd = ('megablast %s -d %s -f -R -i %s -o %s'
                     % (self.MBparams,self.dbName,f,f+'.MBout'))
            screenCmd = ('./HSPscreen.py %s %s '  % (f+'.MBout',f))
            
            self.submitThread('%s' , [mbCmd,screenCmd],
                              qrshArgs='-cwd -l arch=fbsd-amd64')
                     



class IterativeScreen (object):
    """framework for iteratively screening through a fa
    """

    def __init__ (self,rpycServer):
        """
        """
        self.server = rpycServer


    def setup(self,inFile):
        """
        """
        

    def screen(self,screenFasta):
        """screen and return the longest sequence after screening
        """
        pass
    
