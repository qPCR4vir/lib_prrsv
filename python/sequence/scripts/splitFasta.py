#!/usr/local/bin/python 
#
# splitFasta.py
#
# Take a fasta file and split it up in to n sections
#

"""USAGE
        splitFasta.py <[-h?] | -r <max records> <fasta file> | [<count> <fasta file>]> 

SUMMARY
        Take a fasta file and split it up into count sections.


        
OPTIONS
        -h      print help message
        -?
"""
__version__ = "$Revision: 1.5 $".split(':')[-1].strip(' $') # don't edit this line.
                                                            # cvs will set it.

__author__ = "Kael Fischer"


import os,sys
import datetime
import shutil
from math import ceil, log10
import os.path
from sequence import fasta

# modules for command line handling
import sys
import getopt


def main(args=None):
    shortOptions="h?r:"  
    longOptions = []   
    try:
        opts, args = getopt.getopt(args, shortOptions, longOptions)
    except getopt.error, msg:
        # there is an unknown option!
            print msg      # prints the option error
            print __doc__  # prints the usage message from the top
            return (-2)

    maxRec=None

    # process options
    for option,optionArg in opts:
        if option=='-h' or option=='-?':
            print __doc__
            return(0)     # '0' = no error in UNIX

        elif option == '-r':
            maxRec=int(optionArg)
        else:
            print "%s option not implemented" % option

    # check arguments
    # correct #, etc
    # the remaining command line arguments (after the
    # option processing) are in 'args'

    minArgs = 1 ## SET THAT
    maxArgs = 2 ## SET THAT
    #print args
    argProblem = False
    if len(args) > maxArgs:
        print ("Wrong number of arguments: %s found (expected max of %s)"
               % (lan(args),maxArgs))
        argProblem = True
    elif len(args) < minArgs:
        print ("Wrong number of arguments: %s found (expected min of %s)"
               % (len(args),minArgs))
        argProblem = True
    # put in other argument checks here
    # print help set argProbem if there is a problem

    ## PUT MORE ARGUMENT CHECKING HERE

    if argProblem:
        print __doc__
        return(-1)

    if maxRec != None:
        fn,=args
        recCt=fasta.fastaCount(fn)
        count = recCt/maxRec + int(recCt%maxRec!=0)
    else:
        count,fn= args
        count = int(count)

    nG = niceName(fn,count)

    fasta.splitFasta(fn,splitCt=count,nameGenerator=nG)
        
    return(0)  # we did it!

def niceName(rootName,maxN):
    
    n=0
    fmt = '_%%0%dd' % int(ceil(log10(maxN)))
    base,ext = os.path.splitext(rootName)
    while True:
        name=base+ fmt% (n) + ext
        yield file(name,'w'),name
        n+=1

    

####### LEAVE THIS ALONE ###########
# If run directly as a program
# this block calls the main function.
# 
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

# end program - beep!
