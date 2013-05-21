#!/usr/local/bin/python 
#
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.5 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

import sys
import optparse
from sequence import fasta
import utils

def main(sysargs=[]):

    oneLineUsage = "Usage: %prog [options] <input files>"

    # set a long description for the help message
    description="Expand sequences in input files so that all combinations of ambiguious positions are represented"

    op = optparse.OptionParser(
        oneLineUsage,description=description,
        version="%prog " + '.'.join([str(x) for x in __version__]))

    #
    # OPTION DEFINITIONS
    #

    #  Below are some examples of setting up options
    #  Note the above line gives you -h, --help, --version
    #  and an error message if an unspecified option is given.
    #  for more see: http://docs.python.org/library/optparse.html#defining-options
## 
##     op.add_option('-t','--tilesizes',dest="tileSizesL", default="60",
##                   help="comma seperated list of tile sizes to generate. e.g. '50,70'.")
##     op.add_option('-s','--shift',default=25,type='int',dest="shift",
##                   help="Start position offset" )
    op.add_option('-c','--clobber',action="store_true", default=False,
                  dest="clobber", help="allow output file to be overwritten.")
##     op.add_option('-E','--endtile',action="store_false", default=False,
##                   dest="endtile", help="Suppress generation of an tile that covers the"
##                   " last base in each recrod regardless of the -s and -t settings.")
    op.add_option('-o','--outputfile',dest="outFileName",default='-',
                  help="name of output file. '-' means: STDOUT")


    #
    # OPTION PATCHING
    #
    
    # all options should be defined above here

    # add defaults to help messages
    for o in op.option_list:
        if o.type == None:
            continue
        if o.help == None:
            o.help = "Default: %default"
        else:
            o.help+= " Default: %default"
    #
    # OPTION PARSING
    #

    # call the parser
    (opts,args) = op.parse_args(sysargs)

    # option values are available as attributes of the opts object
    # e.g. opts.myFavoriteDest


    #
    # OPTION CHECKING
    #

    #   Example of option checking beyond what parse_args does:
    #   (Note that parse_args can do some basic type checking
    #   and conversion.)
    
## 
##     try:
##         tileSizes = [int(x) for x in opts.tileSizesL.split(',') if len(x)>0]
##         outFile = utils.safeOFW(opts.outFileName,clobber=opts.clobber)
##     except Exception, eData:
##         print >> sys.stderr, ("\nOption Error: %s\n" %eData.message)
##         print >> sys.stderr, op.format_help()
##         return 1
## 

    #
    # ARGUMENT CHECKING
    #

    # you may like this simple way of requiring
    # >1 argument, assuming the one or more files needs
    # to be processed

    try:
        if len(args) == 0:
            raise RuntimeError, "No input files specified."
    except Exception, eData: 
        print >> sys.stderr, ("\nUsage Error: %s\n" %eData.message)
        print >> sys.stderr, op.format_help()
        return 1

##     #
##     # Here is a more complicated argument checking example
##     #
##     minArgs = 0    ## SET THAT
##     maxArgs = None ## SET THAT (can be None for no limit)
    
##     argProblem = False
##     if maxArgs != None and len(args) > maxArgs:
##         print >> sys.stderr, ("\nWrong number of files or arguments: "
##                               "%s found (expected max of %s)"
##                               % (len(args), maxArgs))
##         argProblem = True
        
##     elif len(args) < minArgs:
##         print >> sys.stderr, ("\nWrong number of files or arguments: "
##                               "%s found (expected min of %s)"
##                               % (len(args),minArgs))
##         argProblem = True

##     if argProblem:
##         print >> sys.stderr, op.format_help()
##         return 1


    oFile=utils.safeOFW(opts.outFileName,clobber=opts.clobber)
    
    somethingUsefull(args,oFile)
    
    return(0)  # we did it!

    
def somethingUsefull(args,oFile):
    for f in args:
        for rec in fasta.generalIterator(f):
            if fasta.isAllDRNA(rec.sequence):
                print >> oFile, rec
            else:
                n=0
                for uaSeq in fasta.expandAmbiguousSequence(rec.sequence):
                    uaTitle=rec.title+'|ua%s'%n
                    n+=1
                    print >> oFile, rec.__class__(title=uaTitle,sequence=uaSeq)
    

####### LEAVE THIS ALONE ###########
# If run directly as a program
# this calls the main function.
# 
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


