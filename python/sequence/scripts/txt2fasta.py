#!/usr/local/bin/python 
#
# Fischer Lab new script template
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.6 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

import sys
import types
import utils
import optparse

import os.path

from sequence import fasta

defaultFastaWidth=50
defaultOutName='-'
defaultSpaceChar=' '

def main(sysargs=[]):

    # short usage message
    oneLineUsage = "Usage: %prog [options] <sequence files>"

    # set a long desccription for the help message
    description=(
        "Takes files containing reads (raw sequence with no title info)"
        " and combines them into a fasta file.  The fasta titles are"
        " derived from the file names."
        )

    op = optparse.OptionParser(
        oneLineUsage,description=description,
        version="%prog " + '.'.join([str(x) for x in __version__]))

    # OPTION DEFINITIONS
    op.add_option('-o','--outputfile',dest="outFileName",default=defaultOutName,
                  help="name of output file. '-' means: STDOUT")
    op.add_option('-c','--clobber',action="store_true", default=False,
                  dest="clobber", help="allow output file to be overwritten.")

    op.add_option('-w','--width',dest="fmtWidth", default=defaultFastaWidth,
                  help="formatted width of output fasta records.")

    op.add_option('-s','--spaceChar',dest="spaceChar",default=defaultSpaceChar,
                 help="replace title spaces with this character/string")

    op.add_option('-P','--useFullPaths',action="store_true",
                  dest="useFullPaths",default=False,
                  help="include pathnames in titles")

    # all options should be defined above here

    # OPTION PATCHING
    # add defaults to help messages
    for o in op.option_list:
        if o.type == None:
            continue
        if o.help == None:
            o.help = "Default: %default"
        else:
            o.help+= " Default: %default"

    # OPTION PARSING
    (opts,args) = op.parse_args(sysargs)  # call the parser

    # option values are available as attributes of the opts object
    # e.g. opts.myClobber
    # args is a list of the balance of the command line

    # OPTION CHECKING
    #   Example of option checking beyond what parse_args does:
    #   (Note that parse_args can do some basic type checking
    #   and conversion.)
    

    try:
        # (check an opt) tileSizes = [int(x) for x in opts.tileSizesL.split(',') if len(x)>0]
        outFile = utils.safeOFW(opts.outFileName,clobber=opts.clobber) # open the output file
    except Exception, eData:
        print >> sys.stderr, ("\nOption Error: %s\n" %eData.message)
        print >> sys.stderr, op.format_help()
        return 1

    # ARGUMENT CHECKING
    # these are the non-option command line things
    #
    #   If this is the command line:
    #   myScript.py -o myOutput.txt myInput0.txt myInput1.txt
    #
    #   opts.outputFileName = 'myOutput.txt'
    #   args = ['myInput0.txt','myInput1.txt']

    
    # you may like this simple way of requiring
    # 1 or more arguments, assuming that one or more
    # (e.g.) files need to be processed
    try:
        if len(args) == 0:
            raise RuntimeError, "No input files specified."
    except Exception, eData: 
        print >> sys.stderr, ("\nUsage Error: %s\n" %eData.message)
        print >> sys.stderr, op.format_help()
        return 1


    # Here is a more complicated argument checking case
    #
    # probably just want to use the simple example (above)
    # or something like this (below).
    #
    # In this example they are funtionaly the same
    # but in the #args = 0 case the error is already handled
    #
    argProblem = False
    minArgs = 1    ## SET THAT
    maxArgs = None ## SET THAT (can be None for no limit)
    
    if maxArgs != None and len(args) > maxArgs:
        print >> sys.stderr, ("\nWrong number of files or arguments: "
                              "%s found (expected max of %s)"
                              % (len(args), maxArgs))
        argProblem = True
    elif len(args) < minArgs:
        print >> sys.stderr, ("\nWrong number of files or arguments: "
                              "%s found (expected min of %s)"
                              % (len(args),minArgs))
        argProblem = True

    if argProblem:
        print >> sys.stderr, op.format_help()
        return 1

    # basic command line interface errors and usage cases
    # have been handled.

    # now do some work.

    for fn in args:
        print >> outFile, seqFile2FastaRecord(fn,
                                              fmtWidth=opts.fmtWidth,
                                              spaceChar=opts.spaceChar,
                                              useFullPaths=opts.useFullPaths)
        print >> outFile
        

    return(0)  # we did it!

#
# end of main
#

    

def seqFile2FastaRecord(inFile,fmtWidth=defaultFastaWidth,
                        spaceChar=defaultSpaceChar,
                        useFullPaths=False):
    """
    """
    if type(inFile) in types.StringTypes:
        fn = inFile
        inFile=file(fn)
    else:
        fn = inFile.name

    rv=fasta.Record()
    pathRt,pathTail=os.path.split(fn)

    if useFullPaths:
        rv.title=(os.path.splitext(fn)[0]).replace(' ',spaceChar)
    else:
        rv.title=(
            os.path.splitext(os.path.split(fn)[-1])[0]).replace(
            ' ',spaceChar)
    
    rv.sequence=fasta.dropNonLetters(inFile.read())
    return rv
    


####### LEAVE THIS ALONE ###########
# If run directly as a program
# this calls the main function.
# 
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


