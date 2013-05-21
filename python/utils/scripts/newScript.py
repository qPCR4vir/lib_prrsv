#!/usr/local/bin/python 
#
# Fischer Lab new script template
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.7 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

import sys
import utils
import optparse

def main(sysargs=[]):

    # short usage message
    oneLineUsage = "Usage: %prog [options] <input files>"

    # set a long desccription for the help message
    description=None

    op = optparse.OptionParser(
        oneLineUsage,description=description,
        version="%prog " + '.'.join([str(x) for x in __version__]))

    # OPTION DEFINITIONS
    op.add_option('-o','--outputfile',dest="outFileName",default='-',
                  help="name of output file. '-' means: STDOUT")
    op.add_option('-c','--clobber',action="store_true", default=False,
                  dest="clobber", help="allow output file to be overwritten.")

    #  Below are some examples of setting up options
    #  Note the above line gives you -h, --help, --version
    #  and an error message if an unspecified option is given.
    #  for more see: http://docs.python.org/library/optparse.html#defining-options
## 
##     op.add_option('-t','--tilesizes',dest="tileSizesL", default="60",
##                   help="comma seperated list of tile sizes to generate. e.g. '50,70'.")
##     op.add_option('-s','--shift',default=25,type='int',dest="shift",
##                   help="Start position offset" )
##     op.add_option('-E','--endtile',action="store_false", default=False,
##                   dest="endtile", help="Suppress generation of an tile that covers the"
##                   " last base in each recrod regardless of the -s and -t settings.")
##

    # all options should be defined above here

    ### Don't Change In Box ###################
    #                                         #
    # OPTION PATCHING                         #
    # add defaults to help messages           # 
    #                                         #
    for o in op.option_list:                  #
        if o.type == None:                    #
            continue                          #
        if o.help == None:                    #
            o.help = "Default: %default"      #
        else:                                 #
            o.help+= " Default: %default"     #
                                              #
                                              #
    #                                         #
    # OPTION PARSING                          #
    (opts,args) = op.parse_args(sysargs)      #
    #                                         #
    # option values are available as attrs    #
    # of the opts object, e.g. opts.myClobber #
    #                                         #
    # args is a list of the balance of the    #
    # command line                            # 
    ###########################################


    # OPTION CHECKING
    #   Option checking beyond what parse_args does:
    #   (Note that parse_args can do some basic 
    #   type checking and conversion.)    
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

    somethingUsefull(args,outFile=outFile)
    
    return(0)  # we did it!

#
# end of main
#

    
def somethingUsefull(args,outFile=sys.stdout):
    """This is where the magic happens
    """
    print >> outFile, "now doing an important thing"
    somethingImportant = 1
    return somethingImportant

####### LEAVE THIS ALONE ###########
# If run directly as a program
# this calls the main function.
# 
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


