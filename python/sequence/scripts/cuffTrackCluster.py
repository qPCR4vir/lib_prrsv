#!/usr/local/bin/python 
#
# Fischer Lab new script template
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.1 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

import sys
import utils
import optparse

from sequence import cuff

def main(sysargs=[]):

    # short usage message
    oneLineUsage = "Usage: %prog [options] <input file>"

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
    op.add_option('-L','--no-log2-transform',action='store_false',
                  dest='log2transform',
                  default=True,help='turn off log transformation')

    # all options should be defined above here
    # (you don't have to define version or help)

    # add defaults to help messages     
    for o in op.option_list:            
        if o.type == None:              
            continue                    
        if o.help == None:              
            o.help = "Default: %default"
        else:                           
            o.help+= " Default: %default"
                                       
    # OPTION PARSING                    
    (opts,args) = op.parse_args(sysargs)

    # OPTION CHECKING
    try:
        # (check an opt) tileSizes = [int(x) for x in opts.tileSizesL.split(',') if len(x)>0]
        outFile = utils.safeOFW(opts.outFileName,
                                clobber=opts.clobber) # open the output file
    except Exception, eData:
        print >> sys.stderr, ("\nOption Error: %s\n" %eData.message)
        print >> sys.stderr, op.format_help()
        return 1

    # ARGUMENT CHECKING
    # something like
    try:
        if len(args) != 1:
            raise RuntimeError, "Must specify one input file"
        else:
            inFile=args[0]
            
    except Exception, eData: 
        print >> sys.stderr, ("\nUsage Error: %s\n" %eData.message)
        print >> sys.stderr, op.format_help()
        return 1

    # now do some work.

    return (doIt(inFile,outFile,opts.log2transform))

#
# end of main
#

    
def doIt(inFile,outFile,log2transform):
    """This is where the magic happens
    """
    try:
        gt=cuff.Gene_Tracking(inFile)
        gt.cluster(outFile,log2transform)
        return(0)
    except Exception, eData: 
        print >> sys.stderr, ("\nError processing %s:\n%s\n" %
                              (inFile,eData.message))
        return 1
    

####### LEAVE THIS ALONE ###########
# If run directly as a program
# this calls the main function.
# 
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


