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

from sequence.fastq import phdQualIterator
from sequence.fasta import Record

def main(sysargs=[]):

    # short usage message
    oneLineUsage = "Usage: %prog [options] <Phred Score> <input fasta> <input quality>"

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
    op.add_option('-a', '--append', action="store_true", default=False,
                  dest="append", help="append to output file..")

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
                                append=opts.append,clobber=opts.clobber) # open the output file
    except Exception, eData:
        print >> sys.stderr, ("\nOption Error: %s\n" %eData.message)
        print >> sys.stderr, op.format_help()
        return 1

    # ARGUMENT CHECKING
    # something like
    try:
        if len(args) == 0:
            raise RuntimeError, "No input files specified."
    except Exception, eData: 
        print >> sys.stderr, ("\nUsage Error: %s\n" %eData.message)
        print >> sys.stderr, op.format_help()
        return 1

    # or look for other problems
    argProblem = False
    if 1 == -1:
        # look for your favorite problems here
        argProblem = True

    if argProblem:
        print >> sys.stderr, op.format_help()
        return 1
    print >> sys.stderr, args
    # now do some work.
    qScore,inFa,inQ=args
    
    somethingUsefull(qScore,inFa,inQ,outFile=outFile)
    
    return(0)  # we did it!

#
# end of main
#

    
def somethingUsefull(qScore,inFa,inQ,outFile=sys.stdout):
    """This is where the magic happens
    """
    for rec in phdQualIterator(inFa,inQ):
        rec.lowQ2N(qScore)
        #print >> outFile, rec.fasta()
        print >> outFile, Record(title=rec.title,sequence=rec.longestStretch('ACGTU'))


####### LEAVE THIS ALONE ###########
# If run directly as a program
# this calls the main function.
# 
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


