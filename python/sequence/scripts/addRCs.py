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

from sequence import fasta

def main(sysargs=[]):

    # short usage message
    oneLineUsage = "Usage: %prog [options] <input files>"

    # set a long desccription for the help message
    description=("Add reverse complement fasta records to Fasta file(s). "
                 "Use '-' to read from STDIN.")

    op = optparse.OptionParser(
        oneLineUsage,description=description,
        version="%prog " + '.'.join([str(x) for x in __version__]))

    # OPTION DEFINITIONS
    op.add_option('-o','--outputfile',dest="outFileName",default='-',
                  help="name of output file. '-' means: STDOUT")
    op.add_option('-c','--clobber',action="store_true", default=False,
                  dest="clobber", help="allow output file to be overwritten.")


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
    try:
        # (check an opt) tileSizes = [int(x) for x in opts.tileSizesL.split(',') if len(x)>0]
        outFile = utils.safeOFW(opts.outFileName,clobber=opts.clobber) # open the output file
    except Exception, eData:
        print >> sys.stderr, ("\nOption Error: %s\n" %eData.message)
        print >> sys.stderr, op.format_help()
        return 1

    # ARGUMENT CHECKING
    # these are the non-option command line things
    try:
        if len(args) == 0:
            raise RuntimeError, "No input files specified."
    except Exception, eData: 
        print >> sys.stderr, ("\nUsage Error: %s\n" %eData.message)
        print >> sys.stderr, op.format_help()
        return 1


    processFile(args,outFile=outFile)
    
    return(0)  # we did it!

#
# end of main
#

def rcString(fRec):
    """
    """
    return '\n'.join((fRec.title+'_rc',fRec.reverseComplement()))

    
def processFile(args,outFile=sys.stdout):
    """This is where the magic happens
    """
    for f in args:
        for rec in fasta.FastaIterator(f):
            print >> outFile, rec
            print >> outFile, rcString(rec)
    return True

####### LEAVE THIS ALONE ###########
# If run directly as a program
# this calls the main function.
# 
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


