#!/usr/local/bin/python 
#
# 
#
"""Usage: insertGridItBlastResults.py [options] <tab delimited blast results>

Synopsis:
    Insert blast -m8 or megablast -D3 results in many files into the
    Blast_Results db.

    It may be best to parallelize your job by using gridIt to run several
    of thsese at the same time.

    '-' may be specified as the input to specify STDIN.
"""

__version__ = tuple([int(x) for x in
                     '$Revision: 1.2 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

import os
import subprocess
import sys
import optparse

import utils

def main(sysargs):

    usageStr =  __doc__
    epilog = ""

    op = optparse.OptionParser(usageStr,
                         version=__version__)

    #
    # OPTION DEFINITIONS
    #
    op.add_option('-t','--tablesuffix', default="",
                  help="Hsp table suffix")
    op.add_option('-c','--clobber',action="store_true", default=False,
                  help="allow output file to be overwritten.")
    op.add_option('-l','--logfile',default='-',
                  help="name of log file. '-' means: STDOUT")
 
    #
    # OPTION PATCHING
    #
    # add "defaults" to help messages
    for o in op.option_list:
        if o.type != None or o.action.startswith('store_'):
            if o.help == None:
                o.help = "Default: %default"
            else:
                o.help+= " Default: %default"

    # call the parser
    (opts,m8files) = op.parse_args(sysargs)

    #
    # OPTION CHECKING
    #

    #   Example of option checking beyond what parse_args does:
    #   (Note that parse_args can do some basic type checking
    #   and conversion.)
    try:
        outFile = utils.safeOFW(opts.logfile,clobber=opts.clobber)
    except Exception, eData:
        print >> sys.stderr, ("\nOption Error: %s\n" %eData.message)
        print >> sys.stderr, op.format_help()
        return 1


    #
    # ARGUMENT CHECKING
    #

    # you may like this simple way of requiring
    # >1 argument, assuming the one or more files needs
    # to be processed

    try:
        if len(m8files) == 0:
            raise RuntimeError, "No input files specified."
    except Exception, eData: 
        print >> sys.stderr, ("\nUsage Error: %s\n" %eData.message)
        print >> sys.stderr, op.format_help()
        return 1


    from sequence import blast
    mySearch = blast.Search(Hsp_Table_Suffix=opts.tablesuffix,
                            insertIfNeeded=True)
    myHsps = mySearch.hspClass()

    myHsps._table.disableKeys()
    try:
        myHsps.insertM8records(m8files)
    finally:
        myHsps._table.enableKeys()
    
    return(0)  # we did it!

  
####### LEAVE THIS ALONE ###########
# If run directly as a program
# this calls the main function.
# 
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


