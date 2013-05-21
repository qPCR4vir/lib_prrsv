#!/usr/local/bin/python 
#
# Fischer Lab new script template
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.2 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

import sys
import utils
import optparse

from sequence.fasta import removeDuplicates

def main(sysargs=[]):

    # short usage message
    oneLineUsage = "Usage: %prog [options] <input files>"

    # set a long desccription for the help message
    description=None

    op = optparse.OptionParser(
        oneLineUsage,description=description,
        version="%prog " + '.'.join([str(x) for x in __version__]))

    # OPTION DEFINITIONS
    op.add_option('-d','--debug',action="store_true", default=False,
                  dest="debug", help="print progress information on stderr.")


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

    # now do some work.

    removeDuplicates(args,debug=opts.debug)
    
    return(0)  # we did it! 

#
# end of main
#

####### LEAVE THIS ALONE ###########
# If run directly as a program
# this calls the main function.
# 
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))



