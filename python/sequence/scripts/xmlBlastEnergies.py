#!/usr/local/bin/python 
#
# Fischer Lab new script template
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.1 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

import sys
from utils import safeOFW
import optparse

from sequence import blastNoSQL as B

def main(sysargs=[]):

    # short usage message
    oneLineUsage = "Usage: %prog [options] <input files>"

    # set a long desccription for the help message
    description=(
"""Calculate free energies from BLAST hits. Hits should be in ncbi XML format
(-m7 (blastall) -D2 -m7 (megablast)). Reports are put in <filename>.dG.
Multiple files can be specified, stdin can not be used to provide the hits.

Output format is: <query> <subject> <energy> 

""")
        

    op = optparse.OptionParser(
        oneLineUsage,description=description,
        version="%prog " + '.'.join([str(x) for x in __version__]))

    # OPTION DEFINITIONS
    op.add_option('-c','--clobber',action="store_true", default=False,
                  dest="clobber", help="allow existing output .dG files "
                  "to be overwritten.")


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

    return(processFiles(args,clobber=opts.clobber))
    
#
# end of main
#

def processFiles(files,clobber=False):
    for f in files:
        outF=safeOFW(f+'.dG',clobber=clobber)
        for hsp in B.xmlBlast2energy(f):
            print >>outF,'\t'.join([str(x) for x in [hsp[0],hsp[1],hsp[-1]]])
    return 0



####### LEAVE THIS ALONE ###########
# If run directly as a program
# this calls the main function.
# 
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


