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
import gzip

def main(sysargs=[]):

    # short usage message
    oneLineUsage = "Usage: %prog [options] <input files>"

    description="""
Checks for binary zeros in files. Exit status is non zero if zeros
are detected. Checking a file is stopped as soon as a zero is found.
Gzip'ed files are accepted, and unzipped on the fly. '-' indicates
that standard input should be checked. 
"""

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
        # (check an opt) 
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

    # now do some work.

    errCt=0
    for fileName in args:
        errCt+=checkFile(fileName,outFile=outFile,verbose=True)
        
    return(0)  # we did it!

#
# end of main
#

def checkFile(fileName,outFile=sys.stderr,verbose=False):
    """returns True if binary zero character is found in file.
    If verbose is set to True, an error message is sent to
    the output file.
    """
    if fileName == '-':
        f=sys.stdin
        fileName='/dev/stdin'
    elif fileName.endswith('.gz'):
        f=gzip.open(fileName)
    else:
        f = file(fileName)
    n=0
    for l in f:
        n+=1
        if l.count('\000') > 0:
            if verbose:
                #print f.tell()
                print 'Red Alert: catastrophuck on %s' %  fileName
                print 'First offending line: %d' % n
            f.close()
            return True
    if verbose:
        print 'Condition Green: all is well with %s' %  fileName
    f.close()
    return False

####### LEAVE THIS ALONE ###########
# If run directly as a program
# this calls the main function.
# 
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


