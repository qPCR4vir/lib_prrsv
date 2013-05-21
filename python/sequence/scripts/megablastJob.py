#!/usr/local/bin/python 
#
# Runs megablast and reads the output file to check the tail
# corresponding to a successful completion of megablast.
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.1 $'.split()[1].split('.')])
__author__ = "Julio Carlos Menendez"

import sys
import utils
import optparse
import commands
import re


def main(sysargs=[]):
    args =  ' '.join(sysargs)
    cmd = 'megablast %s' % args
    status, output = commands.getstatusoutput(cmd)
    if len(output) > 0:
        print >> sys.stdout, output

    # Get the input and output filenames from the args.
    inputFilename = None
    outputFilename = None
    regexp = r'-%s\s*(\S+)'
    match = re.search(regexp % 'i', args)
    if match:
        inputFilename = match.group(1)

    match = re.search(regexp % 'o', args)
    if match:
        outputFilename = match.group(1)
        
    if inputFilename and outputFilename:
        # Check megablast reported the same number of entries existing
        # in the input file.
        content = open(inputFilename).read()
        nEntries = len(re.findall(r'^>', content, re.M))
        tailOut = commands.getoutput('tail -n1 %s' % outputFilename)
        match = re.match(r'^#(?:.*)finished(?:.*) (\d+) (?:.*)',
                         tailOut)
        if match is None or int(match.group(1)) != nEntries:
            return 1
    return status

####### LEAVE THIS ALONE ###########
# If run directly as a program
# this calls the main function.
# 
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


