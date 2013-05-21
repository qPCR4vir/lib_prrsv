#!/usr/local/bin/python
#
# Generates a report with the energies for each query from a D1 blast
# results file.
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.1 $'.split()[1].split('.')])
__author__ = "Julio C. Menendez"

import sys
import re
from optparse import OptionParser

from sequence import blastNoSQL as B

def main(args):
    parser = OptionParser(('usage: %prog [options] '
                           '<D1 file> <query path> '
                           '<database path>'),)
    parser.add_option('-o','--output',dest='outFileName',
                  help='name of output file')
    (options, args) = parser.parse_args(args)

    if len(args) < 3:
        print >> sys.stderr, ('Usage error: all arguments are '
                              'required')
        print >> sys.stderr, parser.format_help()
        return 1
    
    inD1 = args[0]
    qPath = args[1]
    dbPath = args[2]

    if not options.outFileName:
        oPath = inD1 + '.dG'
    else:
        oPath = options.outFileName
    
    B.d1alignmentEnergiesReport(inD1, qPath, dbPath, oPath)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
