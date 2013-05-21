#!/usr/local/bin/python
#
# Creates THD from sequening reads.
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.2 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

import sys
import re
from optparse import OptionParser

from sequence import blastNoSQL as B

def main(args):
    description = ('Calculates dG for alignmnets in results in '
                   'mbr file add floating point energy values to '
                   'the end of each m8record. Results are saved '
                   'to <mbr file>.dG')
    parser = OptionParser(('usage: %prog [options] '
                           '<mbr file> <query path> '
                           '<database path>'),
                          description=description)
    (options, args) = parser.parse_args(args)

    if len(args) < 3:
        print >> sys.stderr, ('Usage error: all arguments are '
                              'required')
        print >> sys.stderr, parser.format_help()
        return 1
    
    inMBR = args[0]
    qPath = args[1]
    dbPath = args[2]

    oPath = inMBR + '.dG'
    
    B.m8alignmentEnergies(inMBR,qPath,dbPath,oPath)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
