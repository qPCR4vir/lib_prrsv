#!/usr/local/bin/python 
"""Usage: m8familyFilter.py <family NCBI TaxonID>

This program reads from stdin and writes to stdout.
"""
__version__ = tuple([int(ver) for ver in "$Revision: 1.2 $".split()[1].split('.')])
__author__ =  "Kael Fischer"

import sys
from sequence import blastNoSQL as B
from ncbi import giInfo
import utils

def main(args):
    try:
        famID=int(args[-1])
    except:
        print >> sys.stderr,  __doc__
        return 1
    famGis=utils.noneDict(giInfo.GiInfo._table.Gi.valueList(Family_Tax_ID=famID))

    if len(famGis) ==0:
        print >> sys.stderr, "No gis found for family taxID: %s" % famID 
        return 1

    for l in sys.stdin:
        try:
            gi = int(B.giFromM8name(l.split()[1]))
        except:
            continue
        if gi in famGis:
            print l,
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
