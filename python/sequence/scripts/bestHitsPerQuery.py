#!/usr/local/bin/python

from sequence import blastNoSQL
import utils

import sys
import os
import optparse

__version__ = tuple([int(x) for x in
                     '$Revision: 1.2 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

def main(sysargs):

    oneLineUsage = "Usage: %prog [-h]|[input files]"
    description = """Considering all input files (and/or stdin if '-' is specified)
as a single group of m8 results.  return only those results for each query that
have the highest score observed for that query.  All alignments
to that query sequence with that score are returned."""

    op = optparse.OptionParser(
        oneLineUsage,description=description,
        version="%prog " + '.'.join([str(x) for x in __version__]))

    try:
        (opts,args) = op.parse_args(sysargs)
        if len(args) == 0:
            raise RuntimeError, "No input files specified."
    except Exception, eData: 
        print >> sys.stderr, ("\nUsage Error: %s\n" %eData.message)
        print >> sys.stderr, op.format_help()
        return 1

    files = utils.readArgsOrFiles(args)
    extractHSPs(files)


def extractHSPs(inFile):
    rv={}
    for m8rec in blastNoSQL.m8generator(inFile):
        if m8rec['query'] not in rv:
            rv[m8rec['query']]=[m8rec]
        elif  m8rec['score'] > rv[m8rec['query']][0]['score']:
            rv[m8rec['query']]=[m8rec]
        elif m8rec['score'] == rv[m8rec['query']][0]['score']:
            rv[m8rec['query']].append(m8rec)

    kList = rv.keys()
    kList.sort()
    for k in kList:
        for m8 in rv[k]:
            print m8['_str_']
    



if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


    

    
