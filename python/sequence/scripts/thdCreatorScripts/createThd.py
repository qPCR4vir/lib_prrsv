#!/usr/local/bin/python
#
# Runs the alignment command on the input filename and creates the THD
# from the result.
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.2 $'.split()[1].split('.')])
__author__ = "Julio Menendez"

import sys, os
from optparse import OptionParser

from viroinfo.megaseq import runCommand, THD, AlignmentProtocol

def main(args):
    parser = OptionParser('usage: %prog [options]')
    optionsAvailable = (
        ('-i', '--input', 'inputfile', 'string'),
        ('-s', '--sequencefileid', 'Sequence_File_ID', 'int'),
        ('-x', '--filterindex', 'filterindex', 'int'),
        ('-a', '--alignment', 'Alignment_Protocol_ID', 'int'),
        ('-t', '--taxonomy', 'Taxonomy_ID', 'int'),
        ('-e','--expect', 'Alignment_Score_Limit','float'),
        ('-f', '--filter', 'Blast_Filter_Protocol_ID', 'int'),
        )
    for oa in optionsAvailable:
        parser.add_option(oa[0], oa[1], dest=oa[2], type=oa[3])
    
    (options, args) = parser.parse_args(args)
    try:
        THD.fromMbrFile(options.inputfile, options.Alignment_Score_Limit,
                        options.Sequence_File_ID,
                        options.Alignment_Protocol_ID,
                        options.Taxonomy_ID,
                        filterProtocolID=options.Blast_Filter_Protocol_ID,
                        filterOutputIndex=options.filterindex)
    except Exception, ex:
        print >> sys.stderr, str(ex)
        import traceback
        print >> sys.stderr, traceback.format_exc()
        return 1
    return 0

if __name__ == '__main__':
    try:
        retCode = main(sys.argv[1:])
    except:
        retCode = 1
    print >> sys.stderr, retCode
    sys.exit(retCode)
