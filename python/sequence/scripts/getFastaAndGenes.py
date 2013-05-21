#!/usr/bin/env python

import sys
import re
from optparse import OptionParser

from ncbi import giInfo, genbank
from utils import safeOFW
from sequence import fasta


def main(args):
    parser = OptionParser(('usage: %prog [options] '
                          '<qualifier to use as gene name>'))
    parser.add_option('-g', '--gi', dest='gi', type='int',
                      help='gi of the genbank record to analyze.')
    parser.add_option('-r', '--record', dest='record', type='string',
                      help='read genbank record from FILE',
                      metavar='FILE')
    parser.add_option('-o', '--output', dest="output",
                      help='name of output file.')
    parser.add_option('-c', '--clobber', dest="clobber",
                      help='Clobber the output files.',
                      action='store_true', default=False)
    (options, args) = parser.parse_args(args)

    if not options.gi and not options.record:
        print >> sys.stderr, 'Usage error: a record file (-r) or ' + \
              'a gi (-g) is required'
        print >> sys.stderr, parser.format_help()
        sys.exit(1)

    if options.output:
        filenameOut = options.output
    else:
        if options.record:
            filenameOut = '%s' % options.record
        else:
            filenameOut = '%s' % options.gi

    if len(args) < 2:
        print >> sys.stderr, 'Qualifier name is required'
        print >> sys.stderr, parser.format_help()
        sys.exit(1)

    if options.record:
        record = genbank.Record(file(options.record))
    else:
        record = giInfo.GiRecord(options.gi, True)
        
    qualifierName = args[1]

    fastaFh = safeOFW('%s.fasta' % filenameOut,
                      clobber=options.clobber)
    print >> fastaFh, record.fasta()
    fastaFh.close()

    genesFh = safeOFW('%s.genes.fasta' % filenameOut,
                      clobber=options.clobber)
    genes = []
    for feature in record.features():
        if len(feature.location.regions) != 1:
            continue

        if feature.type != 'gene':
            continue

        region = feature.location.regions[0]
        seq = record.sequence[region.start - 1:region.end]
        seq = seq.upper()
        if region.complement:
            seq = genbank.reverseComplement(seq)
        
        if feature.qualifiers.has_key(qualifierName):
            feat_name = feature.qualifiers[qualifierName]
        else:
            print >> sys.stderr, ("Feature doesn't have qualifier "
                                  "'%s'.\nQualifiers: %s" %
                                  (qualifierName, feature.qualifiers))
            continue

        print >> genesFh, fasta.Record(title=feat_name, sequence=seq)
        print >> genesFh, ''
    
    genesFh.close()


if __name__ == '__main__':
    main(sys.argv)
