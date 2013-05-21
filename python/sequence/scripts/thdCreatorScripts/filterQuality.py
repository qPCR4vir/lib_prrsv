#!/usr/local/bin/python
#
# Filters a fasta file based on the quality and
# stores the value in the ReadCount table of MegaChipDB.
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.2 $'.split()[1].split('.')])
__author__ = "Julio Menendez"

import sys, os
from optparse import OptionParser

from viroinfo.megaseq import (runCommand, ReadCount, countFasta,
                              FASTA_FILTER_PATH, ReadCountType)

def main(args):
    parser = OptionParser('usage: %prog [options]')
    optionsAvailable = (
        ('-i', '--input', 'inputfile', 'string'),
        ('-o', '--output', 'outputfile', 'string'),
        ('-s', '--sequencefileid', 'Sequence_File_ID', 'int'),
        )
    for oa in optionsAvailable:
        parser.add_option(oa[0], oa[1], dest=oa[2], type=oa[3])

    (options, args) = parser.parse_args(args)

    outfn = options.outputfile.replace('.fasta','')
    out, rc = runCommand((FASTA_FILTER_PATH, '-r0.4', '-H18', '-a:4',
                          '-A0', '--output=%s' % outfn,
                          options.inputfile), True) 
    if rc != 0:
        raise Exception

    readsC, sequencesC = countFasta(options.outputfile)
    if int(readsC) == 0 and int(sequencesC) == 0:
        print readsC, sequencesC
        raise Exception
    rcType = ReadCountType(Name='qualityFiltered',
                           insertIfNeeded=True)
    ReadCount(Sequence_File_ID=options.Sequence_File_ID,
              Reads=readsC, Sequences=sequencesC,
              Read_Count_Type_ID=rcType, insertIfNeeded=True)

if __name__ == '__main__':
    try:
        retCode = main(sys.argv[1:])
    except:
        retCode = 1
    sys.exit(retCode)
