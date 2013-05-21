#!/usr/local/bin/python
#
# Removes duplicates from a fasta file and counts a fasta file and
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
    runCommand((FASTA_FILTER_PATH, '-Q', '-d', options.inputfile))
    output = '.'.join((os.path.splitext(options.inputfile)[0],
                      'noDups', 'fastq'))
    readsC, sequencesC = countFasta(options.outputfile)
    if readsC == 0 and sequencesC == 0:
        raise Exception
    rcType = ReadCountType(Name='noDups', insertIfNeeded=True)
    ReadCount(Sequence_File_ID=options.Sequence_File_ID,
              Reads=readsC, Sequences=sequencesC,
              Read_Count_Type_ID=rcType, insertIfNeeded=True)

if __name__ == '__main__':
    try:
        retCode = main(sys.argv[1:])
    except:
        retCode = 1
    sys.exit(retCode)
