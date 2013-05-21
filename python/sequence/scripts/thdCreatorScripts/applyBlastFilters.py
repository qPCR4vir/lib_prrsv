#!/usr/local/bin/python
#
# Apply a blast filter to a fasta file and
# stores the count values of each output file
# in the ReadCount table of MegaChipDB.
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.2 $'.split()[1].split('.')])
__author__ = "Julio Menendez"

import sys, os
from optparse import OptionParser

from viroinfo.megaseq import (ReadCount, countFasta,
                              BlastFilterProtocol,
                              applyFilterProtocol,
                              ReadCountType)

def main(args):
    parser = OptionParser('usage: %prog [options]')
    optionsAvailable = (
        ('-i', '--input', 'inputfile', 'string'),
        ('-f', '--filterid', 'Blast_Filter_Protocol_ID', 'int'),
        ('-s', '--sequencefileid', 'Sequence_File_ID', 'int'),
        )
    for oa in optionsAvailable:
        parser.add_option(oa[0], oa[1], dest=oa[2], type=oa[3])
    
    (options, args) = parser.parse_args(args)
    outputs = applyFilterProtocol(options.inputfile,
                                  options.Blast_Filter_Protocol_ID)
    protocol = BlastFilterProtocol(options.Blast_Filter_Protocol_ID)
    stepsC = len(list(protocol.steps()))
    for fn in outputs:
        parts = fn.split('.')[-(stepsC+1):][:-1]
        rcType = ReadCountType(Name=''.join(parts),
                               insertIfNeeded=True)
        readsC, sequencesC = countFasta(fn)
        if readsC == 0 and sequencesC == 0:
            raise Exception
        ReadCount(Sequence_File_ID=options.Sequence_File_ID,
                  Reads=readsC,
                  Sequences=sequencesC,
                  Read_Count_Type_ID=rcType,
                  insertIfNeeded=True)
        
if __name__ == '__main__':
    try:
        retCode = main(sys.argv[1:])
    except:
        retCode = 1
    sys.exit(retCode)
