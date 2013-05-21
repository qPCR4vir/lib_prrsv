#!/usr/local/bin/python
#
# Creates THD from sequening reads.
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.5 $'.split()[1].split('.')])
__author__ = "Julio Menendez, Kael Fischer"

import sys
from optparse import OptionParser

from viroinfo.megaseq import THD,THDCreator

def main(args):
    parser = OptionParser('usage: %prog [options]')
    optionsAvailable = (
        ('-i', '--input', 'Sequence_File_ID', 'int',
         'Sequence_File_ID from MegaChip.Sequence_File table.'),
        ('-a', '--alignment', 'Alignment_Protocol_ID', 'int',
         ('Alignment_Protocol_ID from MegaChip.Alignment_Protocol '
          'table.')),
        ('-t', '--taxonomy', 'Taxonomy_ID', 'int',
         'Taxonomy_ID from MegaChipDB.Taxonomy table.'),
        ('-e','--expect','Alignment_Score_Limit','float',
         'Score cutoff for alignment.'),
        ('-f', '--filter', 'Blast_Filter_Protocol_ID', 'int',
         ('Blast_Filter_Protocol_ID from MegaChipDB.'
          'Blast_Filter_Protocol table. If None no filter will be '
          'used.')),
        ('-l', '--likeTHD', 'templateTHD_ID','int',
         'use exsiting THD as template for unspecified values'),
        )

    ids =[ x[2] for x in optionsAvailable if x[2] != 'templateTHD_ID' ]

    for opt in optionsAvailable:
        if opt[3] is None:
            otype = 'string'
        else:
            otype = opt[3]
        parser.add_option(opt[0], opt[1], dest=opt[2], type=otype,
                          help=opt[4])
    (options, args) = parser.parse_args(args)


    defined = lambda x: getattr(options, x) is not None

    if options.templateTHD_ID != None:
        templateTHD=THD(options.templateTHD_ID)
        for opt in ids: 
            if not defined(opt):
                setattr(options,opt,getattr(templateTHD,opt))
        
    if not all(map(defined, ids)):
        print >> sys.stderr, ('Usage error: all arguments are required, '
                              'or derived from template THD')
        print >> sys.stderr, parser.format_help()
        sys.exit(1)

    thdc = THDCreator(options.Sequence_File_ID,
                      options.Alignment_Protocol_ID,
                      options.Alignment_Score_Limit,
                      options.Taxonomy_ID,
                      options.Blast_Filter_Protocol_ID)
    thdc.run()

if __name__ == '__main__':
    main(sys.argv[1:])
