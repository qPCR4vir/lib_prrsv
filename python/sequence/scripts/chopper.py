#!/usr/local/bin/python
__version__ = "$Revision: 1.4 $".split(':')[-1].strip(' $') # don't edit this line.
__author__ = "Kael Fischer"

import utils
from sequence import fasta
import optparse
import sys


def main(sysargs):

    op = optparse.OptionParser("Usage: %prog [options] <fasta files>",
                         version=__version__)
    op.add_option('-t','--tilesizes',dest="tileSizesL", default="60",
                  help="comma seperated list of tile sizes to generate. e.g. '50,70'.")
    op.add_option('-s','--shift',default=25,type='int',dest="shift",
                  help="Start position offset" )
    op.add_option('-c','--clobber',action="store_true", default=False,
                  dest="clobber", help="allow output file to be overwritten.")
    op.add_option('-E','--endtile',action="store_false", default=True,
                  dest="endtile", help="Suppress generation of an tile that covers the"
                  " last base in each recrod regardless of the -s and -t settings.")
    op.add_option('-o','--outputfile',dest="outFileName",default='-',
                  help="name of output file. '-' means: STDOUT")
    op.add_option('-P','--permissiveTitles',action="store_true", default=False,
                  dest="nonGiTitles",
                  help="allow non-gi fasta titles when calculating tile titles")


    for o in op.option_list:
        if o.help == None:
            o.help = "Default: %default"
        else:
            o.help+= " Default: %default"

    
    (opts,args) = op.parse_args(sysargs)

    if opts.nonGiTitles:
        tileTitleCallback = fasta.generalTileTitle
    else:
        tileTitleCallback = fasta.giTileTitle

    try:
        tileSizes = [int(x) for x in opts.tileSizesL.split(',') if len(x)>0]
        outFile = utils.safeOFW(opts.outFileName,clobber=opts.clobber)
    except Exception, eData:
        print >> sys.stderr, ("\nOption Error: %s\n" %eData.message)
        print >> sys.stderr, op.format_help()
        return 1

    try:
        if len(args) == 0:
            raise RuntimeError, "No input files specified."
    except Exception, eData: 
        print >> sys.stderr, ("\nUsage Error: %s\n" %eData.message)
        print >> sys.stderr, op.format_help()
        return 1
    

    for fn in args:
        fi = fasta.FastaIterator(fn)

        for rec in fi:
            for subRec in rec.DRNAsubrecords():
            #for subRec in [rec]:
                for size in tileSizes:
                    for tile in subRec.tile(size,startOffset=opts.shift,
                                            returnEndTile=opts.endtile,
                                            colwidth=size,
                                            titleCallback=tileTitleCallback):
                        tile.title+='.%s' % size

                        print >> outFile, tile
                        print >> outFile



        return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
