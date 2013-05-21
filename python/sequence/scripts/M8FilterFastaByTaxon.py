#!/usr/local/bin/python 
#
# Fischer Lab new script template
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.2 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

import sys
import optparse
from sequence import blastNoSQL as B
from sequence import fasta as F
from utils import *

def parseListInt(option, opt_str, value, parser):
    setattr(parser.values, option.dest,
            [int(x) for x in value.split(',')])

def parseListStr(option, opt_str, value, parser):
    setattr(parser.values, option.dest,
            [str(x) for x in value.split(',')])

def main(sysargs=[]):

    oneLineUsage = (
        "Usage: %prog [options] -s <fasta/q file(s)> "
        "-a <BLAST result file(s)>" )

    # set a long desccription for the help message
    description=(
        """Filter sequence files based on alignment matches based on various
alignment characteristics.

Sequences in fasta (or fastq) formatted files that have hits in the
BLAST results are output to a specified output file (-o) or STDOUT.
BLAST results can be filtered by Tax_ID, score, e-value and length.
        """
        )
    
    op = optparse.OptionParser(
        oneLineUsage,description=description,
        version="%prog " + '.'.join([str(x) for x in __version__]))


    op.add_option('-t','--taxa',dest="taxa",type='string',
                  callback=parseListInt,action='callback',
                  help="Comma seperated list of taxIDs. Default: all taxa")
    op.add_option('-s','--seqfiles',dest="seqFiles",default='',
                  callback=parseListStr,action='callback',type='string',
                  help="Comma seperated paths to fasta/fastq files. [REQUIRED]" )
    op.add_option('-a','--alignments',dest="mbrFiles",default='',
                  callback=parseListStr,action='callback',type='string',
                  help=("Comma seperated paths BLAST alignemnts "
                        "(tabular format).  [REQUIRED]"))    
    op.add_option('-o','--outputfile',dest="outFileName",default='-',
                  help="Name of output file. '-' means: STDOUT.")
    op.add_option('-c','--clobber',action="store_true", default=False,
                  dest="clobber", help="allow output file to be overwritten.")
    op.add_option('-p','--append',action="store_true", default=False,
                  dest="append", help="Append output to output file.")    
    op.add_option('-l','--minlength',dest="minLen",default=0,type='int',
                  help="Minimum alignment length.")
    op.add_option('-b','--minscore',dest="minScore",default=0.0,type='float',
                  help="Minmum alignmet bitscore.")
    op.add_option('-e','--maxe',dest="maxE",default=100.0,type='float',
                  help="Maximum alignment e-value.")
    op.add_option('-S','--usesubject',dest="subjMatch",default=False,
                  action='store_true',
                  help="Match sequence titles to alignment subjects,"
                  "rather than the queries.")
    

    #
    # OPTION PATCHING FOR PRETTYER HELP
    #
    
    # all options should be defined above here

    # add defaults to help messages
    for o in op.option_list:
        if o.type == None:
            continue
        if o.help == None:
            o.help = "Default: %default"
        elif o.help.count('[REQUIRED]') or o.help.count('Default:'):
            continue
        else:
            o.help+= " Default: %default"

    # call the parser
    (opts,args) = op.parse_args(sysargs)

    try:
        #
        # OPTION CHECKING
        #
        print opts
        if opts.mbrFiles == None:
            raise RuntimeError, "No alignment files specified."
        elif opts.seqFiles == None:
            raise RuntimeError, "No sequence files specified."

        #
        # ARGUMENT CHECKING
        #

        if len(args) > 0:
            raise RuntimeError, "Unparsed arguments: %s." % (' '.join(args))
    
    except Exception, eData: 
        print >> sys.stderr, ("\nUsage Error: %s\n" %eData.message)
        print >> sys.stderr, op.format_help()
        return 1

    return (gameOn(args,opts))
    
    
def gameOn(args,opts):
    if opts.subjMatch:
        m8subjMatch = False
        m8qryMatch = True
    else:
        m8subjMatch = True
        m8qryMatch = False

    if opts.taxa == None or len(opts.taxa) == 0:
        ignoreT = True
    else:
        ignoreT = False

    oFile = safeOFW(opts.outFileName)
            

    fObj = B.M8TaxonFilter(taxonIDs=opts.taxa,minLen=opts.minLen,
                           minScore=opts.minScore, maxE=opts.maxE,
                           subjMatch=m8subjMatch, qryMatch=m8qryMatch,
                           ignoreTaxaAndGis=ignoreT)
    
    fastaG=fObj.extractFastaRecords(opts.mbrFiles,opts.seqFiles,outFile=oFile,
                                    append=opts.append,clobber=opts.clobber)

    return 0
    



    

####### LEAVE THIS ALONE ###########
# If run directly as a program
# this calls the main function.
# 
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


