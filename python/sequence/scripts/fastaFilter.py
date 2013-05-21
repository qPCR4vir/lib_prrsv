#!/usr/local/bin/python
"""fastaFilter.py
By Kael Fischer, 2008

USAGE
   fastaFilter.py [-d] [-l <LZW size>] [-r <LZW ratio>] [-m <regexp>]
       [-p <regexp>] [-e<regexp name> -s|-t ] [-H <max homo stretch>]
       [-v] [-Q[q=#]] [-q n:s] [-B AAA,AAC,AAG,...] [-b <# barcodes>]
       [-P] [-a <[5 prime count]:[3 prime count]>]
       [-A <bad quality limit>] [--output=<output file name>] <file names ... > 

SUMMARY
     Perform several filtering steps on one or more fasta or fastq files.
     Output files are named by inserting the filtering options before
     the extension (of any) of the input file name.

     Output file is systematically named based on filtering options.
     
     If a regex (compared to sequence) is required or prohibited, a
     descriptive regexRequireName and/or regexProhibitName, can and should
     be given as well.

     Possible operations:
        -Filter on LZW parameters
        -remove duplicate records
        -require sequence regex match to the sequence, title, or both. 
        -require sequence regex no match to the sequence, title, or both.
        -remove records with more that a specified number of low quality
         positions.

OPTIONS
     -d        Remove duplicate sequences, '|#' is added to each title
               indicating # of times the sequence appears.
     -l #      Drop records with LZW compressed sizes lower than #.
     -r #.#    Drop records with (LZW Size / record length) lower than #.#.
     -m <re>   Regular expression to required to match.
     -p <re>   Regular expression that is prohibited.
     -e <str>  String to add to file name to indicate regex filtering .
     -s        Apply regex rules to fasta sequences.
     -t        Apply regex rules to fasta titles.
     -H <#>    Maximum length of longest homopolymeric stretch. 
     -v        verbose debugging on all input sequences.
     -Q        Output fastq files.
     -F <#>    Force quality to be this value.  Required if input is fasta
               and output is fastq.
     -q n:s    Drop records with n or more positions with quality value <= s.
     -5 <#>    Trim # of bases from the 5' end of the read before filtering.
     -3 <#>    Trim # of bases from the 3' end of the read before filtering.
     -B "X,Y"  Split records by indicated barcodes, separated by commas
               (e.g. "AAA,CCA,GGC,...").
     -b <#>    Split records by barcodes of length # (excluding parity bit).
     -P        Use parity bit.
     -a n5:n3  Trims the sequence from the 5' end when n5 adjacent
               qualities lower than the value specified with -A are
               found. Does the same thing with the 3' end using
               n3. Note that one of n5, n3 can be left empty to avoid
               trimming that end.
     -A #.#    Limit to consider a quality as good. Required if using -a.
     -L #      Drop records with length lower than #
     --output=<fn>   Output filename. If not specified will be generated
               based in the rest arguments.
     
     -h        print this message
     -?        print this message

BUGS
     Report suspect features to kael.fischer@gmail.com

"""

__version__ = "$Revision: 1.15 $".split(':')[-1].strip(' $') # don't edit this line.
                                                            # cvs will set it.
import sys
import re
import getopt
from sequence import fasta




def main(args=None):

    #default runtime options
    DEBUG = False

    MAXHOMO = 0
    FIVETRIM = 0
    THREETRIM = 0

    BARCODES = []
    NUMBARCODEBASES = 0
    PARITYBIT = False

    NODUPS = False
    LZWSIZE = 0
    LZWRATIO = 0.0

    REQRE = None
    PRORE = None

    RENAME = None

    ARESEQ = False
    ARETI = False

    FASTQOUT = False
    FORCEQUAL = None
    QUALFILTER = (None,None)

    QUALTRIM = (None, None)
    QUALTRIML = None

    LENGTHFILTER = 0

    OUTPUT = None

    DEFAULTS = (MAXHOMO,FIVETRIM,THREETRIM,
                BARCODES,NUMBARCODEBASES,PARITYBIT,
                NODUPS,LZWSIZE,LZWRATIO,
                REQRE,PRORE,RENAME,
                ARESEQ,ARETI,
                FASTQOUT,FORCEQUAL,QUALFILTER,QUALTRIM,
                QUALTRIML,LENGTHFILTER, OUTPUT)
    
    
    # parse command line options like '-h'
    # see pydoc getopt for option formats
    shortOptions="h?dstl:r:m:p:e:vH:QF:q:5:3:B:b:P:a:A:L:o"
    longOptions = ['output='] 
    try:
        opts, args = getopt.getopt(args, shortOptions, longOptions)
    except getopt.error, msg:
        # there is an unknown option!
            print msg      # prints the option error
            print __doc__  # prints the usage message from the top
            return (-2)

    try:
        # process options
        for option,optionArg in opts:
            if option=='-h' or option=='-?':
                print __doc__
                return(0)     # '0' = no error in UNIX
            elif option =='-d':
                NODUPS = True
            elif option=='-l':
                LZWSIZE = int(optionArg)
            elif option=='-r':
                LZWRATIO = float(optionArg)
            elif option=='-m':
                REQRE=re.compile(optionArg)
            elif option=='-p':
                PRORE=re.compile(optionArg)
            elif option =='-e':
                RENAME = optionArg
            elif option=='-s':
                ARESEQ = True
            elif option=='-t':
                ARETI = True
            elif option == '-v':
                DEBUG = True
            elif option == '-H':
                MAXHOMO = int(optionArg)
            elif option == '-Q':
                FASTQOUT = True
            elif option == '-F':
                FORCEQUAL=int(optionArg)
            elif option =='-q':
                QUALFILTER = tuple([int(x) for x in optionArg.split(':')])
            elif option == '-3':
                THREETRIM = int(optionArg)
            elif option == '-5':
                FIVETRIM = int(optionArg)
            elif option == '-B':
                BARCODES = optionArg.split(',')
            elif option == '-b':
                NUMBARCODEBASES = int(optionArg)
            elif option == '-P':
                PARITYBIT = True
            elif option == '-a':
                parts = map(lambda x: x.strip(), optionArg.split(':'))
                n5, n3 = None, None
                if len(parts[0]) > 0:
                    n5 = int(parts[0])
                if len(parts[1]) > 0:
                    n3 = int(parts[1])
                QUALTRIM = (n5, n3)
            elif option == '-A':
                QUALTRIML = float(optionArg)
            elif option == '-L':
                LENGTHFILTER = int(optionArg)
            elif option == '--output':
                OUTPUT = optionArg
            else:
                print "%s option not implemented" % option

    except StandardError, msg:
        print msg      # prints the option error
        print __doc__  # prints the usage message from the top
        return (-2)
       
    # check arguments
    # correct #, etc
    # the remaining command line arguments (after the
    # option processing) are in 'args'

    minArgs = 1 ## SET THAT
    maxArgs = None ## SET THAT
    
    argProblem = False
    if maxArgs != None and len(args) > maxArgs:
        print "Wrong number of arguments: %s found (expected max of %s)" % (lan(args),
                                                                            maxArgs)
        argProblem = True
    elif len(args) < minArgs:
        print "Wrong number of files: %s found (expected min of %s)" % (len(args),
                                                                            minArgs)
        argProblem = True

    elif ( ((REQRE or PRORE) and not (ARESEQ or ARETI or RENAME == None)) or
           ((ARESEQ or ARETI) and not (REQRE or PRORE or RENAME == None)) ):
        print ("If a regular expression is specified it has to be explicitly applied to\n"
               "either the sequence (-s) or title (-t) explicitly. -e must also be used to\n"
               "specify a name fragment for the output.\n")
        argProblem = True

    if FASTQOUT and FORCEQUAL == None and (

        True in [arg.endswith('.fasta') or arg.endswith('.fa') for arg in args]):

        print ("\nIf fastq output is desired (-Q) and input sequences are fasta format\n"
               "you must use -F to force a quality value.\n")
        argProblem=True


    if DEFAULTS == (MAXHOMO,FIVETRIM,THREETRIM,
                    BARCODES,NUMBARCODEBASES,PARITYBIT,
                    NODUPS,LZWSIZE,LZWRATIO,
                    REQRE,PRORE,RENAME,
                    ARESEQ,ARETI,
                    FASTQOUT,FORCEQUAL,QUALFILTER,QUALTRIM,
                    QUALTRIML,LENGTHFILTER, OUTPUT):
        print ("\nNO FILTERING OPTIONS SPECIFIED.\n"
               "NOTHING TO DO.\n")
        argProblem=True
        

    if argProblem:
        print __doc__
        return(-1)



    if NUMBARCODEBASES > 0:
        import sequence
        BARCODES = sequence.makeBarcodes(NUMBARCODEBASES, PARITYBIT)
            
    
    fasta.filterFasta(args,removeRepeats=NODUPS,minLZWsize=LZWSIZE,
                      minLZWratio=LZWRATIO,regexRequire=REQRE,
                      regexProhibit=PRORE,regexName=RENAME,
                      applyRegexSequence=ARESEQ,applyRegexTitle=ARETI,
                      maxHomo=MAXHOMO,fastqOutput=FASTQOUT,
                      forceQuality=FORCEQUAL,rejectQuality=QUALFILTER,
                      fivePTrim=FIVETRIM, threePTrim=THREETRIM,barcodes=BARCODES,
                      badQualityTrim=QUALTRIM, badQualityLimit=QUALTRIML,
                      minTrimLength=LENGTHFILTER, output=OUTPUT, DEBUG=DEBUG)
    
    return(0)  # we did it!


####### LEAVE THIS ALONE ###########
# If run directly as a program
# this block calls the main function.
# 
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
