#!/usr/local/bin/python
#
# Runs the alignment command on the input filename and creates the THD
# from the result.
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.2 $'.split()[1].split('.')])
__author__ = "Julio Menendez"

import sys, os, re
from optparse import OptionParser

from viroinfo.megaseq import runCommand, AlignmentProtocol

def main(args):
    parser = OptionParser('usage: %prog [options]')
    optionsAvailable = (
        ('-i', '--input', 'inputfile', 'string'),
        ('-o', '--output', 'outputfile', 'string'),
        ('-a', '--alignment', 'Alignment_Protocol_ID', 'int'),
        )
    for oa in optionsAvailable:
        parser.add_option(oa[0], oa[1], dest=oa[2], type=oa[3])
    
    (options, args) = parser.parse_args(args)
    ap = AlignmentProtocol(options.Alignment_Protocol_ID)
    cArgs = ap.Arguments % {'database': ap.Database_Path}
    cArgs = cArgs.split(' -')
    cArgs = [cArgs[0],] + map(lambda x: '-' + x, cArgs[1:])
    cArgs.extend(('-i%s' % options.inputfile,
                  '-o%s' % options.outputfile))
    rc = runCommand((ap.Program, ) + tuple(cArgs), True)[1]
    if rc != 0:
        raise Exception

    failed = True
    if os.path.exists(options.outputfile):
        content = open(options.inputfile).read()
        nEntries = len(re.findall(r'^>', content, re.M))
        lastLine = runCommand(('tail', '-n1', options.outputfile))
        match = re.match(r'^#(?:.*)finished(?:.*) (\d+) (?:.*)',
                         lastLine)
        failed = match is None or int(match.group(1)) != nEntries

    if failed:
        return 1
    return 0

if __name__ == '__main__':
    try:
        retCode = main(sys.argv[1:])
    except:
        retCode = 1
    sys.exit(retCode)
