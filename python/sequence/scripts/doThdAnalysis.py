#!/usr/local/bin/python

import sys
from optparse import OptionParser

from sequence.processor import THDAnalyzer


def main(args):
    parser = OptionParser('usage: %prog [options]')
    (options, args) = parser.parse_args(args)
    thda = THDAnalyzer(configFile=args[0], printResults=True)
    thda.run()

if __name__ == '__main__':
    main(sys.argv[1:])
