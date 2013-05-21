#!/usr/local/bin/python 

import sys
from sequence import blastNoSQL as B

def __main__(args):
    label = args[0][:6]
    plts = B.M8Plots(args[0],30,label)
    plts.savePlots(5)



if __name__ == "__main__":
    sys.exit(__main__(sys.argv[1:]))
