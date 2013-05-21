#!/usr/local/bin/python 
#
# plot.py
#
# Plot utilities
# 
#

__verison__ = tuple([int(x) for x in
                     '$Revision: 0.0 $'.split()[1].split('.')])
__author__ = "Kael Fischer"


import pylab as P
import numpy as N

def plotCorMat(cmat,tickLabels=None):
    f=P.figure()
    P.imshow(cmat,interpolation='nearest')
    if tickLabels != None:
        ylim = tuple(reversed(P.xlim()))
        P.yticks(N.arange(len(tickLabels)), tickLabels,
                 family='monospace')
        P.ylim(ylim)
        P.xticks(N.arange(len(tickLabels)), tickLabels,
                 family='monospace',rotation=90)
    return f
