#!/usr/local/bin/python 
#
# Fischer Lab new script template
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.1 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

import sys
import optparse
from types import *
import re


def lineDemangle(lines):
    """Generator of lists coresponding
    to line in a REBASE itype2.<ver> file, per:
    <enzyme name>,<prototype>,<recognition sequence with cleavage site>,
    <methylation site and type>,<commercial source>,<references>
    """
    for l in lines:
        if type(l) in StringTypes:
            l=l.strip('\n').split('\t')
        if len(l) == 6:
            yield l
            
def siteMatch(lines,pat):
    """return generator of lines that
    match pat re.
    """
    for l in lineDemangle(lines):
        if re.search(pat,l[2],re.IGNORECASE):
            yield(l)

def blunt(site):
    """True if site is blunt otherwose False
    """
    try:
        cPos=site.index('^')
    except ValueError:
        return False
    if len(site[:cPos]) == len(site[cPos+1:]):
        return True
    else:
        return False


####### LEAVE THIS ALONE ###########
# If run directly as a program
# this calls the main function.
# 
if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


