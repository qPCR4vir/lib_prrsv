#!/usr/local/bin/python 
#
# makeNewTaxDB.py
#
# Down load current taxonomy dumps from NCBI and make a new
# MySQL database to hold the results.
#

__verison__ = tuple([int(x) for x in
                     '$Revision: 1.1 $'.split()[1].split('.')])
__author__ = "Kael Fischer"
