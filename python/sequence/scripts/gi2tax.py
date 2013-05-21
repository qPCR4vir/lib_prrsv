#!/usr/local/bin/python
import sys
from utils import multiFile
import ncbi

def main(args):
    for gi in multiFile(args):
        gi=int(gi.strip())
        try:
            tid=str(ncbi.giTaxid(gi))
        except:
            tid = 'Lookup Error'

        try:
            name=str(ncbi.giTaxonScientificName(gi))
        except:
            name = 'Lookup Error'

        print '%s\t%s' % (tid,name)

    return 0

if __name__ == '__main__':
    sys.exit( main(sys.argv[1:]))
