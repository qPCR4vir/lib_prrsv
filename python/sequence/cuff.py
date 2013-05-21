#!/usr/local/bin/python 
#
# sequence/cuff.py
#
# module for dealing with all the cufflinks,cuffmerge,
# cuffdiff input and output 
# 

from types import StringTypes
from utils import MultiDictCluster,safeOFW
from math import log

__verison__ = tuple([int(x) for x in
                     '$Revision: 1.1 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

class Gene_Tracking(object):
    """The cuffdiff genes.fpkm_tracking output.
    """

    def __init__ (self,fName):
        """
        """
        self.fName=fName
        self.headerLine=file(fName).next().strip()
        self.fields=self.headerLine.split()
        self.sampleNames=[x[:-5] for x in self.fields if x.endswith("_FPKM")]
        self.sampleIndexes = [i for i,x in enumerate(self.fields) if x.endswith("_FPKM")]
    
    def dataLines(self):
        """
        """
        

    def sampleIdNameData(self,log2transform):
        """
        """
        f=file(self.fName)
        f.next()
        for l in f:
            x=l.strip().split()
            row=[x[6],x[4]]
            for i in self.sampleIndexes:
                if x[i+3] == 'OK':
                    if log2transform:
                        try:
                            row.append(str(log(float(x[i]),2)))
                        except ValueError:
                            row.append('')
                    else:
                        row.append(x[i])
                else:
                    row.append('')
            yield(row)


    def cluster(self,outFile,log2transform=True):
        """
        """
        if type(outFile) in StringTypes:
            oFile=safeOFW(outFileName)
        else:
            oFile=outFile
        print >> oFile, '\t'.join(['LOCUS','NAME']+self.sampleNames)
        for row in self.sampleIdNameData(log2transform):
            if row[2:].count('') < len(row[2:]):
                print >> oFile, '\t'.join(row) 
        oFile.close()
    
    
