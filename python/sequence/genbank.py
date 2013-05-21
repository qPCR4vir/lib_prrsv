#!/usr/local/bin/python
#
#
#
# $Id: genbank.py,v 1.20 2012/10/24 23:07:56 kael Exp $ 
#
from sequence import complement

__version__ =  '$Revision: 1.20 $'

import sys
import re
import time
import types
import urllib

from __init__ import *

from datetime import datetime
from StringIO import StringIO
from traceback import print_exc

class GenBank_Lookup_Error (Exception):
    pass

class GenBank_Insertion_Error (Exception):
    pass


sectionStub = r'^(%s.+?)^\S'
subsectionStub = r'(^ {5}%s.+?)^ {5}\S'

sectionRE = re.compile(r'^(\S+) +',re.MULTILINE)
subsectionRE = re.compile('^ {5}(\S+) +',re.MULTILINE)

giRE =  re.compile(r'^VERSION.* GI:(\d+)',re.MULTILINE)

def txt2GI (recText):
    return int(giRE.search(recText).group(1))
               

def getNextRecordFromOpenFile(fHandle):
    """Look in file for the next record
    return text of the record
    """
    cSize =1024
    recFound = False
    recChunks = []
    try:
        fHandle.seek(-1,1)
    except IOError:
        pass
    sPos = fHandle.tell()
    
    gbr=None
    while True:
        cPos=fHandle.tell()
        c=fHandle.read(cSize)
        if c=='':
            return None
        if not recFound:
            
            locusPos=c.find('\nLOCUS')
            if sPos==0 and c.startswith('LOCUS'):
                locusPos=0
            elif locusPos == -1:
                continue
            if locusPos>0:
                locusPos+=1
            c=c[locusPos:]
            recFound=True
        else:
            locusPos=0

        if (len(recChunks)>0 and
            ((c.startswith('//\n') and recChunks[-1].endswith('\n'))
             or (c.startswith('\n') and recChunks[-1].endswith('\n//'))
             or (c.startswith('/\n') and recChunks[-1].endswith('\n/'))
             )):
            eorPos=0
        else:
            eorPos=c.find('\n//\n',locusPos)
            
        if eorPos == -1:
            recChunks.append(c)
        else:
            recChunks.append(c[:(eorPos+4)])
            gbr=''.join(recChunks)
            fHandle.seek(cPos-locusPos+eorPos)
            #print txt2GI(gbr),fHandle.tell(),cPos,locusPos,eorPos
            return gbr

        
    
        
class Record(StringSequence):
    """Simple GenBank Representation.
    """

    def __init__(self,recordContent=None,parse=True):
        """Make a new genbank record.
        recordContent can be a string or a file like thing.
        if recordContent is file like the files position will be moved 
        forward to the beginning of the next record (if any) of to EOF.

        By default the sequence and record text are placed in the
        'sequence' and 'gbr' instnace attributes.  If parse is not True,
        the the sequence attribute is left as None.  All the base StringSequence
        stuff and self.fasta won't work if you set parse to False, but you may do
        it to save memory in special cases. 
        """
        self.sequence=''
        self.title=''
        if recordContent != None:
            if type(recordContent) in types.StringTypes:
                self.gbr=recordContent
            else:
                self.gbr=getNextRecordFromOpenFile(recordContent)

            # force AttributeError if no record
            if self.gbr == None:
                raise AttributeError
            if parse:
                self.sequence=self.parseSequence()
            else:
                self.sequence=None
                

    def __getitem__(self,key):
        """Dictionary interface to Record's sections.
        """
        for s,c in self.sections():
            if s == key:
                return c
        raise KeyError, key

    def giTitle(self):
        """Returns the definition without the newline at the end.
        """
        return ' '.join(self['DEFINITION'].strip().split())

    
    def GI(self):
        """return the gi (integer) of the record
        """
        return txt2GI(self.gbr)

    gi=GI
                    
    def features(self):
        """generator of Feature 
        """
        for s,c in self.sections():
            if s == 'FEATURES':
                fSplit = subsectionRE.split(c)
                fSplit.pop(0)
                while len(fSplit) > 1:
                    yield(Feature(self,fSplit.pop(0),fSplit.pop(0)))
                


    def sourceNCBITaxID (self,debugLevel=0):
        """get the taxon of the first listed 'source' organism"""

        try:
            for f in self.features():
                 if f.ftype == 'source':
                     for ref in f.qualifiers['db_xref']:
                         (db,id) = ref.split(':')
                         if db == 'taxon':
                             return int(id)
            return None

        except:
            if debugLevel != 0:
                print_exc(None,sys.stdout)
            
            return None

    def sections(self):
        """Generator of sections
        """
        secSplit =  sectionRE.split(self.gbr)[1:]
        while len(secSplit) >1:
            yield(secSplit.pop(0),secSplit.pop(0))

    def section(self,sectionStr):
        """Returns the string content of a particular section.
        """
        sRslt = re.search(r'^(%s.+?)^\S'%sectionStr,
                         self.gbr,re.MULTILINE | re.DOTALL)
        if sRslt == None:
            return None
        else:
            return sRslt.groups()[0]

    def subsection(self,sectionStr,subsectionStr):
        """Returns the string content of a particular subsection.
        """
        mySection = self.section(sectionStr)

    def parseSequence(self):
        o=self.section("ORIGIN")
        if o != None:
            oLines = o.split('\n')[1:]
            return ''.join([''.join(l.split()[1:]) for l in oLines]).upper()
        else:
            return None
    
    def fasta(self):
        """ return a fasta Record
        """
        if self.sequence != None:
            import fasta
            return fasta.Record(title=str(self.GI()),sequence=self.sequence)

    def division(self):
        """return the division string from the LOCUS line
        """
        return self['LOCUS'].split()[-2]

    def gad(self):
        """Returns a string like
        <gi><TAB><accession><TAB><definition>
        for the specified `gi`
        """
        definition = re.sub(r'\s{1,}', ' ', self.giTitle())
        accessions = self['ACCESSION'].strip().split(' ')
        gads = []
        for acc in accessions:
            gads.append('%d\t%s\t%s' % (self.GI(), acc, definition))
        return '\n'.join(gads)

class Feature  (object):
    """GenBank Feature representation.
    """
    def __init__ (self,record,featureType,featureTxt):
        self.record=record
        self.qualifiers={}
        lineOne, remainder = featureTxt.split('\n',1)
        self.type=featureType
        self.ftype=self.type
        self.location = Location(record,lineOne.strip())
        refs =remainder.split('                     /')
        print refs
        for r in refs[1:]:
            print 0.0
            kv = r.strip().split('=',1)
            if len(kv) ==1:
                k=kv[0]
                v=None
                self.qualifiers[k]=v
            else:
                k,v = r.strip().split('=',1)
                k = k.strip('"')
                v = ' '.join(v.strip('"').split())
                if k=='db_xref' and ('db_xref' not in self.qualifiers):
                    print 0.1
                    self.qualifiers[k] = [v]
                elif k in self.qualifiers:
                    if type(k) != type([]):
                        print 0.2
                        try:
                            self.qualifiers[k].append(v)
                        except AttributeError:
                            self.qualifiers[k]=[self.qualifiers[k]].append(v)
                    else:
                        print 0.3
                        self.qualifiers[k].append(v)
                else:
                    self.qualifiers[k.strip('"')]=v.strip('"')
                    
    def sequence(self):
        """return the 
        """
        recSeq = self.record.sequence
        sections = []
        for region in self.location.regions:
            s = recSeq[region.start-1:region.end]
            if region.complement:
                s=reverseComplement(s)
            sections.append(s)
            
        return ''.join(sections)
    
    def regions(self):
        """return location.regions
        """
        return self.location.regions


class Location  (object):
    """Feature Location
    """
    
    complementRE = re.compile(r'^complement\((.*)\)$')
    joinRE = re.compile(r'^join\((.*)\)$')
    orderRE = re.compile(r'^order\((.*)\)$')
    singleBaseRE = re.compile(r'^(\d+)$')
    locRE = re.compile(r'^(<?)(\d+)([<.^>]+)(\d+)(>?)$')
    
    def __init__(self,record,locText):
        """Return a new location
        """
        self.record=record
        self.regions=[]
        self._parseRegions(locText)
        
    def _parseRegions(self,txt,complement=False):
        cMatch = self.complementRE.match(txt)
        jMatch = self.joinRE.match(txt)
        oMatch = self.orderRE.match(txt)
        sbMatch = self.singleBaseRE.match(txt)
        locMatch = self.locRE.match(txt)
              
        #print cMatch,jMatch,oMatch,sbMatch,locMatch,1000
        
        if cMatch != None:
            #print 'complement'
            self._parseRegions(cMatch.group(1),complement=True)
        elif jMatch != None:
            #print 'join'
            self._parseRegions(jMatch.group(1),complement=complement)
        elif oMatch != None:
            #print 'order'
            self._parseRegions(oMatch.group(1),complement=complement)
        else:
            # no parens
            locs = txt.split(',')
            if len(locs) > 1:
                for loc in locs:
                    self._parseRegions(loc,complement=complement)
    
            elif sbMatch != None:
                start = int(sbMatch.group(1))
                self.regions.append(Region(start,start,complement))        
            elif  locMatch != None:
                start=int(locMatch.group(2))
                end=int(locMatch.group(4))
                midA = locMatch.group(3)
                if midA == '..':
                    ambigStr =''
                else:
                    ambigStr = midA
                if locMatch.group(1) != None:
                    ambigStr += locMatch.group(1)
                if locMatch.group(5) != None:
                    ambigStr += locMatch.group(5)    
                    
                self.regions.append(Region(start,end,complement,ambiguity=ambigStr))
            
class Region (object):
    """A subregion of a location
    """        
    def __init__(self,start,end,complement,ambiguity=''):
        self.start=int(start)
        self.end=int(end)
        self.complement=complement
        self.ambiguity=ambiguity
        
    def range(self):
        return (self.start,self.end)

        
        
def GenBankIterator(fh):
    """A GenBank Record generator
    """
    closeOnDone=False
    if type(fh) in StringTypes:
        fh =open(fh)
        closeOnDone=True
    try:
        while True:
            yield Record(fh)
    except AttributeError:
        if closeOnDone:
            fh.close()
        raise StopIteration
iterator=GenBankIterator
        
        
def unparsedIterator(fh):
    """same as iterator (aka GenBankIterator) except
    the objects don't know their sequence.
    """
    closeOnDone=False
    if type(fh) in StringTypes:
        fh =open(fh)
        closeOnDone=True
    try:
        while True:
            yield Record(fh,parse=False)
    except AttributeError:
        if closeOnDone:
            fh.close()
        raise StopIteration
        
        
        
