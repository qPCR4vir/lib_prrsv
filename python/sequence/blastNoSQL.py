"""BLAST Utilities that are free of DB connections.
"""
__version__ = tuple([int(ver) for ver in
                     "$Revision: 1.45 $".split()[1].split('.')])

__author__ = "Kael Fischer"


import copy
import os, os.path
import operator
import commands
import random
from StringIO import StringIO
import sys
import time
import xml.dom.minidom
import utils

from types import IntType, LongType
IntegerTypes=(IntType, LongType)

import numpy as N

from __init__ import *
from utils import *

NTBLASTDB_EXTS = ["nhr","nin","nsq" ]
PROBLASTDB_EXTS = ["phr","pin","psq"]

def ntDbExistsTime(fastaPath):
    """returns the ealiest timestamp (mtime) for a
    nt blast database made from fastaPath. Returns
    None if a required db file is missing.
    """
    fs=['.'.join((fastaPath,x)) for x in NTBLASTDB_EXTS]
    if False in [os.path.exists(f) for f in fs]:
        return None
    else:
        return min([os.path.getmtime(f) for f in fs])

def formatNtDb(fastaPath):
    """format nucleotide blast database.
    """
    os.system('formatdb -pf -i %s' %fastaPath)



def m8formatTuple(m8t):
    """Format an m8tuple (or HSP-like object with m8tuple method) as a
    string similar to a blast outputed m8 HSP  
    """
    if type(m8t) not in (TupleType, ListType):
        try:
            m8t = m8t.m8tuple()
        except AttributeError :
            raise ValueError, "invalid data for m8 formatting: %s" % m8t

    # there may be an additional key to strip from db layer
    if len(m8t) == 13:
        m8t = list(m8t[1:])
    elif len(m8t) == 12:
        m8t = list(m8t)
    else:
        raise ValueError, "invalid data for m8 formatting: %s" % m8t

    # input may contain longs - change to ints
    for i,v in enumerate(m8t):
        if type(v) == LongType:
            m8t[i] = int(v)

    return '\t'.join([str(x) for x in m8t])
    

def giFromM8name(m8name):
    try:
        return int(m8name)
    except:
        return int(m8name.split('|')[1])


def uniqueAlignments(m8inPaths,m8filtPaths,outputFh,allowLen=15):
    """output to outputFh alignments that are in inPaths,
    whose queries do not occour in m8filtPaths.

    inPaths and outPaths must be lists, or not string sequences
    """

    s,o = commands.getstatusoutput(
        "cat %s | "
        "awk '{if (NF==12 && $4>allowLen) {print $1;}}' | "
        "sort -u"
        % ' '.join(m8filtPaths))

    filtQtitles = noneDict(o.split('\n'))

    for f in m8inPaths:
        for l in file(f):
            qt = l.split()[0]
            if qt not in filtQtitles:
                print >> outputFh, l,
    
def l2M8 (line):
    """parse a line to a low overhead m8 tuple, where
    the elements have the correct data
    """
    (qryName,subjName,ident,length,gaps,mismatch,
     q_start,q_end,t_start,t_end,e,score)   =  line.strip().split()

    ident,e,score = map(float,(ident,e,score))
    length,gaps,mismatch,q_start,q_end,t_start,t_end = map(
        int,(length,gaps,mismatch,q_start,q_end,t_start,t_end))

    return (qryName,subjName,ident,length,gaps,mismatch,
            q_start,q_end,t_start,t_end,e,score)



def m8demangleTiles(m8):
    """convert query and subject tiles (if any) to
    gis and adjust coordinates.
    """
    (qryName,subjName,ident,length,gaps,mismatch,
     q_start,q_end,t_start,t_end,e,score) = m8

    if isTileTitle(qryName):
        gi,s,l = splitTileTitle(qryName)
        qryName=int(gi)
        q_start+=(s-1)
        q_end+=(s-1)
    
    if isTileTitle(subjName):
        gi,s,l = splitTileTitle(subjName)
        t_start+=(s-1)
        t_end+=(s-1)

    return (qryName,subjName,ident,length,gaps,mismatch,
            q_start,q_end,t_start,t_end,e,score)
    

def estimateM8rows(m8files):
    """Read a sample of the imput files and return an estimate
    of the total mumber of rows (hsps)
    """
    if len(m8files) < 6:
        m8samp = m8files
    else:
        m8samp = random.sample(m8files,(len(m8files)/20)+1)

    ct=0
    for r in multiFile(m8samp):
        ct+=1
    return ct *len(m8files)/len(m8samp)

  
def m8polarity(m8):
    """return 1 if qry and subj are pointing the same way, otherwise -1.
    """
    if m8[6] < m8[7]:
        # qry forward
        if m8[8] < m8[9]:
            return 1
        else:
            return -1
    else:
        if m8[8] < m8[9]:
            return -1
        else:
            return 1 

def m8TopHitsPerQuery(m8file):
    """return a generator of the top hit (eg the first hit)
    for each query in m8file
    """
    lastQ=''
    for t in m8tupleGenerator(m8file):
        if t[0] != lastQ:
            lastQ=t[0]
            yield t
    



def m8complete(mbrFile):
    """
    """
    if (not os.path.exists(mbrFile) or os.path.isdir(mbrFile)
        or os.stat(mbrFile).st_size == 0):
        return False
    
    if re.match('^#Mega BLAST run finished',lastLine(mbrFile)) != None:
        return True
    else:
        return False


def d1generator(d1File):
    """Returns a generator of records on a D1  BLAST results file.
    Example of record format returned:
    {'gaps': 1,
    'gaps_info': [{'pct_ident': 100,
                'q_end': 14,
                'q_start': 1,
                's_end': 2359442,
                's_start': 2359429}],
    'length': 13,
    'q_end': 14,
    'q_start': 1,
    'query': 'c2019_48',
    's_end': 2359442,
    's_start': 2359429,
    'strand': '+',
    'subject': '26111730'}    
    """
    currentRecord = None
    for l in multiFile(d1File):
        l = l.strip()
        if l.startswith("#'"):
            if currentRecord is not None:
                raise Exception, ('record header when there is a '
                                  'record already defined\n%s' % l)
            parts = l.split()
            names = parts[0].split('==')
            names = map(lambda x: x.strip("#'>"), names)
            subject_id = names[0]
            query_id = names[1]
            strand = query_id[0]
            query_id = query_id[1:]
            positions = map(lambda x: int(x.strip('()')), parts[1:-1])
            length_query = abs(positions[3]-positions[1])
            length_subject = abs(positions[2]-positions[0])
            currentRecord = dict(query=query_id, subject=subject_id, 
                                 length=max(length_subject,
                                            length_query),
                                 q_start=positions[1],
                                 q_end=positions[3],
                                 s_start=positions[0],
                                 s_end=positions[2],
                                 strand=strand,
                                 segments=[],
                                 _str_='')
        elif l.startswith('l'):
            if currentRecord is None:
                raise Exception, ('gap information outside a '
                                  'record\n%s' % l)
            parts = l.split()
            positions = map(lambda x: int(x), parts[1:-1])
            pctIdent = int(parts[-1].strip('()'))
            gap_info = dict(s_start=positions[0],
                            q_start=positions[1],
                            s_end=positions[2],
                            q_end=positions[3],
                            pct_ident=pctIdent)
            currentRecord['segments'].append(gap_info)
        elif l == '}':
            if currentRecord is None:
                raise Exception, 'record end without a record header'
            currentRecord['gaps'] = len(currentRecord['segments']) - 1
            yield currentRecord
            currentRecord = None                                 


class HugeBlast (PickleableObject):
    """Q-query index
    S-Subject index
    M-Sparse matrix of results (bit socres or p values)
    """

    def __init__ (self, m8files=None):
        """
        """
        self.M = None
        self.pairs = None
        self.Q={}
        self.S={}
        if type(m8files) == GeneratorType:
            self.m8gen=copy.deepcopy(m8files)
        else:
            self.m8gen=None
            
        if type(m8files) in StringTypes:
            m8files = (m8files,)
            
        self.m8files = m8files
#        if self.m8paths != None:
#            self.getQandS()


    def getQSpairs(self):
        """build set of Q/S pairs
        """
        self.pairs = set()
        for p in self.m8files:
            hc=0
            ts=time.time()
            lOld = len(self.pairs)
            for t in m8tupleGenerator(p):
                self.pairs.add((t[0],t[1]))
                hc+=1
            print p
            print "%d new hsps" %  hc
            print "%d new pairs" % (len(self.pairs)-lOld)
            print "%d seconds" % (time.time() - ts)
        return len(self.pairs)

    

    def getQandS(self):
        """build Q snd S indicies
        """
        for p in self.m8files:
            for t in m8tupleGenerator(p):
                self.Q[t[0]]=None
                self.S[t[1]]=None
            print len(self.Q),len(self.S)
        self.Q=dict(zip(self.Q.keys(),range(len(self.Q))))
        self.S=dict(zip(self.S.keys(),range(len(self.S))))
        return len(self.Q),len(self.S)


    def loadScores(self):
        """
        """
        import scipy.sparse
        if len(self.Q) + len(self.S) == 0:
            self.getQandS()
        self.M = scipy.sparse.dok_matrix((len(self.Q),len(self.S)),'f')
        
        for p in self.m8files:
            print p
            hi=0
            ts=time.time()
            for t in m8tupleGenerator(p):
                q,s = t[:2]
                bs = t[-1]
                self.M[self.Q[q],self.S[s]]+=bs
                hi+=1
            print "%d hsps" %  hi
            print "%d seconds" % (time.time() - ts)


class M8Set (object):
    """A set of -m 8 formated BLAST results
    """

    def __init__(self):
        self.alignments={}
        pass

    def parse(self,m8file):
        """
        """
        for l in m8file:
            (qryName,subjName,ident,length,gaps,mismatch,
             q_start,q_end,t_start,t_end,e,score)   =  l.split()

            ident,e,score = map(float,(ident,e,score))
            length,gaps,mismatch,q_start,q_end,t_start,t_end = map(
                int,(length,gaps,mismatch,q_start,q_end,t_start,t_end))

            if (qryName,subjName) not in self.alignments:
                self.alignments[(qryName,subjName)] = []
            self.alignments[(qryName,subjName)].append(
                [ident,length,gaps,mismatch,
                 q_start,q_end,t_start,t_end,e,score])

    def insertSearchHits(self,search):
        sequenceMap = search.sequenceMap()
        searchID=search.ID()

        rows = map(lambda x: (searchID,sequenceMap[x[0]],sequenceMap[x[1]]) ,
                   self.alignments.keys())

        
        blastDB.Hit.insertMany(
            ('Search_ID','Query_Sequence_ID','Subject_Sequence_ID'),
            rows,
            disableFKcheck=True,
            ignore=True)

        
    def insertSearchHsps(self,search):
        
        sequenceMap = search.sequenceMap()
        hitMap = {}
        

        def prepRow (alignmentItem):
            keySequences,myAlignments = alignmentItem
            (qryName,subjName) = keySequences
            try:
                hitID = hitMap[(qryName,subjName)]
            except KeyError:
                hitID = Hit._table(
                    Query_Sequence_ID=sequenceMap[qryName],
                    Subject_Sequence_ID=sequenceMap[subjName],
                    Search_ID=search.ID())[0].Hit_ID
                hitMap[(qryName,subjName)]=hitID
                
            return map(lambda x:  [hitID] + x, myAlignments)

        rows = []
        for a in self.alignments.items():
            rows.extend(prepRow(a))


        blastDB.Hsp.insertMany(('Hit_ID','Ident','Length',
                                'Gaps','Mismatch',
                                'Q_Start','Q_End',
                                'T_Start','T_End',
                                'E','Score'),
                               rows,
                               disableFKcheck=True,
                               ignore=True)


def m8LineGenerator(m8files):
     """generator filter for m8files that skips comments and
     non-HSP lines.
     """
     for l in multiFile(m8files):
         if (l.startswith('#') or l.startswith('Mega BLAST run') or
             l.startswith('Warning: no access to tty (Bad file descriptor).') or
             l.startswith('Thus no job control in this shell.')):
             continue
         else:
             yield l

  
def m8dict(m8file):
    """parse an m8 file and return a list of dictionaries of results
    """
    return list(m8generator(m8file))


def m8tupleGenerator(m8file):
    """Iterator over -m 8 BLAST results return tuples of this form: 

    (qryName,subjName,ident,length,gaps,mismatch,
    q_start,q_end,t_start,t_end,e,score)

    Skips comments and common SGE lines of output.
    """
    
    bareQryGi=True
    bareSubjGi=True

    for l in multiFile(m8file):
        if (l.startswith('#') or l.startswith('Mega BLAST run') or
            l.startswith('Warning: no access to tty (Bad file descriptor).') or
            l.startswith('Thus no job control in this shell.')):
            continue
        try:
            (qryName,subjName,ident,length,gaps,mismatch,
             q_start,q_end,t_start,t_end,e,score)   =  l.split()
        except ValueError:
            raise
            continue

        ident,e,score = map(float,(ident,e,score))
        length,gaps,mismatch,q_start,q_end,t_start,t_end = map(
            int,(length,gaps,mismatch,q_start,q_end,t_start,t_end))

        if bareQryGi:
            try:
                qryName = int(qryName)
            except:
                bareQryGi = False
        
        if bareSubjGi:
            try:
                subjName = int(subjName)
            except:
                bareSubjGi = False
              
        yield (qryName,subjName,ident,length,gaps,mismatch,
               q_start,q_end,t_start,t_end,e,score)



def m8generator(m8file):
    """An iterator over a file of -m 8 formatted BLAST results
    """
    for (qryName,subjName,ident,length,gaps,mismatch,
         q_start,q_end,t_start,t_end,e,score)  in m8tupleGenerator(m8file):

        yield {
            'query':qryName,
            'subject':subjName,
            'pctIdent':float(ident),
            'length':int(length),
            'gaps':int(gaps),
            'mismatch':int(mismatch),
            'q_start':int(q_start),
            'q_end':int(q_end),
            's_start':int(t_start),
            's_end':int(t_end),
            'e':float(e),
            'score':float(score),
            '_str_': '\t'.join([str(x) for x in (qryName,subjName,ident,length,gaps,mismatch,
                                q_start,q_end,t_start,t_end,e,score)])
            }    

class disposableMB (object):
    """Run MEGABLAST but do not store results in a database"""

    def __init__ (self,query,db,
                  options='-W24 -E10 -FF',deferJob=False, runOnGrid=True):

        self._tmpFiles=[]

        self.options=options
        if type(query) in (StringType,UnicodeType):
            # is the query a valid path?
            if os.path.exists(query) and os.access(query,os.F_OK):
                self.qPath=query
                self.qFile=file(query)
            else:
                self.qFile,self.qPath=mystemp(suffix=".fasta",dir='')
                self.qFile.write(query)
                #os.path
                self.qFile.write("\n")
                self.qFile.flush()
                self._tmpFiles.append(self.qPath)
                self.qFile.close()
                
        elif type(query) == FileType:
            self.qFile = query
            self.qPath = query.name

        if not (os.path.exists(db) and os.access(db,os.F_OK)):
            raise ArgumentError, "db: %s Not Found" %db
        else:
            self.dPath=db
            
        self.oFile,self.oPath=mystemp(suffix=".dmblast",dir='')
        self._tmpFiles.append(self.oPath)
        self.grid=runOnGrid
        

    def run(self):
        if self.grid:
            sgeCmd = "qrsh -cwd -N megablast "
        else:
            sgeCmd = ''
        self.blastCmd = (
            "%s/usr/local/bin/megablast  -i %s -d %s -o %s -D 3 %s -f -R" 
            % (sgeCmd,self.qPath,self.dPath,
               self.oPath,self.options))
        self.blastStat,self.blastOut = commands.getstatusoutput(self.blastCmd)
                
        if self.blastStat != 0:
            raise RuntimeError, (
                "---MEGABLAST FAILED---\nMEGABLAST command"
                ":%s\nEnd of output: %s\nExit Status: %s"
                % (self.blastCmd,self.blastOut[-25:],self.blastStat))
        
        return m8generator(file(self.oPath))


    def __del__ (self):
        for fn in self._tmpFiles:
            os.unlink(fn)

class DisposableNCBIBLAST (object):
    """Run blast /megablast
    """
    def __init__ (self,query,db,options='-p blastn',deferJob=False):

        self._tmpFiles=[]

        self.options=options
        if type(query) in (StringType,UnicodeType):
            if os.path.exists(query) and os.access(query,os.F_OK):
                self.qPath=query
                self.qFile=file(query)
            else:
                self.qFile,self.qPath=mystemp(suffix=".fasta")
                self.qFile.write(query)
                #os.path
                self.qFile.write("\n")
                self.qFile.flush()
                self._tmpFiles.append(self.qPath)
                
        elif type(query) == FileType:
            self.qFile = query
            self.qPath = query.name

        if not (os.path.exists(db) and os.access(db,os.F_OK)):
            raise ArgumentError, "db: %s Not Found" %db
        else:
            self.dPath=db
            
        self.oFile,self.oPath=mystemp(suffix=".dblast")
        self._tmpFiles.append(self.oPath)
 
    def __del__ (self):
        for fn in self._tmpFiles:
            os.unlink(fn)
    


class disposableBLAST (DisposableNCBIBLAST):
    """Run BLAST but do not store results in a database"""

    def run(self):
        blastCmd = ("/usr/local/bin/blastall  -i %s -d %s -o %s -m 8 %s" 
                    % (self.qPath,self.dPath,
                       self.oPath,self.options))
        blastStat,blastOut = commands.getstatusoutput(blastCmd)
                
        if blastStat != 0:
            raise RuntimeError, ("---BLAST FAILED---\nBLAST command:%s\nEnd of output: %s\nExit Status: %s"
                                 % (blastCmd,blastOut[-500:],blastStat))
        
        return m8generator(file(self.oPath))



DisposableBlast = disposableBLAST

class DisposableMB(DisposableNCBIBLAST):
    """Run MEGABLAST but do not store results in a database"""
    

    def run(self):
        blastCmd = ("/usr/local/bin/megablast  -i %s -d %s -o %s -D 3 %s -f -R" 
                    % (self.qPath,self.dPath,
                       self.oPath,self.options))
        blastStat,blastOut = commands.getstatusoutput(blastCmd)
                
        if blastStat != 0:
            raise RuntimeError, (
                "---MEGABLAST FAILED---\nMEGABLAST command:"
                "%s\nEnd of output: %s\nExit Status: %s"
                % (blastCmd,blastOut[-25:],blastStat))
        
        return m8generator(file(self.oPath))

class TaxHitComparison (object):

    def __init__(self, *distributions):
        self.distributions = distributions


class TaxHitDistribution (PickleableObject):
    """m8 parsing to taxon aggregated counts, various normalizations
    and comparisons.

    The dicts in use to count taxon hits are this form {(NCBI_taxID, count), ..}
    This should probably be subclassed.
    """
    def __init__(self,m8paths=[],m8eMax=10,m8eMin=0,subjCountPaths=[],
                 taxonomyDB=None):
        
        import ncbi, ncbi.giInfo
        reload(ncbi.giInfo)
        self.totalHits=0
        self.giCount = {}
        self.speciesCount={}
        self.genusCount = {}
        self.familyCount = {}
        self.taxonCount = {}
        self.taxonomyDB = taxonomyDB

        self.noFamilyGis=[]

        # these are populated when a reference is
        # specified
        self.speciesExpect = None
        self.genusExpect = None
        self.familyExpect = None



        if type(m8paths) in StringTypes:
            m8paths=(m8paths,)

        for p in m8paths:
            bareSubjGi = True
            for r in m8generator(file(p)):
                if r['e'] > m8eMax or r['e'] < m8eMin:
                    continue
                if type(r['subject']) != IntType:
                    bareSubjGi = False
                    
                if bareSubjGi:
                    gi = r['subject']
                else:
                    gi=int(r['subject'].split('|')[1])

                if gi not in self.giCount:
                    self.giCount[gi] = 1
                else:
                    self.giCount[gi] +=1

                self.totalHits+=1

        for p in subjCountPaths:
            for l in file(p):
                l=l.strip()
                try:
                    gi,n = [int(x) for x in l.split()]
                except ValueError:
                    
                    continue
                if gi not in self.giCount:
                    self.giCount[gi] = n
                else:
                    self.giCount[gi] +=n
            print '2'

        # now using giInfo!
        taxMap=ncbi.giInfo.taxMap(self.giCount.keys(),
                                  taxonomyDb=self.taxonomyDB)
        for gi in self.giCount.keys():
            try:
                t,s,g,f  = taxMap[gi]
                if f== None:
                    #print gi, g
                    self.noFamilyGis.append(gi)
            except KeyError:
                continue
            
            if t not in self.taxonCount:
                self.taxonCount[t] = self.giCount[gi]
            else:
                self.taxonCount[t] += self.giCount[gi]
            
            if s not in self.speciesCount:
                self.speciesCount[s] = self.giCount[gi]
            else:
                self.speciesCount[s] +=self.giCount[gi]

            if g not in self.genusCount:
                self.genusCount[g] = self.giCount[gi]
            else:
                self.genusCount[g] +=self.giCount[gi]

            if f not in self.familyCount:
                self.familyCount[f] = self.giCount[gi]
            else:
                self.familyCount[f] += self.giCount[gi]
                

    def expungeTaxon(self,taxid):
        """remove hits generated by gis beloning to clade
        specified by taxid.  Returns number of hits affected.
        """
        import ncbi, ncbi.giInfo, ncbi.taxonomy
        reload(ncbi.giInfo)
        if type(taxid) in IntegerTypes:
            t = taxid
        else:
            t=int(ncbi.taxonomy.stringLookup(taxid))
        cids = ncbi.taxonomy.Taxon(t).childIDs()
        try:
            n = self.taxonCount[t]
        except KeyError:
            print t
            n=0
        else:
            self.totalHits -= n
            s,g,f = ncbi.giInfo.tToSGF(t)
            self.speciesCount[s] -= n
            self.genusCount[g] -= n
            self.familyCount[f] -= n
            del self.taxonCount[t]

        return sum([n] + [self.expungeTaxon(c) for c in cids])
    
    def setBlastDBSizeExpected(self,DBpath):
        """Set the expected values for each taxa to those based number of nt
        in DB.
        """
        import blast
        ffRec = blast.FastaFile(Filename=DBath)

        sRows = ncbi.ncbiDB.fetchall (
            """SELECT G.Species_Tax_ID, sum(G.Length) FROM Blast_Results.Fasta_File F
            LEFT JOIN Blast_Results.Sequence S ON S.Fasta_File_ID=F.Fasta_File_ID
            LEFT JOIN NCBI.Gi_Info G ON S.Gi=G.Gi
            WHERE F.Fasta_File_ID = %s
            Group By G.Species_Tax_ID
            """ %(ffRec.ID()))
        self.speciesExpect=dict(sRows)

        gRows = ncbi.ncbiDB.fetchall (
            """SELECT G.Genus_Tax_ID, sum(G.Length) FROM Blast_Results.Fasta_File F
            LEFT JOIN Blast_Results.Sequence S ON S.Fasta_File_ID=F.Fasta_File_ID
            LEFT JOIN NCBI.Gi_Info G ON S.Gi=G.Gi
            WHERE F.Fasta_File_ID = %s
            Group By G.Genus_Tax_ID
            """ %(ffRec.ID()))
        self.genusExpect=dict(gRows)

        fRows = ncbi.ncbiDB.fetchall (
            """SELECT G.Family_Tax_ID, sum(G.Length) FROM Blast_Results.Fasta_File F
            LEFT JOIN Blast_Results.Sequence S ON S.Fasta_File_ID=F.Fasta_File_ID
            LEFT JOIN NCBI.Gi_Info G ON S.Gi=G.Gi
            WHERE F.Fasta_File_ID = %s
            Group By G.Family_Tax_ID
            """ %(ffRec.ID()))
        self.familyExpect=dict(fRows)


    def setExpected(self,refDist):
        """Set the expected hit distribution to another TaxHitDistribution.
        """

        if not isinstance(refDist,TaxHitDistribution):
            raise ValueError , "refDist must be a  TaxHitDistribution instance"

        self.speciesExpect = refDist.speciesCount
        self.genusExpect = refDist.genusCount
        self.familyExpect = refDist.familyCount

    def hitCount(self):
        """Return number of hits in distribution.
        """
        return self.totalHits

    def normalize(self,n=100):
        """nomalize to n observations
        """

        (self.giCount,self.taxonCount,
         self.speciesCount,self.genusCount,
         self.familyCount) = [normDict(d,n) for d in (self.giCount,
                                                      self.taxonCount,
                                                      self.speciesCount,
                                                      self.genusCount,
                                                      self.familyCount)]


    def pForStudentsT_P(self,obsDict,expectedDict,key):
        from scipy import stats
        from scipy import array,float64
        try:
            oOfKey =  obsDict[key]
        except KeyError:
            oOfKey = 0

        try:
            eOfKey = expectedDict[key]
        except KeyError:
            eOfKey = 0
   

    @classmethod
    def chiSquared(cls,obsDict,expectedDict,key):
        from scipy import stats
        from scipy import array,float64

        try:
            oOfKey =  obsDict[key]
        except KeyError:
            oOfKey = 0

        try:
            eOfKey = expectedDict[key]
        except KeyError:
            eOfKey = 0
       
        obs = array((oOfKey, sum(obsDict.values())-oOfKey),dtype=float64)
        #print obs
        expt = array((eOfKey, sum(expectedDict.values())-eOfKey),dtype=float64)
        #print expt
        if sum(obs) > sum(expt):
            obs*=float(sum(expt))/float(sum(obs))
        else:
            expt*=float(sum(obs))/float(sum(expt))
        

        x2,p = stats.chisquare(obs,expt)
        #
        if obs[0] > expt[0]:
            overObs = True
            
        else:
            overObs = False

        #if overObs:
        #    print "%s over-observed" % key
            
        return x2,p,obs[0],expt[0],overObs



    def cladeReport (self,obsDict,expectedDict,pLimit=1,
                     overCalledOnly=False,taxIDs=None,
                     familyRestrict=None, genusRestrict=None):
        """calculate chance distributions of clade hits are different by chance.
        Return report.
        """

        from ncbi import taxonomy
        Taxon = taxonomy.Taxon
        if self.taxonomyDB is not None:
            newDB = taxonomy.taxDB.clone(self.taxonomyDB)
            Taxon._table = newDB.Nodes
        
        rows = []
        if taxIDs == None:
            taxIDs = expectedDict.keys()

        restrictTaxon=familyRestrict
        if genusRestrict != None:
            restrictTaxon = genusRestrict

        if restrictTaxon != None:
            rtID = restrictTaxon.ID()
            rTaxIDs = []
            for tid in taxIDs:
                try:
                    tax=Taxon(tid)
                    if tax.idInLineage(rtID):
                        rTaxIDs.append(tid)
                except:
                    pass
            taxIDs = rTaxIDs


        for s in taxIDs:
            x2,p,o,e,overObs = self.chiSquared(obsDict,expectedDict,s)
            if (N.isnan(p) or
                p > pLimit or
                (overCalledOnly and not overObs)):
                continue
            try:
                tax=Taxon(s)
            except:
                tax = s

            rows.append((tax,p,o,e))

        rows.sort(key=lambda x: (x[2]<x[3],x[1],-x[2]))


        return '\n'.join(['\t'.join(
            [str(tax), str(p),"%.2f"%o,"%.2f"%e]) for tax,p,o,e in rows] )


    def topNTaxa(self,numTaxa=10):
        """returns the top N taxa (default = 10) and their counts as
        a list of tuples
        """
        return(sorted(self.taxonCount.items(), lambda x, y: cmp(x[1], y[1]),
                      reverse=True)[0:numTaxa])

        
    def familyReport(self,expectDist=None,**args):
        """like cladeReport but for all families
        """
        if expectDist==None:
            e = self.familyExpect
        else:
            e = expectDist.familyCount
                         
        return self.cladeReport(self.familyCount,e,**args)

    def genusReport(self,expectDist=None,**args):
        """like cladeReport but for all genera
        """
        if expectDist==None:
            e = self.genusExpect
        else:
            e = expectDist.genusCount
                         
        return self.cladeReport(self.genusCount,e,**args)

    def speciesReport(self,expectDist=None,**args):
        """like cladeReport but for all species
        """
        if expectDist==None:
            e = self.speciesExpect
        else:
            e = expectDist.speciesCount
                         
        return self.cladeReport(self.speciesCount,e,**args)

    def familyPReport(self,expectDist=None,**args):
        """like cladeReport but for all families
        """
        if expectDist==None:
            e = self.familyExpect
        else:
            e = expectDist.familyCount
                         
        return self.cladeReport(expectDist,**args)

    def genusPReport(self,expectDist=None,**args):
        """like cladeReport but for all genera
        """
        if expectDist==None:
            e = self.genusExpect
        else:
            e = expectDist.genusCount
                         
        return self.cladeReport(self.genusCount,e,**args)

    def speciesPReport(self,expectDist=None,**args):
        """like cladeReport but for all species
        """
        if expectDist==None:
            e = self.speciesExpect
        else:
            e = expectDist.speciesCount
                         
        return self.cladeReport(self.speciesCount,e,**args)



    def pieChart(self,obsDict,showN=None,showPct=5,
                 outFileName=None,dpi=150,axes=None):
        """Pie chart is ruturned on the given axes of saved in file (
        outFileName), counting things in obsDict { (NCBI taxID,count), ...)
        """

        import matplotlib as M
        from matplotlib import rcParams as Mrc
        Mrc['useafm']=True
        
        import pylab as P
        
        from viroinfo import taxon

        # if taxID is none slice is ignored
        tCounts = [ (ct,tid) for tid,ct in  obsDict.items() if tid != None]
        tCounts.sort()
        tCounts.reverse()

        if showN == None:
            if showPct == None:
                showN=10

            else:
                n=float(sum([c for c,tid in tCounts]))
                ctCutoff = n/100.0*float(showPct)
                showN = sum([c>ctCutoff for  c,tid in tCounts])
        
        #debugging
        print tCounts

        labels=[]
        counts=[]
        
        for c,s in tCounts[:showN]:
            counts.append(c)
            try:
                # should be changed to use ncbi.giInfo
                tax=taxon.Taxon(NCBI_TaxID=s)
            except:
                # if that doesn't compute use number
                tax = s
                
            labels.append(str(tax)+ "\n" + str(c))

        # is this right?
        counts.append(sum(obsDict.values())-sum(counts))
        labels.append('other\n' + str(counts[-1]))

        if axes == None:
            # make figure and save it in a file
            fig=P.figure(figsize=(5,5)).add_axes([0.1, 0.1, 0.8, 0.8])
            axes=P.gca()
            axes.pie(counts,labels=labels)
            if outFileName==None:
                P.show()
            else:
                fig.savefig(outFileName,dpi=dpi)
        else:
            # work on user supplied axes
            axes.apply_aspect(1)
            axes.pie(counts,labels=labels)
            return axes

    def familyPie(self,**args):
        rv = self.pieChart(self.familyCount,**args)
        return rv

    def genusPie(self,**args):
        rv = self.pieChart(self.genusCount,**args)
        return rv

    def speciesPie(self,**args):
        rv = self.pieChart(self.speciesCount,**args)
        return rv

    def __add__(self,other):
        """Return a new TaxHitDistribution instance. The observatoin
        counts are the sum of this and the 'other' instance's
        """
        rObj = copy.deepcopy(self)

        rObj.totalHits = rObj.totalHits+other.totalHits

        rObj.giCount = addDict(rObj.giCount,other.giCount)
        rObj.taxonCount = addDict(rObj.taxonCount,other.taxonCount)
        rObj.speciesCount = addDict(rObj.speciesCount,other.speciesCount)
        rObj.genusCount = addDict(rObj.genusCount,other.genusCount)
        rObj.familyCount = addDict(rObj.familyCount,other.familyCount)

        return rObj

    def __sub__(self,other):

        """Return a new TaxHitDistribution instance. The observatoin
        counts are the difference of this and the 'other' instance's
        """
        rObj = copy.deepcopy(self)

        rObj.totalHits = rObj.totalHits+other.totalHits

        rObj.giCount = subtractDict(rObj.giCount,other.giCount)
        rObj.taxonCount = subtractDict(rObj.taxonCount,other.taxonCount)
        rObj.speciesCount = subtractDict(rObj.speciesCount,other.speciesCount)
        rObj.genusCount = subtractDict(rObj.genusCount,other.genusCount)
        rObj.familyCount = subtractDict(rObj.familyCount,other.familyCount)

        return rObj

#
# Functional report interface
#
def fRpt(obs,expt,**kwArgs):
    return obs.familyReport(expt,**kwArgs)


def gRpt(obs,expt,**kwArgs):
    return obs.genusReport(expt,**kwArgs)


def sRpt(obs,expt,**kwArgs):
    return obs.speciesReport(expt,**kwArgs)


class TaxHitDistributionSet (MultiDictCluster,PickleableObject):
    """A set of TaxHitDistributions.
    """

    def __init__(self,files=[], taxonomyDB=None):
        #super(MultiDictCluster,self).__init__(self)
        MultiDictCluster.__init__(self)
        self.taxonomyDB = taxonomyDB
        self.loadDistSet(files)



    def bcFromName(self,filePath):
        """make a key for the Set (account for barcode)
        """
        fn=os.path.split(filePath)[-1]
        mo =re.search(r'\.([A-Z]{3,})\.',fn)
        if  mo == None:
            return fn.split('.')[0]
        else:
            return fn[:4]+mo.group(1)
            

    def familyIDs(self,codes=[]):
        if len(codes) == 0:
            codes = self.keys()
        return union([v.familyCount.keys() for k,v in self.items()
                      if k in codes])


    def genusIDs(self,codes=[]):
        if len(codes) == 0:
            codes = self.keys()
        return union([v.genusCount.keys() for k,v in self.items()
                      if k in codes])


    def speciesIDs(self,codes=[]):
        if len(codes) == 0:
            codes = self.keys()
        return union([v.speciesCount.keys() for k,v in self.items()
                      if k in codes])


    def taxonIDs(self,codes=[]):
        if len(codes) == 0:
            codes = self.keys()
        return union([v.taxonCount.keys() for k,v in self.items()
                      if k in codes])


    def giIDs(self,codes=[]):
        if len(codes) == 0:
            codes = self.keys()
        return union([v.giCount.keys() for k,v in self.items() if k in codes])

    def loadDistSet(self,files):
        for c in files:
            self[self.bcFromName(c)]=TaxHitDistribution.restore("%s"%( c))

    def normalize(self,n):
        for d in self.values():
            d.normalize(n)

    def expungeTaxon(self,t):
        return [d.expungeTaxon(t) for d in self.values()]


    def copy(self):
        return copy.deepcopy(self)


    def nXnCompare(self,rptFcn=fRpt,pLimit=1e-100,eCodes=None,outF=sys.stdout):
        try:
            iter(rptFcn)
        except:
            rptFcn=(rptFcn,)

        if eCodes == None:
            eCodes = self.keys()
        elif type(eCodes) in StringTypes:
            eCodes = (eCodes,)

        ref = reduce(operator.add, [self[x] for x in self.keys()])

        for eCode in eCodes:
            for rf in rptFcn:
                for l in rf(self[eCode],ref,pLimit=pLimit,
                            overCalledOnly=True).split('\n'):
                    print >> outF, "%s vs %s\t%s" % (eCode,'ref',l)
                for l in rf(self[eCode],ref-self[eCode],pLimit=pLimit,
                            overCalledOnly=True).split('\n'):
                    print >> outF, "%s vs ref-%s\t%s" % (eCode,eCode,l)
                for cCode in self.keys():
                    if cCode == eCode:
                        continue
                    for l in rf(self[eCode],self[cCode],pLimit=pLimit,
                                overCalledOnly=True).split('\n'):
                        print >> outF, "%s vs %s\t%s" % (eCode,cCode,l)


    def fCompare(eList,cList,self,pLimit=1e-100,outF=sys.stdout):
        for e in eList:
            for c in cList:
                for l in self[e].familyReport(self[c],pLimit=pLimit,
                                              overCalledOnly=True).split('\n'):
                    print >> outF, "%s vs %s\t%s" % (e,c,l)


    def gCompare(eList,cList,self,pLimit=1e-100,outF=sys.stdout):
        for e in eList:
            for c in cList:
                for l in self[e].genusReport(self[c],pLimit=pLimit,
                                             overCalledOnly=True).split('\n'):
                    print >> outF, "%s vs %s\t%s" % (e,c,l)


    def sCompare(eList,cList,self,pLimit=1e-100,outF=sys.stdout):
        for e in eList:
            for c in cList:
                for l in self[e].speciesReport(self[c],pLimit=pLimit,
                                               overCalledOnly=True).split('\n'):
                    print >> outF, "%s vs %s\t%s" % (e,c,l)


    def makeCluster(self,**kwArgs):
        #from ncbi import taxonomy
        kwArgs['eCallback']=lambda x: x.taxonCount
        rv=TaxHitCluster(self.taxonomyDB)
        rv.addNestedDict(self,**kwArgs)
        #rv.kSqlNameMap(taxonomy.Taxon)
        return rv


    def fgstCluster(self,**kwArgs):
        from ncbi import taxonomy
        rv = self.tCluster(**kwArgs)
        [rv.update(x) for  x in (
            [self.sCluster(**kwArgs),
             self.gCluster(**kwArgs),
             self.fCluster(**kwArgs)])]
        #rv.kSqlNameMap(taxonomy.Taxon,
        #               objectCallback= lambda x: x.fgstStr())
        return rv

    def fgsCluster(self,**kwArgs):
        rv=self.sCluster(**kwArgs)
        [rv.update(x) for  x in (
            [self.gCluster(**kwArgs),
             self.fCluster(**kwArgs)])]
        return rv


    def tCluster(self,**kwArgs):
        return self.makeCluster(eCallback=lambda x: x.taxonCount,**kwArgs)


    def sCluster(self,**kwArgs):
        return self.makeCluster(eCallback=lambda x: x.speciesCount,**kwArgs)


    def gCluster(self,**kwArgs):
        return self.makeCluster(eCallback=lambda x: x.genusCount,**kwArgs)


    def fCluster(self,**kwArgs):
        return self.makeCluster(eCallback=lambda x: x.familyCount,**kwArgs)


    def combine(self,codes=[]):
        if len(codes) == 0:
            codes = self.keys()
        return nanSum([v for k,v in self.items() if k in codes])
    

    def pTTest(self,thdAttr,eCode,catKey,df=1,
               dropCodes=None,reference=None):
        """Helper function for cladeTReport.
        """
        from scipy import stats
        import numpy as N

        #print dropCodes
        
        def thdDict(thd):
            return getattr(thd,thdAttr)
        
        #print dropCodes
        
        try:
            obs = float(thdDict(self[eCode])[catKey])
        except KeyError:
            obs = 0.

        if dropCodes == None:
            dropCodes = []
        if eCode not in dropCodes:
            dropCodes.append(eCode)


        expt = [] 

        if reference == None:
            reference = self
        for k,v in reference.items():
            #print v
            if k in dropCodes or k == eCode:
                continue
            va=thdDict(v)
            if catKey in va:
                expt.append(va[catKey])
            else:
                expt.append(0.)
        expt=N.array(expt,dtype=N.float64)
        #print expt


        t=(N.mean(expt)-obs)/(N.std(expt)/N.sqrt(len(expt)))
        p=1-stats.t.cdf(abs(t),df)


        if obs > N.mean(expt):
            overObs = True
            
        else:
            overObs = False

        #if overObs:
        #    print "%s over-observed" % key
        dropCodes=[]
        return t,p,obs,N.mean(expt),overObs,reference



    def cladeTReport(self,thdAttr,eCodes,pLimit=1, eLimit=0,titleAnnotation='',
                     overCalledOnly=False,taxIDs=None,
                     familyRestrict=None, genusRestrict=None,
                     dropCodes=None, everyLineECode=False,
                     reference=None, enableRefWarn=False,
                     returnValues=False):
        """calculate chance distributions of clade hits are different by chance.
        Return report.
        """
        from ncbi import taxonomy
        Taxon = taxonomy.Taxon
        if self.taxonomyDB is not None:
            newDB = taxonomy.taxDB.clone(self.taxonomyDB)
            Taxon._table = newDB.Nodes

        if type(eCodes) in StringTypes:
            eCodes=[eCodes]
        elif type(dropCodes) in StringTypes:
            dropCodes=[dropCodes]

        eCodes.sort()

        outLines=[]

        if taxIDs == None:
            taxIDs = self.taxonIDs()

        restrictTaxon=familyRestrict
        if genusRestrict != None:
            restrictTaxon = genusRestrict

        if restrictTaxon != None:
            rtID = restrictTaxon.ID()
            rTaxIDs = []
            for tid in taxIDs:
                try:
                    tax=Taxon(tid)
                    if tax.idInLineage(rtID):
                        rTaxIDs.append(tid)
                except:
                    pass
            taxIDs = rTaxIDs

        for eCode in eCodes:
            outLines.append("Student's T %s report: %s" %
                            (titleAnnotation,eCode))
            if restrictTaxon != None:
                outLines.append("Clade restriction: %s" % (restrictTaxon))
            if overCalledOnly:
                outLines.append("Over represented taxa only")
            outLines.append("pLimit: %s\teLimit: %s"% (pLimit,eLimit))
            outLines.append("Refernce Keys: %s" % ','.join([]) )

                                
            rows=[]
            for s in taxIDs:
                x2,p,o,e,overObs,ref = self.pTTest(
                    thdAttr,eCode,s,dropCodes=dropCodes,
                    reference=reference)
                if (N.isnan(p) or
                    p > pLimit or e < eLimit or
                    (overCalledOnly and not overObs)):
                    continue
                try:
                    tax=Taxon(s)
                except:
                    tax = s
                #print ref.keys()
                if eCode in ref and enableRefWarn:
                    refWarn = '***%s*** in reference'% eCode
                else:
                    refWarn=''
                rows.append((tax,p,o,e,refWarn))

            if returnValues:
                return rows
            
            rows.sort(key=lambda x: (x[2]<x[3],x[1],-x[2]))
            outLines.extend([
                '\t'.join([str(eCode), str(tax), str(p),"%.4f"%o,"%.4f"%e,refWarn])
                for tax,p,o,e,refWarn in rows])
            outLines.append('')

        return '\n'.join(outLines)


    def familyReport(self,**args):
        return self.cladeTReport('familyCount',self.keys(),
                                 taxIDs=self.familyIDs(),**args)


    def genusReport(self,**args):
        return self.cladeTReport('genusCount',self.keys(),
                                 taxIDs=self.genusIDs(),**args)


    def speciesReport(self,**args):
        return self.cladeTReport('speciesCount',self.keys(),
                                 taxIDs=self.speciesIDs(),**args)
    

    def taxonCtCor(self):
        from scipy import stats
        codeCts = [ (k,v.taxonCount) for k,v in self.items() ]
        ary=dicts2array(*[x[1] for x in codeCts ])[:,1:]
        ary=ary.astype(N.float)
        cIdx = stats.corrcoef(ary,rowvar=False).argsort(axis=0)[:,0]
        ary=ary[:,cIdx]
        corMat=stats.corrcoef(ary)
        codes=list(N.array([x[0] for x in codeCts])[cIdx])
        
        return corMat, codes


    def plotTxCor(self,labels=None,scaleBar=True,outFile=None,**imshowargs):
        import pylab as P
        P.figure(figsize=(11,9))
        cm,ck = self.taxonCtCor()
        P.imshow(cm,interpolation='nearest',**imshowargs)

        if labels != None:
            ck = list(labels)
        P.yticks(P.arange(len(ck)),ck,size='x-small',family='monospace')
        P.xticks(P.arange(len(ck)),ck,rotation=90,size='x-small',family='monospace')
        P.ylim(list(reversed(P.xlim())))
        if scaleBar:
            P.colorbar()
        if type(outFile) in StringTypes:
            P.gcf.savefig(outFile)
        return P.gcf()
               

    def dissimilarRef(self,excludeCodes,corrMax=0.6,topPctMax=1.):
        """
        """

        rv=self.__class__()
        if corrMax == 1:
            for k in self.keys():
                if k not in excludeCodes:
                    rv[k]=self[k]
            return rv
        #rv=[]
        
        if type(excludeCodes) in StringTypes:
            excludeCodes=[excludeCodes]
        elif type(excludeCodes) == TupleType:
            excludeCodes=list(excludeCodes)
        if topPctMax < 1:
            for code in self.keys():
                if code in excludeCodes:
                    continue
                thd=self[code]
                if float(thd.topNTaxa(1)[0][1])/sum(thd.taxonCount.values()) > topPctMax:
                    print "dropping: " + code
                    excludeCodes.append(code)
                
        cm,ck = self.taxonCtCor()

        keys =ck

        for ec in excludeCodes:
            idx=keys.index(ec)
            idxs = range(len(keys))
            del idxs[idx]
            idxs=N.array(idxs)
            del keys[idx]
            cm=cm[:,idxs][idxs,:]


        while cm.size >0:
            mostDiffIdx = cm.mean(axis=0).argsort()[0]
            rv[keys[mostDiffIdx]] = self[keys[mostDiffIdx]]
            #rv.append(self[keys[mostDiffIdx]])

            idxs = range(len(keys))

            dropIdx=[]
            for i in range(len(keys)):
                if cm[mostDiffIdx,i] > corrMax:
                    dropIdx.append(i)

            for i in reversed(dropIdx):
                del keys[i]
                del idxs[i]

            if len(idxs)==0:
                break
            idxs=N.array(idxs)
            #print idxs
            cm=cm[:,idxs][idxs,:]

        return rv

        
    def topTaxaReport(self,eCodes,pLimit=1, eLimit=0,titleAnnotation='',
                      taxIDs=None, familyRestrict=None, genusRestrict=None,
                      dropCodes=None, everyLineECode=False,
                      topN=10):
                      
        """calculate chance distributions of clade hits are different by chance.
        Return report.
        """

        from ncbi import taxonomy
        Taxon = taxonomy.Taxon

        if type(eCodes) in StringTypes:
            eCodes=[eCodes]
        elif type(dropCodes) in StringTypes:
            dropCodes=[dropCodes]

        eCodes.sort()

        outLines=[]

        if taxIDs == None:
            taxIDs = self.taxonIDs()

        restrictTaxon=familyRestrict
        if genusRestrict != None:
            restrictTaxon = genusRestrict

        if restrictTaxon != None:
            rtID = restrictTaxon.ID()
            rTaxIDs = []
            for tid in taxIDs:
                try:
                    tax=Taxon(tid)
                    if tax.idInLineage(rtID):
                        rTaxIDs.append(tid)
                except:
                    pass
            taxIDs = rTaxIDs


        for eCode in eCodes:
            rows = []
            outLines.append("Top Taxa %s report: %s" %
                            (titleAnnotation,eCode))
            if restrictTaxon != None:
                outLines.append("Clade restriction: %s" % (restrictTaxon))
            
            thd = self[eCode]

            for tid,n in thd.topNTaxa(topN):
                try:
                    tax=Taxon(tid)
                except:
                    tax = tid

                rows.append((tax,n))

            outLines.extend([
                '\t'.join([eCode, str(tax), str(n)])
                for tax,n in rows])
            outLines.append('')

        return '\n'.join(outLines)


class TaxHitCluster (MultiDictCluster):
    """
    """
    def __init__ (self,taxdb,*args,**kwds):
        self.taxonomyDB=taxdb
        MultiDictCluster.__init__(self,*args,**kwds)
        #super(TaxHitCluster,self).__init__(*args,**kwds)
        print self.gWeights

    def __add__ (self,other):
        """combine clusters, return copy
        """
        rv=self.__class__(self.taxonomyDB)
        #rv=copy.deepcopy(self)
        for n,w,e in self.iterExpts():
            rv.append(e,n,eWeight=w)
        
        for n,w,e in other.iterExpts():
            rv.append(e,n,eWeight=w)
        return rv

 
    def pTransform(self,controlCluster,
                   minDAp=None,
                   maxSWp=None,
                   missingValue=1.0):
        """
        
        Arguments:
        - `self`:
        - `controlCluster`:
        """
        import scipy.stats as SS
        import scipy as S
        
        cStats = controlCluster.gStats(missingValue)
        
        def pCallback(gID,vDict):
            if gID in cStats:
                if minDAp is not None:
                    DAp = cStats[gID][-1][-1]
                    if minDAp != None and DAp <= minDAp:
                        return None
                    
                if maxSWp is not None:
                    SWp = cStats[gID][-2][-1]
                    if maxSWp != None and SWp > maxSWp:
                        return None

                tList = list(vDict.items())
                keys=[x[0] for x in tList]
                vals=[x[1] for x in tList]
                
                mean=cStats[gID][2]
                std=cStats[gID][4]
                #print mean,std, vals
                ary=SS.norm(mean,std).sf(vals)

                for i,v in enumerate(vals):
                    if v<mean:
                        ary[i] = 1.0
            else:
                return None
                #ary = S.array([S.nan]*len(values))
            return dict(zip(keys,ary))
            
        self.gTransform(pCallback)


    def taxonFilter(self,taxon,decendants=True):
        """
        """
        from ncbi import taxonomy
        taxClass = taxonomy.Taxon
        taxID=int(taxon)
        if type(taxon) == IntType:
            taxon = taxClass(taxID)

        taxIDs=[taxID]
        if decendants:
            taxIDs.extend(taxon.decendantIDs())

        for taxID in taxIDs:
            if taxID in self:
                del self[taxID]
                #del self.gNames[taxID]


        
class M8TaxonFilter (object):
    """Filter m8 results based on taxa of gis in hits.
    """

    def __init__ (self,taxonIDs=None,gis=None,minScore=0,maxE=1000,minLen=0,
                  qryMatch=False,subjMatch=False,ignoreTaxaAndGis=False):
        """Make a new filter.  blast results involving taxonIDs and gis given
        and with scores cuttoffs specified.

        taxa gis are looked up as follows: directly linked gis are used and
        decendants of family genus and species nodes are used.

        TaxonIDs can be for a single taxon or an iterable (Taxon_ID or
        Taxon like objects where int(thing) is a Taxon_ID are allowed.  

        gis are expected in (or to be part of) the subject field or the
        qry field.  set subjMatch or qryMatch to True as appropriate.
        There is no valid default.
        """
        
        from ncbi import giInfo
        self.scoreOnly=ignoreTaxaAndGis
        self.gis=[]

        if (not ignoreTaxaAndGis and (int(qryMatch) + int(subjMatch) != 1)):
            raise ValueError, "Either qrySet or subjSet must be True"
        tGis=set()
        if type(taxonIDs) in (IntType, LongType):
            taxonIDs=(taxonIDs,)
        for tid in taxonIDs:
            if tid != None:
                taxID=int(tid)
                tGis=set(giInfo.GiInfo._table.Gi.valueList(Family_Tax_ID=taxID))
                if len(tGis)==0:
                    tGis=set(giInfo.GiInfo._table.Gi.valueList(
                        Genus_Tax_ID=taxID))
                if len(tGis)==0:
                    tGis=set(giInfo.GiInfo._table.Gi.valueList(
                        Species_Tax_ID=taxID))
                if len(tGis)==0:
                    tGis=set(giInfo.GiInfo._table.Gi.valueList(
                        Tax_ID=taxID))
                if len(tGis)==0:
                    raise ValueError , "NCBI Taxon ID not found: %s" % taxID

            if gis == None:
                gis=[]
            self.gis=set(gis)
            self.gis.update(tGis)


        if len(self.gis)==0 and not self.scoreOnly:
            raise ValueError , "No valid Gis given or specified by taxa"
        

        self.minLen = minLen
        self.minScore=minScore
        self.maxE = maxE
        if subjMatch:
            self.giK = 'subject'
        elif qryMatch:
            self.giK = 'query'
        else:
            raise ValueError, "Either qrySet or subjSet must be True"
        self.subjMatch=subjMatch
        self.qryMatch=qryMatch
        
        
    def filterFile(self,inFile,outFile,append=False,clobber=False):
        """
        """

        if hasattr(outFile,'write'):
            oFile=outFile
            finishFile=oFile.flush
        else:
            oFile = safeOFW(outFile,append=append,clobber=clobber)
            finishFile=oFile.close
            
        for hsp in m8generator(inFile):
            if (hsp['length'] >= self.minLen and
                hsp['score'] >= self.minScore
                and hsp['e'] <= self.maxE):
                hitGi = None
                if not self.scoreOnly:
                    #print hsp[self.giK]
                    if type(hsp[self.giK]) == IntType:
                        # bare int parsing
                        hitGi=hsp[self.giK]
                    else:
                        # ncbi title parsing
                        hitGi = int(hsp[self.giK].split('|')[1])
                if self.scoreOnly or hitGi in self.gis:
                    oFile.write(hsp['_str_'])
                    oFile.write('\n')
        finishFile()


    def extractFastaRecords(self,m8file,fastaFile,outFile=None,clobber=False,append=False):
        """By default, return an iterator of fasta records with titles
        that match filtered m8 records.

        If outFile is not None, all the records returned by the generator
        are dumped to the file subject to clobber and append settings,
        and the number of records is returned.
        """
        
        if outFile != None:
            # do this first to check for error
            fastaOut = safeOFW(outFile,clobber=clobber,append=append)
        
        import fasta
        reload(fasta)
        fOut = StringIO()
        self.filterFile(m8file,fOut)
        fOut.seek(0)
        
        tSet=m8TitleSet(fOut,qrySet=self.subjMatch,
                        subjSet=self.qryMatch)
        gen = fasta.generalIterator(fastaFile,titleSet=tSet)
        if outFile == None:
            return gen
        else:
            c=0
            for rec in gen:
                print >> fastaOut , rec
                c+=1
            return c


class M8scoreCompare(object):
    """
    """
    def __init__(self,mbrPaths=[],
                 fnTransformer=lambda x: x.split('.')[0]):
        self.scores={}
        self.fDone={}
        self.fnTransformer=fnTransformer
        self.addFiles(mbrPaths)


    def addFiles(self,mbrPaths):
        for f in getIterable(mbrPaths):
            k=self.fnTransformer(f)
            self.scores[k]=[]
            for hsp in m8tupleGenerator(f):
                self.scores[k].append(hsp[11])
            if m8complete(f):
                self.fDone[k]=''
            else:
                self.fDone[k]='*'


    def sKeys(self,*sortArgs,**sortKW):
        return sorted(self.scores.keys(*sortArgs,**sortKW))

    def sLabels(self,*sortArgs,**sortKW):
        return [k+self.fDone[k] for k in self.sKeys(*sortArgs,**sortKW)]

    def sValues(self,*sortArgs,**sortKW):
        """values in file name sort order
        """
        return [self.scores[k] for k in self.sKeys(*sortArgs,**sortKW)]
        

    def plotBS(self,outFile='bs.png'):
        """
        """
        from matplotlib import pyplot
        pyplot.figure()
        pyplot.boxplot(self.sValues())
        k=self.sLabels()
        pyplot.xticks(range(1,len(k)+1), k ,rotation=90)
        
        pyplot.title(os.getcwd())
        pyplot.title(os.path.split(os.getcwd()))
        pyplot.title(os.path.split(os.getcwd())[1])
        pyplot.title(os.path.split(os.getcwd())[1]+' HSP Bit Scores')

        pyplot.subplots_adjust(0.125,0.175,0.9)
        pyplot.subplots_adjust(bottom=0.2)
        pyplot.savefig(outFile)
        pyplot.close()

        
    def plotHspCt(self,outFile='HspCt.png'):
        """
        """
        from matplotlib import pyplot
        pyplot.figure()

        pyplot.bar(range(1,len(self.sValues())+1),
                   [len(x) for x in self.sValues()])
        k=self.sLabels()
        pyplot.xticks([i+0.5 for i in range(1,len(k)+1)], k ,rotation=90)
        
        pyplot.title(os.getcwd())
        pyplot.title(os.path.split(os.getcwd()))
        pyplot.title(os.path.split(os.getcwd())[1])
        pyplot.title(os.path.split(os.getcwd())[1]+' HSP Count')
        pyplot.subplots_adjust(0.125,0.175,0.9)
        pyplot.subplots_adjust(bottom=0.2)
        pyplot.savefig(outFile)
        pyplot.close()


class M8Compare (object):

    def __init__():
        
        self.qryIdx={}
        self.subjIdx={}
        

    def addHSPs (self,rsltPath,name=None):
        self.rsltPath=rsltPath
        if name == None:
            self.name
        



class M8Coverage (dict):
    """coverage of gis representated in a set of m8 results as db hits in
    can be plotted or otherwise examined.
    """

    def __init__ (self,m8paths,minScore=0,parseOccCt=True,
                  giFilterSet=None,taxonFilter=None,
                  qryPath=None,dbPath=None,
                  targetFilter=None):

        from ncbi import giInfo

        if type(m8paths) in StringTypes: 
            m8paths = [m8paths]
            
        self._boundsSorted={}
        self.minScore=minScore
        self.qryPath=qryPath
        self.dbPath=dbPath
        self.taxonFilter=taxonFilter
        if targetFilter == None:
            self.targetFilter = lambda x: True
        else:
            self.targetFilter=targetFilter

        self.m8Paths=[]
        
        self.outStub=None
        
        try:
            self.outStub= os.path.splitext(m8paths[0])[0]
        except:
            self.outStub=''
            #raise
            #pass

        gis=[]
        if taxonFilter != None:
            taxID=int(taxonFilter)
            gis=noneDict(giInfo.GiInfo._table.Gi.valueList(Family_Tax_ID=taxID))
            if len(gis)==0:
                gis=noneDict(giInfo.GiInfo._table.Gi.valueList(
                    Genus_Tax_ID=taxID))
            if len(gis)==0:
                gis=noneDict(giInfo.GiInfo._table.Gi.valueList(
                    Species_Tax_ID=taxID))
            if len(gis)==0:
                gis=noneDict(giInfo.GiInfo._table.Gi.valueList(
                    Tax_ID=taxID))
            if len(gis)==0:
                raise ValueError , "NCBI Taxon ID not found: %s" % taxID

        #print gis
        self.m8paths=m8paths
        self.targetGis = {}
        for m8file in fileIterator(self.m8paths):
            #print m8file
            for hsp in m8generator(m8file):
                #print hsp
                if hsp['score'] >= minScore:
                    if type(hsp['subject']) == IntType:
                        hitGi=hsp['subject']
                    else:
                        hitGi = int(hsp['subject'].split('|')[1])
                    if len(gis)==0 or hitGi in gis:
                        try:
                            s,e = (int(hsp['s_start']),int(hsp['s_end']))
                            qs,qe=(int(hsp['q_start']),int(hsp['q_end']))
                            if hitGi not in self:
                                self[hitGi]={}
                                self[hitGi]['_locSorted']=False
                                self[hitGi]['bounds']=[]
                                self[hitGi]['regions']={}
                                self[hitGi]['locQ']=[]
                                self[hitGi]['qryD']={}
                                try:
                                    self[hitGi]['target']=giInfo.GiInfo(Gi=hitGi)
                                except:
                                     self[hitGi]['target']=str(hitGi)
                                if self.targetFilter(self[hitGi]['target']) == False:
                                    del self[hitGi]
                                    continue
                                     
                            if parseOccCt:
                                try:
                                    oCt = int(hsp['query'].split('|')[-1])
                                except:
                                    raise
                            else:
                                oCt = 1
                            
                            self[hitGi]['locQ'].append(((s,e),hsp['query'],(qs,qe)))
                            try:
                                self[hitGi]['qryD'][hsp['query']][0].append((s,e))
                                self[hitGi]['qryD'][hsp['query']][1].append((qs,qe))
                            except KeyError:
                                self[hitGi]['qryD'][hsp['query']]=[[(s,e)],[(qs,qe)]]
                            self[hitGi]['bounds'].extend([(s,e)]*oCt)
                            if (s not in self[hitGi]['regions'] or
                                self[hitGi]['regions'][s] < e):
                                self[hitGi]['regions'][s]=e
                        except:
                            raise
                            pass


    def qrysAtPos(self,gi,pos):
        """return a generator of query tiltes that
        align to gi at position (pos).  Pos can be a number
        or a (start,end) tuple.
        
        Arguments:
        - `self`:
        - `pos`:
        """
        if not self[gi]['_locSorted']:
            self[gi]['locQ'].sort()

        try:
            if len(pos) !=2:
                raise ArgumentError, "pos must be a number or a 2-tuple"
        except TypeError:
            pos=(pos,pos)

        if pos[0]>pos[1]:
            pos = (pos[1],pos[0])
        

        for loc,q,qPos in self[gi]['locQ']:
            s,e = loc
            if s>pos[1]:
                break
            elif e<pos[0]:
                continue
            else:
                yield q,qPos
        


    def qryFastaRecs(self,qrys,trim=False):
        """
        """
        import fasta
        
        qs=set([])
        qLocs={}
        
        for q,l in qrys:    
            qs.add(q)
            if not trim:
                continue
            if q not in qLocs:
                qLocs[q]=[]
            qLocs[q].append(l)
        
        for rec in fasta.FastaIterator(self.qryPath):
            if len(qs) == 0:
                break
            t=rec.title
            if t not in qs:
                continue
            qs.remove(t)
            if not trim:
                yield rec
                continue
            oRec=copy.deepcopy(rec)
            oRec.title+='|'
            for s,e in qLocs[t]:
                rec=copy.deepcopy(oRec)
                if e < s:
                    tmpS=s
                    s=e
                    e=tmpS
                    rec=rec.subrecord(s,e)
                    rec.title+='_rc'
                else:
                    rec=rec.subrecord(s,e)
                yield rec


    def savePlots(self,numberOfGis=None,hspLabel=''):
        """
        """
        from sequence import plot
        reload(plot)
        if hspLabel == '' and hasattr(self,'hspLabel'):
            hspLabel = self.hspLabel

        stuffToPlot = self.items()
        
        if numberOfGis != None:
            stuffToPlot.sort(key=lambda x: (len(x[1]['regions']),
                                            len(x[1]['bounds'])))
            stuffToPlot.reverse()
            
            stuffToPlot=stuffToPlot[:numberOfGis]

        #return None
        for gi,v in stuffToPlot:
            info=v['target']
            bounds= v['bounds']
            outFileName = "%s.%s.png" %(self.outStub,gi)
            print outFileName
            hp=plot.HspPlot(figSize=(600,600))
            hp.setScaffold(gi)
            hp.addHsps(bounds,hspLabel=hspLabel)
            hp.scaleHspHeight(maxH=25)
            hp.savePlot(outFileName)


    def plotGi(self,gi,hspLabel='',figSize=(600,600),
               outFileName=None):
        """
        
        Arguments:
        - `self`:
        - `gi`:
        
        """
        from sequence import plot
        info = self[gi]['target']
        bounds = self[gi]['bounds']
        print bounds
        if outFileName == None:
            outFileName = "%s.%s.png" %(self.outStub,gi)
            
        hp=plot.HspPlot(figSize=(600,600))
        hp.setScaffold(gi)
        hp.addHsps(bounds,hspLabel=hspLabel)
        hp.scaleHspHeight(maxH=25)
        hp.savePlot(outFileName)
    


    def pctCovered(self,gi,minCov=1):
        """return the overall fraction of covered bases for gi
        """
        cVect=self.coverageVector(gi)
        return float(sum(cVect>=minCov))/float(len(cVect))
        
    def coverageVector(self,gi):
        """returns a numpy vector fit fold coverage at each postion
        """
        if gi not in self:
            raise ValueError , "Gi (gi) not found in target set" %s

        l = self[gi]['target'].Length
        rv = N.zeros((l,),dtype=N.uint32)

        for s,e in self[gi]['bounds']:
            rv[s-1:e]+=1
        return rv

    def mostCoveredNt(self):
        """return the githat has the most covered nucleotides
        """

        rv = None
        maxNt = 0
        for gi in self.keys():
            cnt = N.sum(self.coverageVector(gi) >0)
            if cnt >=maxNt:
                rv = gi
                maxNt=cnt
        return rv

        

class M8Overlay  (object):
    pass


M8Plots=M8Coverage


class BlastCluster(MultiDictCluster):
    """ a cluster bulit from m8 results
    """

    def __init__ (self,m8Path):
        utils.MultiDictCluster.__init__(self,missingValue=0)
        for t in m8tupleGenerator(m8Path):
            o,gi=t[:2]
            bs=t[-1]
            if o not in self:
                self[o]={}
            try:
                if bs>self[o][gi]:
                    self[o][gi]=bs
            except KeyError:
                self[o][gi]=bs
        if len(self.eOrder) == 0:
            self.eOrder=list(self.experiments())
 
 


def m8Queries2Fasta( resultPath,queryPath,outputPath,minBitScore=None,
                     append=False,clobber=False):
    """make a fastafile coresponding to the query sequence of hit in
    blast results.
    """
    import fasta
    #outF = file(outputPath,'w')
    outF = safeOFW(outputPath,clobber=clobber,append=append)
    qryI=fasta.FastaIterator(queryPath)

    qryIDs = set()
    for hsp in m8generator(file(resultPath)):
        if minBitScore == None or hsp['score'] >= minBitScore:
            #print hsp['score']
            qryIDs.add(hsp['query'])
            qryIDs.add(str(hsp['subject']))
    print qryIDs


    for q in qryI:
        if q.title in qryIDs:
            qryIDs.remove(q.title)
            print >>outF ,  q
        if len(qryIDs) ==0 :
            break
    outF.close()


def fillGapsInSequences(qSeq, sSeq, hsp):
    segments = hsp['segments']
    lastSend = None
    lastQend = None

    if hsp['strand'] == '-':
        segments.reverse()
        #qSeq = complement(qSeq)

    tempSseq = list(sSeq)
    tempQseq = list(qSeq)

    sSeqPad = 0
    qSeqPad = 0

    for seg in segments:
        if hsp['strand'] == '-':
            s_start = seg['s_end']
            s_end = seg['s_start']
            q_start = seg['q_end']
            q_end = seg['q_start']
        else:
            s_start = seg['s_start']
            s_end = seg['s_end']
            q_start = seg['q_start']
            q_end = seg['q_end']
            
        if lastQend is None and lastSend is None:
            lastQend = q_end
            lastSend = s_end
            continue

        qGadSize = abs(q_start - lastQend) - 1
        sGadSize = abs(s_start - lastSend) - 1

        if qGadSize > sGadSize:
            size = qGadSize - sGadSize
            dashes = '-'*size

            if hsp['strand'] == '-':
                insertionPos = abs(lastSend - hsp['s_end'])
            else:
                insertionPos = abs(lastSend - hsp['s_start'])

            insertionPos += sSeqPad
            tempSseq.insert(insertionPos, dashes)
            sSeqPad += 1
        elif qGadSize < sGadSize:
            size = sGadSize - qGadSize
            dashes = '-'*size

            if hsp['strand'] == '-':
                insertionPos = abs(lastQend - hsp['q_end'])
            else:
                insertionPos = abs(lastQend - hsp['q_start'])

            insertionPos += qSeqPad
            tempQseq.insert(insertionPos, dashes)
            qSeqPad += 1

        lastQend = q_end
        lastSend = s_end

    return (''.join(tempQseq), ''.join(tempSseq))

def m8alignmentEnergies(resultPath, queryPath, dbPath, outputPath,
                        isD1format=False,lazyDict=True):
    """calculate dG for alignmnets in results in resultPath
    add floating point energy values to the end of each m8record.
    """

    import aos,fasta
    from ncbi import giInfo
    reload(giInfo)
    qrySeq=fasta.fastaDictionary(queryPath)
    if lazyDict:
        subSeq=giInfo.GiFastaCache()
    else:
        subSeq=fasta.fastaDictionary(dbPath)

    outF = file(outputPath,'w')

    if isD1format:
        resultsGen = d1generator(resultPath)
    else:
        resultsGen = m8generator(file(resultPath))

    for hsp in resultsGen:
        q_start = hsp['q_start']
        q_end = hsp['q_end']
        s_start = hsp['s_start']
        s_end = hsp['s_end']

        if isD1format:
            if q_start < q_end:
                qSeq=qrySeq[str(hsp['query'])][q_start-1:q_end]
            else:
                qSeq=qrySeq[str(hsp['query'])][q_end-1:q_start]
                
            if  s_start < s_end:    
                sSeq=subSeq[str(hsp['subject'])][s_start-1:s_end]
            else:
                sSeq=subSeq[str(hsp['subject'])][s_end-1:s_start]
                
            #print qSeq, hsp['strand']
            qSeq, sSeq = fillGapsInSequences(qSeq, sSeq, hsp)
            #print qSeq, hsp['strand']
            if hsp['strand'] == '-':
                pass #qSeq=complement(qSeq)
            else:
                sSeq=reverseComplement(sSeq)
            try:
                e=aos.energy(qSeq,sSeq)
                #print e
            except:
                print >> sys.stderr, '>>', qSeq
                print >> sys.stderr, '::', sSeq
                print >> sys.stderr, ';;', hsp
                continue
        else:
            if q_start < q_end:
                qSeq=qrySeq[str(hsp['query'])][q_start-1:q_end]
            else:
                qSeq=qrySeq[str(hsp['query'])][q_end-1:q_start]
                qSeq=reverse(qSeq)

            if  s_start < s_end:
                sSeq=subSeq[str(hsp['subject'])][s_start-1:s_end]
                sSeq=reverse(sSeq)
                sSeq=complement(sSeq)
            else:
                sSeq=subSeq[str(hsp['subject'])][s_end-1:s_start]

                try:            
                    e=aos.blastEnergy(qSeq,sSeq,findGaps=findGaps)
                except Exception, ex:
                    print >> sys.stderr, hsp['_str_']
                    print >> sys.stderr, '>>', qSeq
                    print >> sys.stderr, '::', sSeq
                    print >> sys.stderr, ';;', hsp
                    print >> sys.stderr, str(ex)
                    continue

        print >>outF, hsp['_str_'] + '\t%.2f' % e
    outF.close()


def d1alignmentEnergiesReport(resultPath, queryPath, dbPath, outputPath):
    """calculate dG for alignmnets in results in resultPath
    add floating point energy values to the end of each d1record.
    """

    outF = file(outputPath,'w')

    import aos,fasta
    from ncbi import giInfo
    reload(giInfo)
    
    qrySeq=fasta.fastaDictionary(queryPath)
    subSeq=giInfo.GiFastaCache() #fasta.fastaDictionary(dbPath)

    resultsGen = d1generator(resultPath)

    queryEnergies = {}
    longestSetEnergies = 0
    
    for hsp in resultsGen:
        q_start = hsp['q_start']
        q_end = hsp['q_end']
        s_start = hsp['s_start']
        s_end = hsp['s_end']

        if q_start < q_end:
            qSeq=qrySeq[str(hsp['query'])][q_start-1:q_end]
        else:
            qSeq=qrySeq[str(hsp['query'])][q_end-1:q_start]
            qSeq=reverse(qSeq)

        if  s_start < s_end:
            sSeq=subSeq[int(hsp['subject'])][s_start-1:s_end]
        else:
            sSeq=subSeq[int(hsp['subject'])][s_end-1:s_start]
            sSeq=reverse(sSeq)
            sSeq=complement(sSeq)

        qSeq, sSeq = fillGapsInSequences(qSeq, sSeq, hsp)
                
        try:            
            e=aos.energy(qSeq,reverseComplement(sSeq))
            if hsp['query'] in queryEnergies:
                queryEnergies[hsp['query']].append(e)
            else:
                queryEnergies[hsp['query']] = [e,]

            l = len(queryEnergies[hsp['query']])
            if l > longestSetEnergies:
                    longestSetEnergies = l
        except Exception, ex:
            print >> sys.stderr, '>>', qSeq
            print >> sys.stderr, '::', sSeq
            print >> sys.stderr, ';;', hsp
            print >> sys.stderr, str(ex)
            continue

    queries = sorted(queryEnergies.keys())
    for query in queries:
        energies = sorted(queryEnergies[query])
        energies = map(lambda x: '%.2f' % x, energies)
        energies.extend(['0',]*(longestSetEnergies - len(energies)))
        print >> outF, query, '\t', '\t'.join(energies)
    outF.close()


def m8TitleSet(m8Files,qrySet=False,subjSet=False,
               minScore=0,maxE=100000,minLen=0,minIdent=0):
    """return a set of query titles present in
    one or more blast reslut files
    """
    if type(m8Files) in StringTypes:
        m8Files = [m8Files,]
    if int(qrySet) + int(subjSet) != 1:
        raise ValueError, "Either qrySet or subjSet must be True"
    idxs=[]
    if qrySet:
        idxs.append(0)
    elif subjSet:
        idxs.append(1)

    rv=set()
    for m8t in m8tupleGenerator(m8Files):
        if (minScore > m8t[-1] or m8t[-2] > maxE
            or minLen > m8t[3] or (minIdent-0.001) > m8t[2]):
            continue
        for i in idxs:
            t = m8t[i]
            rv.add(t)
    return rv


def d1SubjRange(alignments):
    """return (min start,max end) for subject in alignments
    """
    s=[]
    e=[]
    for a in getIterable(alignments):
        s.append(a['s_start'])
        e.append(a['s_end'])
    return (min(s),max(e))

def d1QryRange(alignments):
    """return (min start,max end) for queries in alignments
    """
    s=[]
    e=[]
    for a in getIterable(alignments):
        s.append(a['q_start'])
        e.append(a['q_end'])
    return (min(s),max(e))


def ncbiHitIdenticalRegions(hit):
    """
    """
    rv={'subject':[],'query':[]}

        
    for hsp in hit['hsps']:
        for mo in re.finditer('\|+',hsp['midline']):
            s=mo.start()
            e=mo.end()-1
            #print s,e
            if hsp['query-from'] < hsp['query-to']:
                qs=int(hsp['query-from'])-1+s-(
                    hsp['qseq'][:s].count('-'))
                #print hsp['qseq'][:s].count('-')
                qe=int(hsp['query-from'])-1+e-(
                    hsp['qseq'][:s].count('-'))
                
            else:
                # this is the best tested
                qe=hsp['query-to']-1+hsp['query-from']-(1+s) +(
                    hsp['qseq'][:s].count('-'))

                qs=hsp['query-to']-1+hsp['query-from']-(1+e)+(
                    hsp['qseq'][:(1+e)].count('-'))
                # over
            rv['query'].append((qs,qe))    

            if hsp['hit-from'] < hsp['hit-to']:
                ss=int(hsp['hit-from'])-1+s-(
                    hsp['hseq'][:s].count('-'))
                se=int(hsp['hit-from'])-1+e-(
                    hsp['hseq'][:s].count('-'))
                #print hsp['hseq'][:s].count('-')
            else:                
                ss=int(hsp['hit-to'])-(e+1)-(
                    hsp['hseq'][:s].count('-'))
                se=int(hsp['hit-to'])-(s+1)-(
                    hsp['hseq'][:s].count('-'))
            rv['subject'].append((ss,se))
        
        #rv['query'].sort()
        #rv['subject'].sort()
    return rv

def ncbiXmlHitGenerator(xmlFile):
    """xmlFile must be a file handle or a
    string specifing a path to a xmlFile.

    
    
    """
    if type(xmlFile) in StringTypes:
        xmlFile=open(xmlFile)
    
    # an <ITERATION> = 1 query
    for iterStr in ncbiXmlTagGenerator("Iteration",xmlFile):
        iterDOM=xml.dom.minidom.parseString(iterStr).childNodes[0]
        iterRV={}
        #print 'I'

        for iterNode in iterDOM.childNodes:
            #print iterNode.nodeName
            if (iterNode.nodeName.startswith('#') or
                iterNode.nodeName == "Iteration_stats"):
                continue
            elif (iterNode.nodeName != "Iteration_hits"):
                value = str(iterNode.childNodes[0].wholeText)
                try:
                    value=int(value)
                except:
                    pass
                #print value
                iterRV[str(iterNode.nodeName[10:])]=value
                #print iterRV
            else:
                #print iterRV
                # these are the hits
                hits=iterNode.childNodes
                #print hits
                for hit in hits:
                    if hit.nodeName != 'Hit':
                        continue
                    #print hit
                    #print "HIT"
                    rv={'hsps':[]}
                    rv.update(iterRV)
                    for n in hit.childNodes:
                        #reset the rv & hsps
                        #print rv
                        if n.nodeName.startswith('#'):
                            continue
                        elif n.nodeName=="Hit_hsps":
                            for hsp in n.childNodes[1::2]:
                                #print 'HSP'
                                innerRv={}
                                for hspNode in hsp.childNodes:
                                    if not hspNode.nodeName.startswith('#'):
                                        value = str(hspNode.childNodes[0].wholeText)
                                        try:
                                            value=int(value)
                                        except ValueError:
                                            pass
                                        innerRv[str(hspNode.nodeName[4:])]=value
                                rv['hsps'].append(innerRv)
                        else:
                            try:
                                rv[str(n.nodeName[4:])]=str(n.childNodes[0].wholeText)
                            except:
                                #print rv, n.nodeName[4:],n.childNodes[0].wholeText
                                raise 
                    yield rv
                


def ncbiXmlTagGenerator(tagName,fh):
    """
    """
    while True:
        tagStr=ncbiXmlConsume(tagName,fh)
        if tagStr is None:
            raise StopIteration
        yield tagStr


def ncbiXmlConsume(tagName,fh):
    """return the string enclosed by the tagName,
    as read from fh.

    If openging tag is not found return none and
    reset the file's position
    """
    bufferList=[]

    origPos=fh.tell()
    openTagL=['<']+list(tagName)
    closeTagL=['<','/']+list(tagName)+['>']

    oTagPos=0
    cTagPos=0

    #find opening tag
    while bufferList[(-len(openTagL)):] != openTagL:
        c=fh.read(1)
        if len(c)==0:
            fh.seek(origPos)
            return None
        bufferList.append(c)
    bufferList=bufferList[-len(openTagL):]
    #print ''.join(bufferList)

    #find closing tag
    while bufferList[(-len(closeTagL)):] != closeTagL:
        c=fh.read(1)
        if c is '':
            raise ValueError, "end tag (%s) not found." %(''.join(closeTagL))
        bufferList.append(c)


    return ''.join(bufferList)

        
def xmlBlast2energy(xmlFile):
    """return generator of (qry-def,subj-def,hsp,energy))
    """
    from aos import energy
    for hit in ncbiXmlHitGenerator(xmlFile):
        for hsp in hit['hsps']:
            qSeq=hsp['qseq']
            sSeq=reverseComplement(hsp['hseq'])
            yield (hit['query-def'],hit['def'],hsp,energy(qSeq,sSeq))
            
        
    
