# solexa.py
# Solexa sequencing classes and functions
#
# by Kael Fischer 
#
# 2007-
#
__version__ = "$Revision: 1.24 $"
# $Source: /r1/cvs/sequence/solexa.py,v $


import types
import warnings
import re
import types
import pickle
import os.path
import tabbedfile
import compact
import fasta
import commands
import shutil
import __init__ as SEQ

ElandUniqueMatchTypes = [ 'U0','U1','U2']
ElandRedundantMatchTypes = [ 'R0','R1','R2']
ElandMatchingTypes = ElandUniqueMatchTypes+ElandRedundantMatchTypes
ElandNoMatchTypes = ['NM','QC']
AllElandTypes=ElandMatchingTypes+ElandNoMatchTypes



class SolexaRead(fasta.Record):
    """Container for a Solexa read.
    Based on the fasta.Record class.
    """
    pass



def parseBustard(fileNames, dropDots=True):
    """A generator for reads in Bustard *_seq.txt files that
    returns SolexaRead objects
    """
    if type(fileNames) not in (types.ListType,types.TupleType):
        fileNames = (fileNames,)

    for fn in fileNames:
        f = file(fn)
        for line in f:
            try:
                (undef,undef,undef,undef,seq) = line.strip().split()
            except:
                warnings.warn("Bustard line '%s' did not parse\n" % line)
                continue
            try:
                # sro = SolexaRead(seq[1:])
                sro = seq[1:]
            except:
                warnings.warn("Bustard seq '%s' does not make a SolexaRead\n" % seq)
                continue
            yield sro
            
            
def parseGerald( fileNames, fastAndLoose=False ):
    """Given a list of the names of Gerald output files
    (s_N_sequence.txt) or a single filename,
    returns a list of the reads as strings.

    fastAndLoose = True will result in a very fast
    parser with no error checking. Use with caution."""

    if type(fileNames) not in (types.ListType,types.TupleType):
        fileNames = (fileNames,)

    sequenceLine = False
    qualityLine = False
    sequenceLineRE = re.compile("\@HWI.+")
    qualityLineRE = re.compile("\+HWI.+")
    line = 0
    for fn in fileNames:
        f = file(fn)
        for line in f:
            if fastAndLoose:
                if line % 4 == 1:
                    yield( line.strip() )
                line = line + 1
            else:
                if sequenceLineRE.match( line ):
                    sequenceLine = True
                elif qualityLineRE.match( line ):
                    qualityLine = True
                else:
                    if sequenceLine:
                        yield( line.strip() )
                    sequenceLine = False
                    qualityLine = False

class ElandRecord (object):
    """Represents a single line from an eland output
    file."""

    ##TODO:
    ## add __setattr__
    ## and __repr__ (should return txt that is eland compliant)
    ##

    fieldOrder = [ "title",
                   "sequence",
                   "match",
                   "numExact",
                   "num1Err",
                   "num2Err",
                   "genomeFile",
                   "genomePosn",
                   "direction",
                   "nHandling",
                   "firstError",
                   "secondError" ]

    def __init__( self, line ):
        """Parses the given line and stores the resulting data
        in self."""

        data = line.strip().lstrip('>').split()
        self.data = {}

        for (i, elt) in enumerate( data ):
            self.data[ self.fieldOrder[i] ] = elt

    def __getattr__( self, item ):
        try:
            return self.data[item]
        except:
            raise AttributeError

    def __getitem__( self, item ):
        try:
            return self.data[item]
        except:
            raise IndexError


    def isMatch( self ):
        """Returns True iff this record represents at least
        one match."""
        if self.match == "QC":
            return False
        return sum( map( int, [self.numExact, self.num1Err, self.num2Err] ) ) > 0

    def isSingleMatch( self ):
        """Returns True iff this record represents a single
        match to a target genome."""

        return self.match.startswith("U")


    def __str__( self ):
        """Prints out all internal data"""
        f=[]
        for i in range(len(self.data)):
            f.append(self.data[self.fieldOrder[i]])
        return '>' + '\t'.join(f)

    def fasta(self):
        """returns fasta.Record corresponding to this record
        """
        import fasta
        return fasta.Record(title=self.title,
                            sequence=self.sequence)
    def LZWsize(self):
        """return the size of the LZW compressed sequence
        """
        import aos
        return aos.LZWsize(self.sequence)

    def LZWratio(self):
        """Return the LZW compression ratio.
        """
        return float(self.LZWsize())/float(len(self.data['sequence']))

    def Ux(self):
        """True if the record is match type U0,U1 or U2.
        """
        return self.match in ('U0','U1','U2')
        

    def Rx(self):
        """True if the record is match type R0,R1 or R2.
        """
        return self.match in ('R0','R1','R2')

    def isMatch( self ):
        """Returns True iff this record represents at least
        one match."""
        if self.match == "QC":
            return False
        return sum( map( int, [self.numExact, self.num1Err, self.num2Err] ) ) > 0
    
    def isSingleMatch( self ):
        """Returns True iff this record represents a single
        match to a target genome."""
        
        return self.match.startswith("U")


class ElandSummary (object):
    """Container for Eland Processing.
    Parses Eland results file and calculates summary stats and
    can filter eland records, or process records selectively for
    fasta output (and further elanding).

    (needs to be rewritten for speed - grep encapsulation probably)

    """

    def __init__(self,fileNames,name=None):
        """
        """
        self.seqs = {} # sequence:dupCount
        self.rsltCounts = {'NM':0,
                           'QC':0,
                           'R0':0,
                           'R1':0,
                           'R2':0,
                           'U0':0,
                           'U1':0,
                           'U2':0
                           }
        
        self.noDupRsltCounts = {'NM':0,
                                'QC':0,
                                'R0':0,
                                'R1':0,
                                'R2':0,
                                'U0':0,
                                'U1':0,
                                'U2':0
                                }


        if type(fileNames) in types.StringTypes:
            fileNames=[fileNames]
        
        for er in parseEland(fileNames):
            s=er.sequence
            code=er.match
            self.rsltCounts[code]+=1
            if s not in self.seqs:
                self.noDupRsltCounts[code]+=1
                self.seqs[s]=1
            else:
                self.seqs[s]+=1

        if name==None:
            self.name = os.path.split(fileNames[0])[1]
        else:
            self.name = name

        #self.mostPopularSeq = self.seqRecords.values()
        #self.mostPopularSeq.sort(key=lambda n: n.dupCount)
        #self.mostPopularSeq.reverse()



    def sequenceCount(self):
        """Return number of distinct sequences.
        """
        return len(self.seqs)

    def recordCount(self):
        """Total number of records.
        """
        return sum(self.seqs.values())

    def __len__(self):
        """Total number of records.
        """
        return self.recordCount()

    def matchSummary(self):
        """Returns a string with the match types counts and totals.  
        """

        mat = [['','records','sequences','repeats','%repeats']]
        rc=self.recordCount()
        sc=self.sequenceCount()
        repeatCount=0
        
        mts = ('NM', 'QC','R0','R1','R2','U0','U1','U2')
        for mt in mts:
            mat.append([mt, self.rsltCounts[mt],self.noDupRsltCounts[mt],
                        self.rsltCounts[mt]-self.noDupRsltCounts[mt],
                        float(self.rsltCounts[mt]-self.noDupRsltCounts[mt])/self.rsltCounts[mt]])
            repeatCount+=self.rsltCounts[mt]-self.noDupRsltCounts[mt]
            
        mat.append(['Total',rc,sc,
                    repeatCount,float(repeatCount)/rc])

        for mt in mts:
            mat.append(['%'+mt, float(self.rsltCounts[mt])/rc,
                        float(self.noDupRsltCounts[mt])/sc])

        return tabbedfile.tabbifyMatrix(mat)


class ElandComparision (object):

    def __init__ (self,*summaries):
        self.summaries = summaries




    def plotSummary(self, fileName=None,
                    title="Eland Match Summary",
                    xLabel=""):
        import pylab as P
        import matplotlib

        matchLabels = ('match U0',
                  '      U1',
                  '      U2',
                  '      R0',
                  '      R1',
                  '      R2',
                  'no match',
                  'QC')

        
        P.figure()

        bWidth = 0.25
        
        matchPos = [0.2 + i for i in range(len(self.summaries))]            
        nomatchPos = [0.55 +i for i in range(len(self.summaries))]

        U0=P.array([float(s.rsltCounts['U0']) for s in self.summaries])
        U1=P.array([float(s.rsltCounts['U1']) for s in self.summaries])
        U2=P.array([float(s.rsltCounts['U2']) for s in self.summaries])
        R0=P.array([float(s.rsltCounts['R0']) for s in self.summaries])
        R1=P.array([float(s.rsltCounts['R1']) for s in self.summaries])
        R2=P.array([float(s.rsltCounts['R2']) for s in self.summaries])

        hackBar=P.bar([0,len(self.summaries)*1.35],[0.0,0.0],0.2)

        mBarU0=P.bar(matchPos,U0, width=bWidth,color=(0.0,0.0,1.0))
        mBarU1=P.bar(matchPos,U1 ,width=bWidth,color=(0.3,0.3,1.0),bottom=U0)
        mBarU2=P.bar(matchPos,U2 ,width=bWidth,color=(0.6,0.6,1.0),bottom=U0+U1)
        mBarR0=P.bar(matchPos,R0 ,width=bWidth,color=(0.0,1.0,0.0),bottom=U0+U1+U2)
        mBarR1=P.bar(matchPos,R1 ,width=bWidth,color=(0.3,1.0,0.3),bottom=U0+U1+U2+R0)
        mBarR2=P.bar(matchPos,R2 ,width=bWidth,color=(0.6,1.0,0.6),bottom=U0+U1+U2+R0+R1)

        NM=P.array([float(s.rsltCounts['NM']) for s in self.summaries])
        QC=P.array([float(s.rsltCounts['QC']) for s in self.summaries])
        
        mBarNM=P.bar(nomatchPos,NM ,width=bWidth,color='red')
        mBarQC=P.bar(nomatchPos,QC ,width=bWidth,color='y',bottom=NM)
        

        P.legend((mBarU0[0],mBarU1[0],mBarU2[0],
                  mBarR0[0],mBarR1[0],mBarR2[0],
                  mBarNM[0],
                  mBarQC[0]),matchLabels,
                 prop = matplotlib.font_manager.FontProperties(size='smaller'))

        
        P.title(title)
        P.xlabel(xLabel)

        P.ylabel('read count')
        
        P.xticks((P.array(matchPos)+P.array(nomatchPos)+bWidth)/2.0,[s.name for s in self.summaries])
        if fileName == None:
            P.show()
        else:
            P.savefig(fileName)
        
        
        


def parseEland( fileNames ):
    """Given one or more eland output files, parses
    them and returns a list of ElandRecords."""
    
    if type(fileNames) not in (types.ListType,types.TupleType):
        fileNames = (fileNames,)

    for fn in fileNames:
        f = file(fn)
        for line in f:
            try:
                yield ElandRecord( line )
            except ValueError:
                warnings.warn("Line Unparseable: \n%s\n" % line )

def eland2fastaString (line):
    """Take an eland string and return a fasta string
    Can be a generator or sequence of eland lines
    """
    if type(line) in (types.ListType,types.TupleType,types.FileType):
        return [eland2fastaString(l) for l in line]
    
    title,seq,undef = line.split(None,2)
    return "%s\n%s\n\n" % (title,seq)


class MultiSequenceSquashFile(fasta.Record):
    """An object that handles going back and forth between
    multi-record Fasta files and cat'ed single sequence fasta
    representations sutible for eland squashed formatting
    """

    def __init__(self,fastaOutput,fastaInput=None,title='',
                 clobber=False,linkerLength=50,linkerChar='N',colwidth=60):
        """Create instance, optionally with a fasta
        file object or path.  
        """
        if not clobber and os.path.exists(fastaOutput):
            raise ValueError , "%s exists - beep!!" % fastaOutput

        if title == '':
            title = 'MSSquash'

        if type(fastaInput) in types.StringTypes:
            title += ('_'+os.path.split(fastaInput)[1])

        
        fasta.Record.__init__(self,title=title,colwidth=colwidth)
        self.fastaOutput = fastaOutput
        self.recOffsets=[]  # recOffsets are the starting pos of each sequence
        self._lastRecordLen = 0
        self.recTitles=[]
        self.linker=linkerChar*linkerLength
        

        print >> file(fastaOutput,'w'), ">%s" % title 

        print len(list(fasta.FastaIterator(fastaInput)))

        if fastaInput != None:
            self.appendFastaRecords(fasta.FastaIterator(fastaInput))

    def appendFastaRecords(self,records):
        """Add fasta records from some iterable of such records. 
        The records must have a '.sequence' attribute, such as
        fasta.Record. See also fasta.FastaIterator.
        """

        outFile=file(self.fastaOutput,'a')

        try:
            for rec in records:
            
                self.recTitles.append(rec.title)
                if len(self.recOffsets) == 0:
                    self.recOffsets.append(0)
                else:
                    print >>outFile, self.linker
                    ##TODO: Test RecOffsets
                    self.recOffsets.append(self.recOffsets[-1] 
                                           + self._lastRecordLen  
                                           + len(self.linker)) 
                self._lastRecordLen=len(rec)
                print >>outFile, rec.wrappedSequence()
                print self.recOffsets[-1]
                
        finally:
            outFile.close()
        

    def savePickle(self,pklName):
        """save me out the the file system.  reload later with
        loadMSSF.
        """
        pickle.dump(self,file(pklName,"w"),0)


    def fixEland(self,elandRec):
        """map eland line to the correct sequence record and offset.
        return corrected ElandRecord.
        """
        if not isinstance(elandRec,ElandRecord):
            raise ValueError, "input to fixEland is not an ElandRecord instance."
        if elandRec.genomeFile != self.title and False:
            #TODO: that's  broken dude
            return elandRec
        else:
            eOffset = int(elandRec.genomePosn)
            print eOffset
            for i in range(len(self.recOffsets)):
                if self.recOffsets[i] > eOffset:
                    correctOffset = eOffset-self.recOffsets[i-1]
                    correctTitle = self.recTitles[i-1]
                    break
            elandRec.data['genomeFile']=correctTitle
            elandRec.data['genomePosn']=correctOffset
            return elandRec
        

def loadMSSF(pklName):
    """restore a MultiSequenceSquashFile instance from the disk.
    """
    return pickle.load(file(pklName))


def laneString (l):
    """returns string: s_<l>
    where l is an integer type
    """
    return "s_%d" % l

def tileString (l,t):
    """returns string: s_l_t
    where l and t are integer types.
    string rep. of t is padded with leading zeros to 4 digits.
    """
    return "s_%d_%04d" %(l,t)

def extractNM( elandResults, outputPath=None, forceLowMemory=False ):
    """Given the path to one or more eland results files (as a list or string),
    generates FASTA data for every read which has no matches (NM) in *ALL* of the
    eland results files.

    If an output file name is specified, the FASTA data is written to file, and
    the entire process takes very little memory. (it is all run using shell commands)
    The ( status, output ) of the job are returned, and the output should be empty on
    success.

    With no output file specified, the FASTA data are returned as the output, so they will
    be in memory. Use with caution for large data sets. Once again the return type is
    (status, output)

    Setting forceLowMemory to True uses a different shell command which will run slower
    but use less memory. In a case of 3 eland files with 6M reads each, the memory use
    was halved while the run time doubled.
    """

    #if type(elandResults) != type([]):  #That is wrong b/c it doesn't work for tuples & generators
    if type(elandResults) in types.StringTypes: 
        elandResults = [ elandResults ]         # protect strings from list-ification
    else:
        elandResults = list(elandResults)       # handles generators/iterables
        

    numFiles = len(elandResults)

    if forceLowMemory:
        cmd = (
            """sort -k 1 %s | """
            """awk 'BEGIN {OFS=""} """
            """/NM/ { if( $1 == lastQ ) count++; else count = 1; """
            """if( count >= numFiles ) print $1,"\n",$2,"\n"; lastQ = $1 }' numFiles=%d"""
            % ( " ".join( elandResults ), numFiles ))
    else:
        cmd = (
            """awk 'BEGIN {OFS=""} """
            """/NM/ { missed[$1]++; if( missed[$1] == numFiles ) print $1,"\\n",$2,"\\n" }' numFiles=%d %s"""
            % ( numFiles, " ".join( elandResults ) ))
    
    if outputPath:
        cmd += " > %s" % (outputPath)

    return commands.getstatusoutput( cmd )


def postProcessEland(fileNames,allowMatchTypes=AllElandTypes,
                     minLZWsize=0,minLZWratio=0.0,
                     removeRepeats=False,
                     regexRequire=None, regexRequireName='regex',
                     regexProhibit=None, regexProhibitName='not_regex'):
    """ Perform several postprocessing/filtering steps on one or more Eland results files.
    outputfiles are named by inserting the filtering options before the extention (of any)
    of the input file name.

    If a regex (compared to sequence) is required or prohibited, a descriptive regexRequireName
    and/or regexProhibitName, can and should be given as well.

    Possible opperations:
    -Filter on LZW parameters
    -filter on matchtype
    -remove duplicate records
    -require sequence regex match
    -require sequence regex no match
    """
    
    if type(fileNames) in types.StringTypes:
        fileNames = [fileNames]

    if type(allowMatchTypes) in types.StringTypes:
        allowMatchTypes = [allowMatchTypes]

    for fn in fileNames:
        filesToClean=[]
        try:
            if removeRepeats:
                wrkFile, wrkFn = SEQ.mystemp()
                outFile, outFn = SEQ.mystemp()

                os.system("sort %s > %s" %(fn,wrkFn))
                diff = bool(commands.getstatusoutput("diff -q %s %s" % (fn,wrkFn))[0])

                if not diff:
                    # already non-redundant and sorted!
                    os.unlink(wrkFn)
                    os.unlink(outFn)
                    compFn = fn

                wrkFile=file(wrkFn)
                filesToClean.append(outFn)

                lCt=0
                lastLine = None
                lastSeq = None
                for l in wrkFile:
                    seq = l.split()[1]
                    if seq != lastSeq:
                        if lastLine == None: 
                            lastLine=l
                        else:
                            er = ElandRecord(lastLine)
                            er.data['title']+='|%s'%(lCt+1)
                            print >>outFile, er 

                            lastLine = l
                            lastSeq = seq
                            lCt=0
                    else:
                        lCt+=1
                    
                lCt+=1
                er = ElandRecord(lastLine)
                er.data['title'] +='|%s'%lCt
                print >>outFile, er
                wrkFile.close()
                outFile.close()
                os.unlink(wrkFn)

                compFn = outFn
            else:
                compFn = fn

            # repeat removal done
            
            # calc ulate output file name
            fnParts = os.path.splitext(fn)
            outFn = fnParts[0]

            if allowMatchTypes != AllElandTypes:
                if allowMatchTypes==ElandMatchingTypes:
                    outFn += ".AllMatch"
                else:
                    outFn+='.'+'_'.join(allowMatchTypes)

            if minLZWsize >0:
                outFn+='.LZW%s' % minLZWsize

            if minLZWratio >0:
                outFn+='.LZWr%s' % minLZWratio

            if removeRepeats:
                outFn+='.noDups'

            if regexRequire != None:
                outFn+= '.'+regexRequireName

            if regexProhibit != None:
                outFn+= '.'+regexProhibitName

            

            outFn+=fnParts[1]

            # if there is nothing more to be done
            # and repeats have been removed copy
            # that and go to the next file
            if (allowMatchTypes==AllElandTypes and minLZWsize == 0 and minLZWratio == 0 and
                regexRequire == None and regexProhibit == None):
                if outFn != fn:
                    if compFn != fn:
                        shutil.move(compFn,outFn)
                        filesToClean.remove(compFn)
                        continue
                    
            # filter line by line
            outFile = file(outFn, 'w')
            for l in file(compFn):
                er=ElandRecord(l)

                # filters
                if er.match not in allowMatchTypes:
                    continue
                if minLZWsize>0:
                    if er.LZWsize() < minLZWsize:
                        continue
                if minLZWratio >0:
                    if er.LZWratio() < minLZWratio:
                        continue
                if regexRequire != None:
                    if regexRequire.search(er.sequence) == None:
                        continue
                if regexProhibit != None:
                    if regexProhibit.search(er.sequence) != None:
                        continue
                
                
                # don't bother with the rec
                # just dmp the line
                outFile.write(l)
            outFile.close()

        finally:
            for fn in filesToClean:
                try:
                    os.unlink(fn)
                except:
                    pass


class PathogenDetector (object):
    """Rolls up components of a 'pipeline' that searches for pathogenic sequences
    in a solexa sequencing run.

    Gerald fastq output is required.
    """

    def __init__ (self, fastqPath,PDspec=None, maxGridJobs=5):
        """Make a PD instance.  You have to do this first.
        one argument is required, the fastqPath, which should
        point to the directory with the gerald fastq output 

        PDspec, if specified, indicates a subdir with the name PD_<PDspec>.
        If the subdir exists processing in there will continue, otherwise a new
        working subdir (of fastqpath) will be created. 

        """

        #
        # Note all jobs should be submitted to the grid 
        #
        

        self.fastqPath=fastqPath
        #self.inputFiles=
        # store some stat info about infiles too

        self.primerSequences=[]
        self.primerScreenComplete=False

        self.minQuality=0
        self.maxLowQualityBases=0
        self.qualityFilterd=False

        self.barcodeLength=0
        self.barcodeSplit=False

        self.retainDupes=False
        self.dupsRemoved=False

        self.maxHomo=None
        self.homoFiltered=False

        self.maxLZW=None
        self.LZWfiltered=False

        self.processBarcodes=None # if None process all

        self.preScreens=None # sequence of Screen Objects (human, bacto, etc)
        self.completedScreenIdx=[]

        # the search objects contain the scoring options
        # there can be more than one search per PD and more that one
        # scoring per search
        self.searches=None
        

    def saveState(self):
        """Pickle myself in PDspec subdir so processing can continue later 
        """
        pass


    def primerScreen(self, primerSequences):
        """Remove primer sequences from  
        """
        pass

    


class Search (object):
    """Object specifying search to be performed, and scoring to be applied
    """
    def __init__ (self):
        pass

    def start(self):
        """Start Job
        """
        pass

    def done(self):
        """Jobe done?
        """
        pass


    def hits(self,scoreIdx=0,**kwargs):
        """Return the scores.  The details of what is a hit and returned structure
        are by the ScoreMethod referenced.  kwargs are passed through to
        the ScoreMethod's report method.
        """
        return self.scores[scoreIdx].hits(self,**kwargs)

    def removeHitsOnDisk(self,scoreIdx=0,**kwargs):
        """remove hits from source file.
        """
        self.scores[scoreIdx].removeHitsOnDisk(self,**kwargs)

    def report(self,scoreIdx=0,**kwargs):
        """Return the scores.  The details of the score and returned structure
        specified by the ScoreMethod referenced.  kwargs are passed through to
        the ScoreMethod's report method.
        """
        return self.scores[scoreIdx].report(self,**kwargs)

    def plot(self,scoreIdx=0, **kwargs):
        """Returns a plot of the scores.  The details of the score and returned thing
        (plot object, indication of success on plot IO, etc.) are specified by the
        ScoreMethod referenced. kwargs are passed through to
        the ScoreMethod's plot method.
        """
        return self.scores[scoreIdx].plot(self,**kwargs)
        
        
    
    

            
class Screen (Search):
    """object specifying a search db and method and scoring to remove reads
    from a file
    """
    pass

class ScoreMethod (object):
    """Protocol for scoring search results
    """
    pass

    
    
