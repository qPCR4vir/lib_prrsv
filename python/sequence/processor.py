## #!/usr/local/bin/python 
## #
## # processor.py
## #
## # TaxHitDistribution creation
## #

## __verison__ = tuple([int(x) for x in
##                      '$Revision: 1.11 $'.split()[1].split('.')])
## __author__ = "Julio Menendez, Kael Fischer"
## import os, os.path
## import subprocess
## import threading
## import re
## import ConfigParser
## import cPickle
## import datetime

## from sequence.fasta import m8partition, bowtiePartition
## from sequence.blastNoSQL import TaxHitDistribution, TaxHitDistributionSet
## from grid import GridSubmission, waitForJob
## from viroinfo.megachip import * #BlastFilterProtocol, \
## #     AlignmentProtocol, THD, TaxonCount, ReadCount, ReadCountType, \
## #     THDAnalysis, THDS
## from utils import mystemp
## from ncbi.taxonomy import Taxon
## import ncbi


## def runCommand(command):
##     """Runs the `command` and returns the output.
    
##     Arguments:
##     - `command`: List with the command as first item and the arguments
##         as the rest of the items.
##     """
##     return subprocess.Popen(command,
##                             stdout=subprocess.PIPE).communicate()[0]
    

## def applyFilterProtocol(inputFilename, filterProtocolID, outputFilename=None):
##     """Takes the `inputFilename` file and apply the filter defined by
##     `filterProtocolID` and save the output in `outputFilename`. Notice
##     that the filtering can return more than 1 file so the
##     `outputFilename` will contain numbers before the extension.
    
##     Arguments:
##     - `inputFilename`: Fasta file to filter.
##     - `filterProtocolID`: Blast_Filter_Protocol_ID from
##         MegaChipDB.Blast_Filter_Protocol.
##     - `outputFilename`: File where the output will be stored. If None
##         the output will be saved using the name of the filter protocol
##         used and the `inputFilename`.
##     """
##     filterProtocol = BlastFilterProtocol(int(filterProtocolID))
##     print filterProtocol.SQLfields
##     #stepList = map(lambda x: x.SQLfields, filterProtocol.stepList())
##     files = [(inputFilename, 0),]


##     #print os.environ
    
    
##     for stepPos,step in  enumerate(filterProtocol.steps()):
        
##         keepHits = step.Keep_Hits is not None
##         keepMisses = step.Keep_Misses is not None
##         stepName = step.Name
##         print 'stepName',stepName
##         stepDbPath = step.Database_Path
##         stepClParams = step.CL_Parameters
##         alignmentProgram = step.Alignment_Program
        
##         inFile, level = files.pop(0)
##         outFilePreffix = os.path.splitext(os.path.split(inFile)[1])[0]

##         print 'beep'
##         print outFilePreffix
##         outMatchFile, outNoMatchFile = None, None
        
##         if keepHits:
##             outMatchFile = '.'.join((outFilePreffix,
##                                      stepName,
##                                      'Match',
##                                      'fasta'))
##             print outMatchFile
##             outMatchFile = mystemp(dir='/tmp',suffix=outMatchFile)[1]
##             print outMatchFile
##             files.append((outMatchFile, stepPos + 1))
##         if keepMisses:
##             outNoMatchFile = '.'.join((outFilePreffix,
##                                        stepName,
##                                        'NoMatch',
##                                        'fasta'))
##             outNoMatchFile = mystemp(dir='/tmp',suffix=outNoMatchFile)[1]
##             files.append((outNoMatchFile, stepPos + 1))

##         if alignmentProgram == 'megablast':
##             mbrFile = mystemp(dir='/tmp',suffix='.mbr')[1]
##             cmdOut = runCommand(['megablast',
##                                  '-i%s' % (inFile,),
##                                  '-d%s' % (stepDbPath,),
##                                  '-o%s' % (mbrFile,),
##                                  stepClParams,
##                                  '-D3'])
##             m8partition(inFile, outMatchFile, outNoMatchFile, mbrFile,
##                         qryMatch=True, clobber=True)

##             os.remove(mbrFile)
##             if inFile != inputFilename:
##                 os.remove(inFile)

##             if files[0][1] > stepPos:
##                 stepPos += 1
##         elif alignmentProgram == 'bowtie':
##             btFile = mystemp(dir='/tmp',suffix='.bt')[1]
##             cmdOut = runCommand(['bowtie',
##                                  '-f',
##                                  stepClParams,
##                                  stepDbPath,
##                                  inFile,
##                                  btFile,])

##             bowtiePartition(inFile, outMatchFile, outNoMatchFile, btFile,
##                             qryMatch=True, clobber=True)

##             os.remove(btFile)
##             if inFile != inputFilename:
##                 os.remove(inFile)

##             if files[0][1] > stepPos:
##                 stepPos += 1
##         else:
##             raise RuntimeError, (
##                 'Blast Filter Step, unknown Alignment_Program: "%s"' \
##                 % (alignmentProgram))

##     if outputFilename is None:
##         outputFilename = os.path.splitext(inputFilename)[0]
##     print outputFilename
##     finalFiles = []
##     for i, f in enumerate(files):
##         fname = '.'.join(f[0].split('.')[1:])
##         print fname
##         finalFilename = '.'.join((outputFilename,fname))
##         fOut = open(finalFilename, 'w')
##         fOut.write(open(f[0]).read())
##         fOut.close()
##         finalFiles.append(finalFilename)
##         os.remove(f[0])
##     return finalFiles


## def countFasta(fastaFile):
##     """Counts the reads and sequences in `fastaFile`. `fastaFile` can
##     be a .fasta or a .fastq file.
    
##     Arguments:
##     - `fastaFile`: .fasta or .fastq file with read information for
##     each sequence in it.
##     """
##     command = ('awk',
##                ('BEGIN {FS="|"; OFS=":" ;rct=0;sct=0;} /^@|^>/ '
##                '{ rct+=$2;sct+=1} END {print rct,sct}'),
##                fastaFile)
##     print command
##     output = runCommand(command)
##     return tuple(map(int, output.strip().split(':')))    


## FASTA_FILTER_PATH = '/r1/lib/python/scripts/fastaFilter.py'


## class THDCreator(object):
##     """Given a sequence file will run all the specified filterings and
##     alignment protocols and store in DB the taxon count of the
##     TaxonHitDistribution object resulting.
##     """
    

##     def __init__(self, sequenceFileID,
##                  alignmentProtocolID,
##                  scoreCutoff,
##                  taxonomyDatabaseID,
##                  filterProtocolID,
##                  clobber=False):
##         """Initialize the instance.
        
##         Arguments:
##         - `sequenceFileID`: Sequence_File_ID from
##             MegaChipDB.Sequence_File
##         - `alignmentProtocols`: Alignment_Protocol_ID from
##             MegaChip.Alignment_Protocol to be used to do the
##             alignment.
##         - `taxonomyDatabaseID`: Taxonomy_ID from MegaChipDB.Taxonomy.
##         - `filterProtocolID`: Blast_Filter_Protocol_ID from
##             MegaChipDB.Blast_Filter_Protocol table. If None no
##             filter will be used.
##         - `clobber`: If True will overwrite the THD with this same
##             information.
##         """
##         self._sequenceFileID = sequenceFileID
##         self._alignmentProtocolID = alignmentProtocolID
##         self._scoreCutoff = scoreCutoff
##         self._taxonomyDBID = taxonomyDatabaseID
##         self._filterProtocolID = filterProtocolID
##         self._clobber = clobber

##         # TODO: Change this!!
##         #self._filesPreffix = '-'.join((str(sampleID),
##         #                               str(sequencingRunID)))
##         self._threads = []

##     def sequenceFile(self):
##         return SequenceFile(self._sequenceFileID)

##     def filePrefix (self):
##         return os.path.splitext(str(self.sequenceFile()))[0]

##     def _alignmentDoneSuccesful(self, resultFilename):
##         """Checks if `resultFilename` exists and checks if it contains the
##         marker of succesful finish of the alignment program.
        
##         Arguments:
##         - `resultFilename`: filename used as output of the alignment program.
##         """
##         if not os.path.exists(resultFilename):
##             return False

##         lastLine = runCommand(('tail', '-n1', resultFilename))
##         matches = re.match(r'^#(?:.*)finished(?:.*) (\d+) (?:.*)',
##                            lastLine)
##         return matches is not None

##     def filterQuality(self, inputFilename):
##         """Filters the sample based on the quality. Returns the
##         filename generated.
##         """
##         prefix = os.path.splitext(inputFilename)[0]
##         output = '.'.join((prefix, 'maxH18', 'LZWr0.4',
##                            'Q5_13Pos', 'fasta'))
##         if not os.path.exists(output):
##             runCommand((FASTA_FILTER_PATH, '-r0.4', '-H18', '-q5:13',
##                         inputFilename))
##         return output

##     def removeDuplicates(self, inputFilename):
##         """Removes duplicates from a sequence file
##         """
##         output = '.'.join((self.filePrefix(), 'noDups', 'fastq'))
##         if not os.path.exists(output):
##             runCommand((FASTA_FILTER_PATH, '-Q', '-d', inputFilename))
##         return output
        

##     def applyFilterProtocol(self, inputFilename):
##         """Applies the filter protocol specified to the sample.
##         """
##         return applyFilterProtocol(inputFilename,
##                                    self._filterProtocolID,)
##                                    #self.filePrefix()) 

##     def startAlignmentJobs(self, inputFilenames):
##         """Starts the alignment jobs on the grid using threads to
##         monitor them waiting the jobs to finish to create the THD objects.
        
##         Arguments:
##         - `inputFilenames`: List of filenames to use as input for the
##             alignment jobs.
##         """
##         apObj = AlignmentProtocol(self._alignmentProtocolID)
##         args = apObj.Arguments % {'database': apObj.Database_Path}
##         command = ' '.join((apObj.Program, args))
##         for idx, inFilename in enumerate(inputFilenames):
##             outFilename = '.'.join((self.filePrefix(),
##                                     apObj.Name, str(idx), 'mbr'))
##             print outFilename
##             thread = threading.Thread(target=self._doAlignment,
##                                       args=(command, inFilename,
##                                             outFilename, idx))
##             self._threads.append(thread)
##             thread.start()

##     def waitForAlignmentJobs(self):
##         """Waits for all the alignment jobs running in threads to
##         finish. It's only a way to block the main thread until every
##         child thread have finished
##         """
##         isAlive = lambda x: x.isAlive()
##         while any(map(isAlive, self._threads)):
##             continue

##     def _doAlignment(self, command, inputFilename, outputFilename, index):
##         """Runs the alignment `command` using `inputFilename` as input
##         and `outputFilename` as output and creates the THD object from
##         the result and stores it in the DB.
        
##         Arguments:
##         - `command`: Command to execute for the alignment
##         process.
##         - `inputFilename`: Input filename for the alignment
##         process. Will be passed to the command as the argument '-i'.
##         - `outputFilename`: Output filename for the alignment process
##         and input filename for the THD processing. Will be passed to
##         the command as the argument '-o'.
##         - `index`: Index of the inputFilename in the filtering output.
##         """
##         command = 'echo %s -o%s -i%s | qsub -cwd' % (command,
##                                                      outputFilename,
##                                                      inputFilename)
##         gs = GridSubmission(command)
##         waitForJob(gs.jobID)
##         if not self._alignmentDoneSuccesful(outputFilename):
##             raise Exception(('Alignment command didn\'t '
##                              'finish correctly.\n Command '
##                              'executed: %s' % (command,)))

##         # it's important to use the cloned class here to
##         # avoid connection clashes
##         # between threads
##         m8eMax=1
##         thd = THD.cloneAndReconnect().fromMbrFile(
##             outputFilename,
##             self._scoreCutoff,
##             self._sequenceFileID,
##             self._alignmentProtocolID,
##             self._taxonomyDBID,
##             filterProtocolID=self._filterProtocolID,
##             filterOutputIndex=index)
##         print >> sys.stdout, 'THD created. THD_ID: %d' % (thd,)

##     def run(self):
##         """Initiate all the process of creating the THDs
##         """

##         inputFilename = self.sequenceFile().Path
##         if not os.access(inputFilename,os.R_OK):
##             raise IOError, "can't open file: %s" % inputFilename
##         # noDups
##         self._filesPreffix = os.path.splitext(inputFilename)[0]
##         noDupFile = self.removeDuplicates(inputFilename)
##         readsC, sequencesC = countFasta(noDupFile)
##         readCountType = ReadCountType(Name='noDups',
##                                       insertIfNeeded=True)
##         ReadCount(Sequence_File_ID=self._sequenceFileID,
##                   Reads=readsC,
##                   Sequences=sequencesC,
##                   Read_Count_Type_ID=readCountType,
##                   insertIfNeeded=True)
##         print >> sys.stdout, ('noDups counts: reads: %i, '
##                               'sequences: %i' % (readsC, sequencesC))
        
##         noDupPrefix = os.path.splitext(noDupFile)[0]        
##         qualityFilteredFile = self.filterQuality(noDupFile)
##         readsC, sequencesC = countFasta(qualityFilteredFile)
##         readCountType = ReadCountType(Name='qualityFiltered',
##                                       insertIfNeeded=True)
##         ReadCount(Sequence_File_ID=self._sequenceFileID,
##                   Reads=readsC,
##                   Sequences=sequencesC,
##                   Read_Count_Type_ID=readCountType,
##                   insertIfNeeded=True)
##         print >> sys.stdout, ('qualityFiltered counts: reads: %i, '
##                               'sequences: %i' % (readsC, sequencesC))
        
##         self._filesPreffix = os.path.splitext(inputFilename)[0]
        
##         filterProto = BlastFilterProtocol(self._filterProtocolID)
##         stepsC = len(filterProto.stepList()) 
##         filteredFiles = self.applyFilterProtocol(qualityFilteredFile)
##         for fname in filteredFiles:
##             parts = fname.split('.')[-(stepsC+1):][:-1]
##             readCountName = ''.join(parts)
##             readCountType = ReadCountType(Name=readCountName,
##                                           insertIfNeeded=True)
##             readsC, sequencesC = countFasta(fname)
##             ReadCount(Sequence_File_ID=self._sequenceFileID,
##                       Reads=readsC,
##                       Sequences=sequencesC,
##                       Read_Count_Type_ID=readCountType,
##                       insertIfNeeded=True)
##             print >> sys.stdout, ('%s counts: reads: %i, '
##                                   'sequences: %i' % (readCountName,
##                                                      readsC,
##                                                      sequencesC))

##         self.startAlignmentJobs(filteredFiles)
##         self.waitForAlignmentJobs()


## class TaxHitDistributionSetSQL(TaxHitDistributionSet):
##     """Subclasses sequence.blastNoSQL.TaxHitDistributionSet to add
##     support to SQL serialization
##     """

##     @classmethod
##     def volatilSet(cls, THDlist):
##         """Returns a TaxHitDistributionSetSQL with
##         TaxHitDistributionSQL already loaded.
##         """
##         thds = cls()
##         for thd in THDlist:
##             thds[thd] = THD(thd)
##         return thds
    
##     def __init__(self, THDSID=None, *args, **kwargs):
##         """If `THDSID` is not None loads the THDS from the DB. Else pass
##         the rest of the arguments to the parent constructor.
        
##         Arguments:
##         - `THDSID`: THDS_ID from MegaChipDB.THDS table.
##         """
##         self._THDSID = THDSID
##         TaxHitDistributionSet.__init__(self, *args, **kwargs)
        
##         if THDSID is not None:
##             rowTHDS = THDS(int(THDSID))
##             for thdId in rowTHDS.thdList():
##                 thdId = thdId[0]
##                 self[thdId] = THD(thdId)


## class THDAnalyzer(object):
##     """Reads a configuration file or a THD_Analysis row from the DB,
##     do all the analysis following that configuration and store it in DB.
##     """
    
##     def __init__(self, thdAnalysisID=None, configFile=None, printResults=True):
##         """Initialize the instance. Requires one of thdAnalysis or
##         configFile to be specified.
        
##         Arguments:
##         - `thdAnalysisID`: THD_Analysis_ID from
##             MegaChipDB.THD_Analysis.
##         - `configFile`: INI file type with the configuration for this analysis.
##         """
##         if thdAnalysisID is not None:
##             self._thdAnalysisID = thdAnalysisID
##             confType = 'DB'
##         elif configFile is not None:
##             self._configFile = configFile
##             confType = 'File'
##         self._conf = dict(map(lambda x: (x, None),
##                                        ('pLimit eLimit expunged '
##                                         'topPercent THDID '
##                                         'IniRefTHDSID ExclRefTHDSID '
##                                         'FinalRefTHDSID obsNorm '
##                                         'refNorm refCC '
##                                         'overCalledOnly').split()))
##         self._loadConfiguration(confType)
##         self._printResults = printResults

##     def _loadConfiguration(self, configurationType):
##         """Loads all the configuration either from a INI file or from
##         the DB.
        
##         Arguments:
##         - `configurationType`: One of 'DB' or 'File'. This define from
##             where is going to be read the configuration.
##         """
##         if configurationType not in ('DB', 'File'):
##             raise Exception('Invalid configuration type.')
        
##         if configurationType == 'DB':
##             rowTHDA = THDAnalysis(int(self._thdAnalysisID))
##             self._conf['pLimit'] = rowTHDA.pLimit
##             self._conf['eLimit'] = rowTHDA.eLimit
##             taxons = rowTHDA.expungedTaxons()
##             if taxons is not None:
##                 self._conf['expunged'] = map(lambda x: x.Taxon_ID,
##                                              taxons)
##             self._conf['topPct'] = rowTHDA.Ref_Top_Tax_Pct
##             self._conf['THDID'] = rowTHDA.THD_ID
##             self._conf['IniRefTHDSID'] = \
##                                        rowTHDA.Initial_Reference_THDS_ID
##             self._conf['ExclRefTHDSID'] = rowTHDA.Ref_Excluded_THDS_ID
##             self._conf['FinalRefTHDSID'] = \
##                                          rowTHDA.Final_Reference_THDS_ID
##             self._conf['obsNorm'] = rowTHDA.Obs_Norm
##             self._conf['refNorm'] = rowTHDA.Ref_Norm
##             self._conf['overCalledOnly'] = rowTHDA.Over_Obs_Limit == 1
##             self._conf['refCC'] = rowTHDA.Ref_CC
##         else:
##             parser = ConfigParser.ConfigParser()
##             parser.read(self._configFile)
##             if parser.has_option('parameters', 'pLimit'):
##                 self._conf['pLimit'] = parser.getfloat('parameters',
##                                                         'pLimit')
##             if parser.has_option('parameters', 'eLimit'):
##                 self._conf['eLimit'] = parser.getfloat('parameters',
##                                                         'eLimit')

##             if parser.has_option('parameters', 'topPercent'):
##                 self._conf['topPct'] = \
##                                          parser.getfloat('parameters',
##                                                          'topPercent')
##             if parser.has_option('parameters', 'overCalledOnly'):
##                 self._conf['overCalledOnly'] = \
##                                              parser.getboolean('parameters',
##                                                                'overCalledOnly')
##             if parser.has_option('parameters', 'obsNorm'):
##                 self._conf['obsNorm'] = parser.getfloat('parameters',
##                                                         'obsNorm')
##             if parser.has_option('parameters', 'refNorm'):
##                 self._conf['refNorm'] = parser.getfloat('parameters',
##                                                         'refNorm')
##             if parser.has_option('parameters', 'refCC'):
##                 self._conf['refCC'] = parser.getfloat('parameters',
##                                                       'refCC')
##             if parser.has_option('parameters', 'sampleTHD'):
##                 self._conf['THDID'] = parser.getint('parameters',
##                                                     'sampleTHD')
##             else:
##                 raise Exception(('The analysis needs a THD to use as '
##                                  'sample.'))
##             if parser.has_option('parameters', 'refTHDS'):
##                 val = parser.get('parameters', 'refTHDS')
##                 try:
##                     val = int(val)
##                 except ValueError:
##                     if parser.has_section(val):
##                         desc = parser.get(val, 'description')
##                         thdids = parser.get(val, 'items').split(',')
##                         for thdid in thdids:
##                             try:
##                                 THD(int(thdid))
##                             except KdbomLookupError:
##                                 raise Exception(('THD %s don\'t '
##                                                  'exists.') % (thdid,))
##                         thdsObject = THDS.builder(desc, thdids)
##                         self._conf['IniRefTHDSID'] = \
##                                                      thdsObject.THDS_ID
##                     else:
##                         raise Exception(('Cannot find setion %s in the '
##                                          'configuration file %s.') % \
##                                         (val, self._configFile))
##                 else:
##                     self._conf['IniRefTHDSID'] = val
                    
##             if parser.has_option('parameters', 'exclTHDS'):
##                 val = parser.get('parameters', 'exclTHDS')
##                 try:
##                     val = int(val)
##                 except ValueError:
##                     if parser.has_section(val):
##                         desc = parser.get(val, 'description')
##                         thdids = parser.get(val, 'items').split(',')
##                         for thdid in thdids:
##                             try:
##                                 THD(int(thdid))
##                             except KdbomLookupError:
##                                 raise Exception(('THD %s don\'t '
##                                                  'exists.') % (thdid,))
##                         thdsObject = THDS.builder(desc,thdids)
##                         self._conf['ExclRefTHDSID'] = \
##                                                      thdsObject.THDS_ID
##                     else:
##                         raise Exception(('Cannot find setion %s in the '
##                                          'configuration file %s.') % \
##                                         (val, self._configFile))
##                 else:
##                     self._conf['ExclRefTHDSID'] = val

##             initID = self._conf['IniRefTHDSID']
##             initRef = TaxHitDistributionSetSQL(THDSID=initID)

##             exclID = self._conf['ExclRefTHDSID']
##             if exclID is not None:
##                 iniSet = set(initRef.keys())
##                 thdsExcl = TaxHitDistributionSetSQL(THDSID=exclID)
##                 excSet = set(thdsExcl.keys())
##                 finSet = iniSet.difference(excSet)
##                 finalRef = TaxHitDistributionSetSQL.volatilSet(finSet)
##             else:
##                 finalRef = initRef                

##             frDesc_dissRef = ''
##             if self._conf['refCC'] is not None and \
##                    self._conf['topPct'] is not None:
##                 finalRef = finalRef.dissimilarRef([],
##                                                   corrMax=self._conf['refCC'],
##                                                   topPctMax=self._conf['topPct'])
##                 frDesc_dissRef = 'refCC-%.2f-topPct-%.2f' % \
##                                  (self._conf['refCC'], self._conf['topPct'])

##             iniDesc = THDS(self._conf['IniRefTHDSID']).Description
##             excDesc = THDS(self._conf['ExclRefTHDSID']).Description
##             frDesc = '%s using %s. Params: %s' % (iniDesc, excDesc,
##                                                   frDesc_dissRef)
##             finalTHDS = THDS.builder(frDesc, finalRef.keys())

##             self._conf['FinalRefTHDSID'] = finalTHDS.THDS_ID

##             thda = THDAnalysis(THD_ID=self._conf['THDID'],
##                                Initial_Reference_THDS_ID=self._conf['IniRefTHDSID'],
##                                Ref_Excluded_THDS_ID=self._conf['ExclRefTHDSID'],
##                                Final_Reference_THDS_ID=self._conf['FinalRefTHDSID'],
##                                Ref_Top_Tax_Pct=self._conf['topPercent'],
##                                Ref_CC=self._conf['refCC'],
##                                Obs_Norm=self._conf['obsNorm'],
##                                Ref_Norm=self._conf['refNorm'],
##                                pLimit=self._conf['pLimit'],
##                                eLimit=self._conf['eLimit'],
##                                Over_Obs_Limit=self._conf['overCalledOnly'],
##                                Date=datetime.datetime.now(),
##                                insertIfNeeded=True)
##             self._thdAnalysisID = thda.THD_Analysis_ID
##             if parser.has_option('parameters', 'expunged'):
##                 expunged = parser.get('parameters', 'expunged')
##                 txs = map(int, expunged.split(','))
##                 thdaid = thda.THD_Analysis_ID
##                 THDAExpunged._table.insertMany(('THD_Analysis_ID',
##                                                 'Taxon_ID'),
##                                                ((thdaid, txid)
##                                                 for txid in txs))
##                 self._conf['expunged'] = txs
##             else:
##                 self._conf['expunged'] = []

##     def run(self):
##         """Performs the analysis using all the parameters previously
##         loaded to the instance of the class. Will store the results in
##         the table THDA_Results.
##         """
##         thdssClass = TaxHitDistributionSetSQL
##         sampleSet = thdssClass.volatilSet([self._conf['THDID'],])
##         refSet = thdssClass(THDSID=self._conf['FinalRefTHDSID'])

##         if len(self._conf['expunged']) > 0:
##             for expTaxon in self._conf['expunged']:
##                 try:
##                     sampleSet.expungeTaxon(expTaxon)
##                 except KeyError:
##                     pass
##                 try:
##                     refSet.expungeTaxon(expTaxon)
##                 except KeyError:
##                     pass

##         if self._conf['obsNorm'] is not None:
##             sampleSet.normalize(self._conf['obsNorm'])
##         if self._conf['refNorm'] is not None:
##             refSet.normalize(self._conf['refNorm'])

##         results = sampleSet.speciesReport(reference=refSet,
##                                           pLimit=self._conf['pLimit'],
##                                           eLimit=self._conf['eLimit'],
##                                           overCalledOnly=True,
##                                           returnValues=True)
##         thdarTable = THDAResults.table()
##         thdarTable.insertMany(('THD_Analysis_ID', 'Taxon_ID', 'P',
##                                'Obs_Pct', 'Ref_Pct'),
##                               ((self._thdAnalysisID, tax, p, o, e)
##                                for tax, p, o, rw in results)
##                               )
##         thdarTable.db.commit()

##         if self._printResults:
##             r = sampleSet.speciesReport(reference=refSet,
##                                         pLimit=self._conf['pLimit'],
##                                         eLimit=self._conf['eLimit'],
##                                         overCalledOnly=True,
##                                         returnValues=False)
##             print >> sys.stdout, r
            
