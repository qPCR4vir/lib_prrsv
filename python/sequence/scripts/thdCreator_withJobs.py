#!/usr/local/bin/python
#
# Creates THD from sequening reads.
#
__version__ = tuple([int(x) for x in
                     '$Revision: 1.4 $'.split()[1].split('.')])
__author__ = "Julio Menendez, Kael Fischer"

import sys
from optparse import OptionParser

from viroinfo.megaseq import THD, THDCreator
from grid.jobs import Job, JobDescription, JobManager


SCRIPTS_FOLDER = '/r1/home/julio/lib/python/sequence/scripts/thdCreatorScripts'


def main(args):
    parser = OptionParser('usage: %prog [options]')
    optionsAvailable = (
        ('-i', '--input', 'Sequence_File_ID', 'int',
         'Sequence_File_ID from MegaChip.Sequence_File table.'),
        ('-a', '--alignment', 'Alignment_Protocol_ID', 'int',
         ('Alignment_Protocol_ID from MegaChip.Alignment_Protocol '
          'table.')),
        ('-t', '--taxonomy', 'Taxonomy_ID', 'int',
         'Taxonomy_ID from MegaChipDB.Taxonomy table.'),
        ('-e','--expect','Alignment_Score_Limit','float',
         'Score cutoff for alignment.'),
        ('-f', '--filter', 'Blast_Filter_Protocol_ID', 'int',
         ('Blast_Filter_Protocol_ID from MegaChipDB.'
          'Blast_Filter_Protocol table. If None no filter will be '
          'used.')),
        ('-l', '--likeTHD', 'templateTHD_ID','int',
         'use exsiting THD as template for unspecified values'),
        )

    ids =[ x[2] for x in optionsAvailable if x[2] != 'templateTHD_ID' ]

    for opt in optionsAvailable:
        if opt[3] is None:
            otype = 'string'
        else:
            otype = opt[3]
        parser.add_option(opt[0], opt[1], dest=opt[2], type=otype,
                          help=opt[4])
    (options, args) = parser.parse_args(args)


    defined = lambda x: getattr(options, x) is not None

    if options.templateTHD_ID != None:
        templateTHD=THD(options.templateTHD_ID)
        for opt in ids: 
            if not defined(opt):
                setattr(options,opt,getattr(templateTHD,opt))
        
    if not all(map(defined, ids)):
        print >> sys.stderr, ('Usage error: all arguments are required, '
                              'or derived from template THD')
        print >> sys.stderr, parser.format_help()
        sys.exit(1)

    thdc = THDCreator(options.Sequence_File_ID,
                      options.Alignment_Protocol_ID,
                      options.Alignment_Score_Limit,
                      options.Taxonomy_ID,
                      options.Blast_Filter_Protocol_ID)

    outputFiles = thdc.allOutputFilenames()
    print outputFiles
    
    script = ('python %s/removeDups.py -i%%(INinitial)s '
              '-o%%(OUTnoDups)s -s%%(seqFileId)s' % SCRIPTS_FOLDER)
    remDupsDesc = JobDescription(Command=script, insertIfNeeded=True)
    remDupsJob = Job(Job_Description_ID=remDupsDesc,
                     INinitial=str(thdc.sequenceFile()),
                     OUTnoDups=outputFiles[0],
                     seqFileId=thdc._sequenceFileID)
    script = ('python %s/filterQuality.py -i%%(INnoDups)s '
              '-o%%(OUTqualFiltered)s -s%%(seqFileId)s' %
              SCRIPTS_FOLDER)
    filQualDesc = JobDescription(Command=script, insertIfNeeded=True)
    filQualJob = Job(Job_Description_ID=filQualDesc,
                     INnoDups=remDupsJob.outputs[0],
                     OUTqualFiltered=outputFiles[1],
                     seqFileId=thdc._sequenceFileID)

    outputVars = {}
    tmpOutVars = []
    for idx, outFn in enumerate(outputFiles[2]):
        n = 'OUTfile%d' % idx
        outputVars[n] = outFn
        tmpOutVars.append('%%(%s)s' % n)
    script = ('python %s/applyBlastFilters.py -i%%(INqualFiltered)s '
              '-s%%(seqFileId)s -f%%(filterId)s %s' %
              (SCRIPTS_FOLDER, ' '.join(tmpOutVars)))
    blastFilterDesc = JobDescription(Command=script,
                                     insertIfNeeded=True)
    blastFilterJob = Job(Job_Description_ID=blastFilterDesc,
                         INqualFiltered=filQualJob.outputs[0],
                         seqFileId=thdc._sequenceFileID,
                         filterId=thdc._filterProtocolID,
                         **outputVars)

    script = ('python %s/doAlignment.py -i%%(INfilteredFile)s '
              '-o%%(OUTmbrFile)s -a%%(alignProtocolId)s' %
              SCRIPTS_FOLDER)

    mbrDesc = JobDescription(Command=script, insertIfNeeded=True)

    script = ('python %s/createThd.py -i%%(INmbrFile)s '
              '-s%%(seqFileId)s -x%%(filterindex)s '
              '-a%%(alignProtocolId)s -t%%(taxonomyId)s '
              '-e%%(scoreLimit)s -f%%(filterId)s' % SCRIPTS_FOLDER)
    thdDesc = JobDescription(Command=script, insertIfNeeded=True)


    addJobs = []
    for idx, filteredFile in enumerate(outputFiles[2]):
        mbrJob = Job(Job_Description_ID=mbrDesc,
                     INfilteredFile=filteredFile,
                     OUTmbrFile=outputFiles[3][idx],
                     alignProtocolId=thdc._alignmentProtocolID)
        
        thdJob = Job(Job_Description_ID=thdDesc,
                     INmbrFile=mbrJob.outputs[0],
                     seqFileId=thdc._sequenceFileID, filterindex=idx,
                     alignProtocolId=thdc._alignmentProtocolID,
                     taxonomyId=thdc._taxonomyDBID,
                     scoreLimit=thdc._scoreCutoff,
                     filterId=thdc._filterProtocolID)
        
        addJobs.append(mbrJob)
        addJobs.append(thdJob)

    jm = JobManager([remDupsJob, filQualJob, blastFilterJob] + \
                    addJobs)
    jm.buildTree()
    jm.run()

    print 'Root JobManager: %d' % int(jm)
        
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
