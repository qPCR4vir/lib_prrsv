from kdbom2 import kdbom
from sequence import fasta as F
from sequence import blastNoSQL as B
from sequence import aos as A
import numpy as N

SQLservers = (('kf-db1',3306),)
prrsvDB = kdbom.tryDBconnect('PRRSV',SQLservers)

class Sequence(kdbom.KSqlObject):
    _table = prrsvDB.Sequence
    
    def generateFasta(self):
        '''generates a fasta file using Sequence_ID as
        the title and Sequence as the sequence
        '''
        rec=F.Record(title=self.Sequence_ID, sequence=self.Sequence)
        return rec

    def taxa(self):
        """return generator of taxa that I belong to
        """
        return Taxon.generator(Sequence.Sequence_ID==self)

    fasta=generateFasta
        

class Candidate(kdbom.KSqlObject):
    _table = prrsvDB.Candidate

    energyStats=dict(
        zip(("HUN4",
             "VR2332",
             "MN184",
             "JA142",
             "TYPE1",
             "3UTR",
             "5UTR",
             "ORF7",
             "GP5",
             "NSP2",),
            ((-63.11740260763215,14.43184756330354),
             (-60.988265420785304,14.272450369924032),
             (-57.83889727434975,14.060091651474204),
             (-62.40689888325441,12.334703817668343),
             (-67.51765860871174,16.190421772570268),
             (-67.1952420477331,12.8864831008811),
             (-74.75313356348771,11.55206525697237),
             (-65.50831532127117,9.058938049132694),
             (-60.275681619582265,9.656863251585118),
             (-76.88749323924402,12.485088710519856),)))

    def energyAnnotations(self,zCutoff=-1):
        """
        """
        possible=('HUN4', 'VR2332', 'MN184',
                  'JA142', 'TYPE1', '3UTR', '5UTR',
                  'ORF7', 'GP5', 'NSP2')
        q=self._table.db.Energy_Crosstab.query()
        q.select([t+'_Mean' for t in possible])

        q.where(q._table.Candidate_ID==self.ID())

        e= q.fetchall()[0]

        vTypeE=sorted(zip(e[:5],possible[:5]))
        vTypeE=[x for x in vTypeE if x[0] is not None]
        vTypes=[x[1] for x in vTypeE if x[0] <=vTypeE[0][0]]
        
        vRegionE=sorted(zip(e[-5:],possible[-5:]))
        vRegionE=[x for x in vRegionE if x[0] is not None]
        vRegions=[x[1] for x in vRegionE if x[0] <=vRegionE[0][0]]

        return vTypes+vRegions

    

    @classmethod
    def insertManyFromFile(cls,iFile):
        '''
        '''
        def rows():
            '''yield the value of name and value of sequence using
            FastaIterator
            '''
            iFasta=F.FastaIterator(iFile)
            for rec in iFasta:
                yield [rec.title,rec.sequence]

        
        cls._table.insertMany(('Name','Sequence'),rows())

    def sequence(self):
        return Sequence(int(self.Name.split('_')[0]))

    def fasta(self):
        return F.Record(sequence=self.Sequence,title=self.Name)

    def taxaMinE(self,cond=None):
        """return a list of:
        (Candidate_Energy.Taxon_ID,Candidate_Energy.Min_E)

        any cond should try to restrict itsef to the
        Candidate_Energy table or you may have problems.
        """
        q=CandidateEnergy.query()
        q.select([CandidateEnergy.Taxon_ID,CandidateEnergy.Min_E])
        q.where((CandidateEnergy.Candidate_ID == self.ID()) & cond)
        return q.fetchall()





class CandidateEnergy(kdbom.KSqlObject):
    _table = prrsvDB.Candidate_Energy
    @classmethod
    def insertManyFromBlastFile(cls,iFile):
        '''
        '''
        def rows():
            '''
            '''
            canCache={}
            pairCache={}
            
            lastCanName=''
            hg=B.ncbiXmlHitGenerator(iFile)
            i=0
            for rec in hg:
                i+=1
                seqID=int(rec['def'])
                #1st variable found
                
                canName=rec['query-def']
                #2nd variable found

                if canName!=lastCanName:
                    for tID,canID in pairCache.iterkeys():
                        min_E=min(pairCache[tID,canID])
                        max_E=max(pairCache[tID,canID])
                        mean_E=N.mean(pairCache[tID,canID])
                        SD_E=N.std(pairCache[tID,canID])

                        N_E=0
                        N_E30=0
                        N_E40=0
                        N_E50=0
                        N_E100=0
                        for ener in pairCache[tID,canID]:
                            N_E+=1
                            if ener<=-30:
                                N_E30+=1
                            if ener<=-40:
                                N_E40+=1
                            if ener<=-50:
                                N_E50+=1
                            if ener<=-100:
                                N_E100+=1
                        rv=[canID,tID,min_E,max_E,mean_E,SD_E,N_E,N_E30,N_E40,
                            N_E50,N_E100]
                        yield rv
                    pairCache={}
                    lastCanName=canName
                #stores and refreshes dictionary
                #once we have moved onto next pair
                    
                  
                try:
                    canID=canCache[canName]
                except KeyError:
                    canCache[canName]=Candidate(canName).ID()
                    canID=canCache[canName]
                #timesaver so we only have to look it up once

                if canID%4!=0:
                    continue
                    
                #hit=str(rec)
                #3rd variable found
                
                qseq=rec['hsps'][0]['qseq']
                hseq=rec['hsps'][0]['hseq']
                rqseq=B.reverseComplement(qseq)
                en=A.energy(hseq,rqseq)
                #3rd variable found
                
                for t in Sequence(seqID).taxa():
                    tID=t.ID()
                    try:
                        pairCache[tID,canID].append(en)
                    except KeyError:
                        pairCache[tID,canID]=[en]
                #checks to see if energy is most negative for the pair
                #replaces old result if it is



            for tID,canID in pairCache.iterkeys():
                        min_E=min(pairCache[tID,canID])
                        max_E=max(pairCache[tID,canID])
                        mean_E=N.mean(pairCache[tID,canID])
                        SD_E=N.std(pairCache[tID,canID])

                        N_E=0
                        N_E30=0
                        N_E40=0
                        N_E50=0
                        N_E100=0
                        for ener in pairCache[tID,canID]:
                            N_E+=1
                            if ener<=-30:
                                N_E30+=1
                            if ener<=-40:
                                N_E40+=1
                            if ener<=-50:
                                N_E50+=1
                            if ener<=-100:
                                N_E100+=1
                        rv=[canID,tID,min_E,max_E,mean_E,SD_E,N_E,N_E30,N_E40,
                            N_E50,N_E100]
                        yield rv
            #empty dictionary one last time
        
        cls._table.insertMany(('Candidate_ID','Taxon_ID','Min_E','Max_E',
                               'Mean_E','SD_E','N_E','N_E30','N_E40','N_E50',
                               'N_E100')
                              ,rows(),ignore=True)
        

class Taxon(kdbom.KSqlObject):
    _table = prrsvDB.Taxon



class TaxonHasSequence(kdbom.KSqlObject):
    _table = prrsvDB.Taxon_has_Sequence



class CandidateSequenceEnergy(kdbom.KSqlObject):
    _table = prrsvDB.Candidate_Sequence_Energy
    _table.addRelationship(_table.Sequence_ID,
                           TaxonHasSequence._table.Sequence_ID)

    
    @classmethod
    def insertManyFromBlastFile(cls,iFile):
        '''
        '''
        def rows():
            '''takes xml blast file and returns a generator of best
            pairwise results
            '''
            hg=B.ncbiXmlHitGenerator(iFile)
            pairCache={}
            canCache={}
            i=0
            for rec in hg:
                i+=1
                seqID=int(rec['def'])
                #1st variable foundandidate_Sequence_Energy
                
                canName=rec['query-def']
                #2nd variable found

                qseq=rec['hsps'][0]['qseq']
                hseq=rec['hsps'][0]['hseq']
                rqseq=B.reverseComplement(qseq)
                en=A.energy(hseq,rqseq)
                #3rd variable found

                try:
                    canID=canCache[canName]
                except KeyError:
                    canCache[canName]=Candidate(canName).ID()
                    canID=canCache[canName]
                #timesaver so we only have to look it up once


                try:
                    if en < pairCache[canID,seqID]:
                        pairCache[canID,seqID]=en
                except KeyError:
                    pairCache[canID,seqID]=en
                #store the most negative value for the pair


            for canID,seqID in pairCache.iterkeys():
                rv=[canID,seqID,pairCache[canID,seqID]]
                yield rv

        cls._table.insertMany(('Candidate_ID','Sequence_ID','Energy'),
                              rows(),ignore=True)

chipTitles = ['03_VR2332_2',
              '06_MN184_3',
              '07_MN184_1',
              '09_SDSU73_0',
              '100_Leylstad_7',
              '102_MN184_5',
              '106_MN184_0',
              '107_SDSU73_7',
              '108_SDSU73_3',
              '109_MN184_2',
              '110_SDSU73_4',
              '112_Leylstad_6',
              '113_SDSU73_6',
              '115_VR2332_5',
              '118_VR2332_6',
              '13_VR2332_1',
              '19_VR2332_0',
              '23_MN184_4',
              '45_Leylstad_3',
              '53_SUSD73_1',
              '73_Leylstad_5',
              '74_Leylstad_1',
              '75_Leylstad_0',
              '76_Leylstad_2',
              '78_Leylstad_4',
              '92_SRD80-nasal_0',
              '93_Leylstad_8',
              '95_Leylstad_9',
              '98_SDSU73_2',
              '99_SRD80-serum_0']

chipMap=dict(((x.split('_')[0],x) for x in chipTitles))

def gprFN2chipNumber(fn):
    return fn.split('_')[1]
    
