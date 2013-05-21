"""BLAST utilities
BLAST result storage and retrieval.
"""

import commands
import copy
import copy_reg
import os
import os.path
import datetime
from types import *

import utils

from __init__ import *
from fasta import *
from blastNoSQL import *

from kdbom2 import kdbom


blastDB = kdbom.tryDBconnect(db='Blast_Results',
                             serversPorts='kf-db1',
                             getTables=False,
                             fatal=True)


SQL = {'createHspTable':
       """
       CREATE TABLE IF NOT EXISTS %s (
       `Hsp_ID` bigint unsigned NOT NULL AUTO_INCREMENT,
       `Q_Gi` varchar(255) NOT NULL,
       `T_Gi` varchar(255) NOT NULL,
       `Ident` float unsigned NOT NULL,
       `Length` int unsigned NOT NULL,
       `Gaps` int unsigned NOT NULL,
       `Mismatch` int unsigned NOT NULL,
       `Q_Start` int unsigned NOT NULL,
       `Q_End` int unsigned NOT NULL,
       `T_Start` int unsigned NOT NULL,
       `T_End` int unsigned NOT NULL,
       `E` double unsigned NOT NULL,
       `Score` float NOT NULL,
       KEY (`Hsp_ID`),
       KEY `qgi` (`Q_Gi`,`T_Gi`),
       KEY `tgi` (`T_Gi`,`Q_Gi`)
       ) ENGINE=MyISAM DEFAULT CHARSET=latin1
       """
    }


MAXRECORDS = 10000000



class Search (kdbom.KSqlObject):
    """
    """
    _table = blastDB.Search
    _strField=_table.Hsp_Table_Suffix

    def __post_init__(self):
        if not self.hspTableExists():
            self.createHspTable()

    def hspTableName(self):
        """returns the name of the associated Hsp table.
        """
        return 'Hsp' + self.Hsp_Table_Suffix


    def hspDbQualName(self):
        """returns the quoted, datatbase qualified name of the Hsp Table
        """
        return '`%s`.`%s`' % (self.db().name, self.hspTableName())
 

    def hspTableExists(self):
        """True if Hsp table exists, otherwise false
        """
        return len(
            self.db().fetchall('SHOW TABLES LIKE %s',self.hspTableName()))==1
        

    def hspClass(self):
        """Return the KSqlObject class of the Hsps from this search.
        """
        return Hsp.makeHspClass(self.Hsp_Table_Suffix)


    def createHspTable(self):
        """Create HSP mysql table.
        """
        self.db().execute(SQL['createHspTable']%self.hspDbQualName())
        self.db().refresh()
        return self.db()[self.hspTableName()]

    def hspClass (self):
        """Return Hsp KSqlObject class
        """
        return Hsp.makeHspClass(self.Hsp_Table_Suffix)


    

class Hsp (kdbom.KSqlObject):
    """
    """
    _tableSuffix=''
    _table = blastDB.Hsp
    

    @classmethod
    def insertM8records(cls,m8paths):
        if type(m8paths) in StringTypes:
            m8paths = [m8paths]

        curRowCt = cls._table.count()
        insertEst = estimateM8rows(m8paths)
        cls.repartitionTable(curRowCt+insertEst)

        cls._table.insertMany(
            ('Q_Gi','T_Gi',
             'Ident','Length',
             'Gaps','Mismatch',
             'Q_Start','Q_End',
             'T_Start','T_End',
             'E','Score'),
            m8tupleGenerator(m8paths),chunkSize=1000000,
            partialCommits=True)

    @classmethod
    def repartitionTable(cls,rowCt):
        """
        """
        partCt = rowCt/MAXRECORDS + 1
        
        cls._table.db.execute(
            """ALTER TABLE %s
            PARTITION BY KEY(`Q_Gi`,`T_Gi`)
            PARTITIONS %s
            """ % (cls._table.dbQualName(),partCt))

    @classmethod
    def makeHspClass(cls,tableSuffix):
        """Return an Hsp class corresponding to a subset of Hsps,
        held in particular table.
        """
        import copy
        tblName= cls._table.name+tableSuffix
        try:
            t = cls.db()[tblName]
        except kdbom.MySQLdb.ProgrammingError, err:
            if err[0] == 1146:
                cls.db().execute(SQL['createHspTable']%tblName)
                cls.db().refresh(getTables=False)
                t=cls.db()[tblName]
            else:
                raise
        
        class HspPartition (cls):
            _table=t

        
        
        return HspPartition


    def m8tuple(self):
        """return a tuple of the m8fields
        """
        return (self.Q_Gi,self.T_Gi,
                self.Indent,self.Length,
                self.Gaps,self.Mismatch,
                self.Q_Start, self.Q_End,
                self.T_Start, self.T_End,
                self.E, self.Score)


    def __repr__(self):
        return '\t'.join(self.ID(),
                         self.Q_Gi,self.T_Gi,
                         self.Indent,self.Length,
                         self.Gaps,self.Mismatch,
                         self.Q_Start, self.Q_End,
                         self.T_Start, self.T_End,
                         self.E, self.Score)


    def __str__(self):
        return m8formatTuple(self)

    
    def m8dict(self):
        """Returns a dictionary representation of the hit
        like those returned by m8generator.
        """
        (qryName,subjName,ident,length,gaps,mismatch,
         q_start,q_end,t_start,t_end,e,score) = self.m8tuple()
        
        return{
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
            'score':float(score)
            }    

    @classmethod
    def pickler(cls):
        return self._tableSuffix


def makeHspClass(tableSuffix):
    """Return an Hsp subclass corresponding to a subset of Hsps,
    held in particular table.
    """

    tblName= Hsp._table.name+tableSuffix
    try:
        t = Hsp.db()[tblName]
    except kdbom.MySQLdb.ProgrammingError, err:
        if err[0] == 1146:
            Hsp.db().execute(SQL['createHspTable']%tblName)
            Hsp.db().refresh(getTables=False)
            t=Hsp.db()[tblName]
        else:
            raise

    class HspPartition (Hsp):
        _table=t
        _tableSuffix=tableSuffix

    return HspPartition

def insertBlastResults(dbFilePath,
                       qryFilePath,
                       program,
                       blastRunParameters,
                       m8resultPaths,
                       hspTableSuffix):
    """Insert results from a blast run in to the Blast_Results
    database.
    """

    def fileMtime(path):
        mts = os.stat(path).st_mtime
        return datetime.datetime.fromtimestamp(mts)


    batchSize = 20000

    today = datetime.datetime.today()

    if type(m8resultPaths) in StringTypes:
        m8resultPaths = [m8resultPaths]

    # check paths
    dataPaths = map(os.path.abspath,[ dbFilePath,qryFilePath] + m8resultPaths)
    for fpath in dataPaths:
        if not os.access(fpath,os.R_OK):
            raise IOError, "%s is unreadable" % fpath

    #
    # Blast Search Record
    #
    search = Search({'Program': program,
                     'Parameters':blastRunParameters,
                     'Query_Path':qryFilePath,
                     'Db_Path':dbFilePath},
                    insertIfNeeded=True)
    searchID=search.ID()
    hsps = search.hspClass()

    #
    # For each result file from this search
    # insert hits and hsps
    for m8resultPath in m8resultPaths:
        
        hsps._table.insertMany(
            ('Q_Gi','T_Gi',
             'Ident','Length',
             'Gaps','ismatch',
             'Q_Start','Q_End',
             'T_Start','_End',
             'E','Score'),
            m8tupleGenerator(mresultPath),chunkSize=50000,
            callback=m8demangleTiles,
            partialCommits=True)
