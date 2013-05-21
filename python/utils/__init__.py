#!/usr/local/bin/python
#
# utilities
#
# DeRisi Lab 2007-2008
# Fischer Lab 2008-
#
#
# Note: itertools has replaced some of this
#
__version__  =[int(x) for x in "$Revision: 1.49 $".split()[1].split('.')]


import copy
import cPickle
import sys, os, re
from types import *
import operator
import tempfile
from StringIO import StringIO
import math

from tabbedfile import tabbifyMatrix,txtToMatrix
from iterators import *
from decorators import *
from CommandLine import *
from timeprofile import *


__all__=["CommandLine","timeprofile",'readArgsOrFiles','flatten','unique',
         'addDict','subtractDict','normDict','union','intersection','shortest',
         'tabbifyMatrix','txtToMatrix','batchList','mystemp','iterCount','dict2optStr',
         'time2datetime','noneDict','compress','decompress','safeOFW',
         'MultiDictCluster','nanSum','multiFile','PickleableObject',
         'BitFlagConstants','dicts2array','vec2col','vec2row','fileIterator',
         'getIterable','unzip','lastLine','benchmark','wrap','SafeDict']


class BitFlagConstants(dict):
    """
    """
    def __init__(self,attrNames,startingBit=0,unsetName=None):
        """
        
        Arguments:
        - `self`:
        - `*attrNames`:
        """

        self[unsetName]=0
        self.mask=0
        for i,n in enumerate(attrNames):
            self[n]=(1 << (i+startingBit))
            self.mask = self.mask | (1 << (i+startingBit))


        #setup reverse lookup
        for k,v in self.items():
            self[v]=k

            
    def __getattr__ (self,thing):
        """
        """
        try:
            return self[thing]
        except KeyError:
            raise AttributeError

    def masked(self,n):
        """returns the bits of n that the object knows about. 
        """
        return n & self.mask

    def anySet(self,n):
        """True if any flags are set in n (or  False)
        """
        return bool(self.masked(n))

    def numberSet(self,n):
        """returns the count of the number of set bits in n
        """
        pass

def dicts2array(*dicts):
    """return a numpy array where the first column is the union of the
    keys in all dicts and the values of each following column are the
    values in the given dictionaries, in turn. 
    """
    import numpy as N
    keySets = [set(d.keys()) for d in dicts]
    keys=N.array(list(reduce(lambda x,y: x.union(y), keySets)))
    data=[N.array(list(safeDictLookup(d,0,*keys))) for d in dicts]
    #print data

    keys=vec2col(keys)
    return reduce(lambda x,y: N.c_[x,vec2col(y)],[keys]+data)
    


    
def safeDictLookup(d,default,*keys):
    """return a generator of values from dict over keys,
    if a key is not in d, return default (None by default).
    """
    for k in keys:
        try:
            yield d[k]
        except KeyError:
            yield default
        
        
def nanSum(*things):
    """Return what you get when you add all the things together.  The things don't
    have to be numbers; but must have __add__ methods.
    """
    return reduce(operator.add, *things)


def union(*lists):
    """return the union of input lists of records"""
    return noneDict(nanSum(*lists)).keys()


def unique(iteratable,remove=[]):
    '''Given an iteratable of hashables return a list of
    unique members

    things in remove will be removed from the list.
    '''

    l = list(iteratable)
    rd=dict(zip(l,[None]*len(l)))
    for thing in remove:
        if thing in rd:
            del rd[thing]
    return rd.keys()


def shortest(input):
    """Return shortest list/string/tuple in input list"""
    weeLength = False
    shortest = None

    for thing in input:
        if not weeLength:
            shortest = thing
            weeLength = len(thing)
        elif len(thing) < weeLength:
            shortest = thing
            weeLength = len(thing)
    return shortest


def intersection(*lists):
    """return the intersection of input lists of records"""

    # this should be rewritten using dictionaries
    rv = []
    rvPKs = []
    notCommonPKs = []
    strippedLists = []

    # if this is a single list remove
    # duplicate records
    if len(lists) == 1:
        return union(lists[0])

    #strip duplicates
    shortest = -1
    shortest_i = -1
    for i in range(len(lists)):
        if type(lists[i]) == LIST:
            strippedLists.append(union(lists[i]))
            if len(lists[i]) < shortest or shortest == -1:
                shortest = len(lists[i])
                shortest_i = i
    
    for elem in strippedLists[shortest_i]:
        # for every element in the shortest list
        # that is not already in the return vector
        # or found to be absent in a list
        if elem.PKvalue not in rvPKs and elem.PKvalue not in notCommonPKs:
            for bList in strippedLists:
                # for every elem in each following list
                # check equality
                for elem2 in bList:
                    foundin_bList = NO
                    if elem == elem2:
                        foundin_bList = YES
                        break
                # ok thats all the elems in the next list
                # or we found it
                if foundin_bList == NO:
                    notCommonPKs.append(elem.PKvalue)
                    break

            if elem.PKvalue not in notCommonPKs:
                rvPKs.append(elem.PKvalue)
                rv.append(elem)
    return rv


def safeOFW(fn,clobber=False,append=False):
    """Open a file for writing or appending and avoid clobbering.
    If a file named (fn) already exists, ValueError is raised unless
    either clobber or append is True.  They have the usual meanings.
    If both are true, ValueError is also raised.

    If fn = '-', sys.stdout is returned.
    If a writeable open file is passed, it is returned unmodified.
    """

    if fn=='-':
        return sys.stdout
    elif hasattr(fn,'write'):
        return fn
    else:
        if clobber and append:
            raise ValueError, ("At most one of 'clobber' and "
                               "'append' may be True.")
        fExists=os.access(fn,os.F_OK)
        if not fExists or (fExists and clobber):
            return file(fn,'wb')
        if fExists and append:
            return file(fn,'ab')
        raise ValueError, "file exists: %s" % fn
    

def time2datetime(timeObj):
    """convert time (from time module)
    to datetime from datetime module
    """
    import datetime
    return datetime.datetime(*timeObj[0:6])


def dict2optStr(oDict):
    """transform a dictionary to a command line options string.
    e.g. {'-o': 'short', '-x': None} => "-o short -x"
    """
    return ' '.join(['%s %s' % (opt,arg) for opt,arg in oDict.items() ])

def flatten(l,ltypes=(ListType,TupleType)):       
    """flatten a n-dimentional sequence.

    by Mike C. Fletcher
    """
    ltype = type(l)
    l=list(l)
    i=0
    while i<len(l):
        while isinstance(l[i],ltypes):
            if not l[i]:
                l.pop(i)
                i-=1
                break
            else:
                l[i:i+1]=l[i]
        i+=1
    return ltype(l)
            
def addDict(d1,d2):
    """Given 2 dictionaries with values that support addition
    numbers, strings, lists, tuples returns a dictionary
    with shared keys' vaules which are the sum of the like-keyed values.
    Unique keys in d1, d2 are also in the result.
    """

    if len(d1) < len(d2):
        dTemp=d1
        d1=d2
        d2=dTemp

    # copy largest dict
    rd = copy.deepcopy(d1)
    
    #iterate over smaller dict
    for k,v in d2.items():
        if k not in rd:
            rd[k]=v
        else:
            rd[k]+=v
    return rd


def subtractDict(d1,d2):
    """Given 2 dictionaries with values that support subtraction
    returns a dictionary with shared keys' vaules which are the difference
    of the like-keyed values. Unique keys in d1, d2 are also in the result.
    """

    if len(d1) < len(d2):
        dTemp=d1
        d1=d2
        d2=dTemp
        invert = True
    else:
        invert = False

    # copy largest dict
    rd = copy.deepcopy(d1)
    
    #iterate over smaller dict
    for k,v in d2.items():
        if k not in rd:
            if invert:
                rd[k]=v
            else:
                rd[k]=-v
        else:
            if invert:
                rd[k]=v-rd[k]
            else:
                rd[k]-=v
    return rd


def normDict(d,s=1):
    """Given a dictionary with values that support sum, and multiplication
    returns a dictionary the values normalized such that their sum is s.
    """
    
    s=float(s)
    t=float(sum(d.values()))
    if t != 0.0:
        nf=s/t
        return dict((x,y*nf) for x,y in d.items())
    else:
        return d

@deprecated
def noneDict(sequence):
    """UPDATE: use set() instead (from stdlib).
    For a sequence of hashable things, return
    a dictionary keyed on the things, with the values = None.
    Any iterable is allowed.

    Useful when you need a hashed set for faster set like opperations,
    e.g. unique, intersection, etc.

    """
    sequence=tuple(sequence)
    return dict(zip(sequence,[None]*len(sequence)))


def unzip(items):
    """Performs the reverse of zip, i.e. produces separate lists from tuples.
  
    items should be list of k-element tuples. Will raise exception if any tuples
    contain more items than the first one.
  
    Conceptually equivalent to transposing the matrix of tuples.
  
    Returns list of lists in which the ith list contains the ith element of each
    tuple.
  
    Note: zip expects *items rather than items, such that unzip(zip(*items))
    returns something that compares equal to items.
  
    Always returns lists: does not check original data type, but will accept
    any sequence.
    """
    if items:
        return map(list, zip(*items))
    else:
        return []
  



def readArgsOrFiles (args,
                     CL_INPUT_ALLOWED = True,
                     STDIN_ALLOWED = True,
                     PATHS_ALLOWED = True,
                     CGI_OVERRIDE=False,
                     UNIVERSAL_NL=True,
                     cgiParam='data',
                     extraCGIDelimChars = ''):
    """Takes a list of strings, normally command line arguments.
    Each may be a path to a file, '-' indicating standard input,
    or if neither of these, it is taken to be a direct input to
    the program.  It the argument is stdin or a file, each line
    is returned, one at at time to the program.  Only for text
    based data.

    A generator is returned, yielding lines of the file like things
    (whitespace stripped), or values for args if they can not be
    evaluated as files.

    These flags control the fine points:
        CL_INPUT_ALLOWED - If False, args which are not valid paths
                           or '-' will raise an error.  If False
                           and STDIN_ALLOWED is True if there are
                           zero args stdin will be read.
        STDIN_ALLOWED    - If False, '-' will not open sys.stdin.
                           '-' will be interpreted as a non-file arg
        PATHS_ALLOWED    - If False, args will not be interpreted as
                           path names.
        UNIVERSAL_NL     - True by default, open files in universal
                           newline mode. 
    Caller can set CGI_OVERRIDE and cgiParam to allow data to come
    via the CGI interface rather than file or stdin.  When the input
    is from cgiParam, extraCGIDelimChars can be set to allow parsing
    on non-whitspace.
    """

    if CGI_OVERRIDE:
        import cgi,cgitb; cgitb.enable()
        import re
        fieldExp = re.compile('[\s%s]' % extraCGIDelimChars)
        fields = cgi.FieldStorage()
        cgiValues= fields.getlist(cgiParam)
        for value in cgiValues:
            for item in fieldExp.split(value):
                yield item


    else:
        if len(args)==0 and STDIN_ALLOWED and not CL_INPUT_ALLOWED:
            args = ('-',)
        for arg in args:
            if arg == '-' and STDIN_ALLOWED:
                while True:
                    l = sys.stdin.readline()
                    if l == '':
                        break
                    else:
                        yield l.rstrip().lstrip()

            elif (os.path.exists(arg) and PATHS_ALLOWED) or \
                 (PATHS_ALLOWED and not CL_INPUT_ALLOWED) :
                # that is kind of ugly, bc I want an
                # exception here if command line 
                # arguments aren't allowed but the arg 
                # doesn't translate to a path name that works
                if UNIVERSAL_NL:
                    f = file(arg,'rU')
                else:
                    f= file(arg,'r')
                    
                while True:
                    l = f.readline()
                    if l == '':
                        break
                    else:
                        yield l.rstrip().lstrip()

            elif CL_INPUT_ALLOWED :
                yield arg
            else:
                raise Exception


# Lempel-Ziv-Welch compression algorithm
 
def compress(uncompressed):
    """Compress a string to a list of output symbols."""
 
    # Build the dictionary.
    dict_size = 256
    dictionary = {}
    for i in range(dict_size):
        dictionary[chr(i)] = chr(i)
 
    w = ''
    result = []
    for c in uncompressed:
        wc = w + c
        if wc in dictionary:
            w = wc
        else:
            result.append(dictionary[w])
            # Add wc to the dictionary.
            dictionary[wc] = dict_size
            dict_size += 1
 
            w = c
    # Output the code for w.
    result.append(dictionary[w])
    return result
 
def decompress(compressed):
    """Decompress a list of output ks to a string."""
 
    # Build the dictionary.
    dict_size = 256
    dictionary = {}
    for i in range(dict_size):
        dictionary[chr(i)] = chr(i)
 
    w = result = compressed[0]
    for k in compressed:
        if k in dictionary:
            entry = dictionary[k]
        elif k == len(dictionary):
            entry = w + w[0]
        else:
            raise ValueError, 'Bad compressed k: %s' % k
        result += entry
 
        # Add w+entry[0] to the dictionary.
        dictionary[dict_size] = w+entry[0]
        dict_size += 1
 
        w = entry
    return result

    
    
#
# Eaiser temp files
#
def mystemp(**args):
    """Like tempfile.mkstemp but returns:
    (fileobject, path) not (filedescriptor,path)

    tempfile.mkstemp documentation:

    mkstemp(suffix='', prefix='tmp', dir=None, text=False)
    mkstemp([suffix, [prefix, [dir, [text]]]])

    User-callable function to create and return a unique temporary
    file.  The return value is a pair (fd, name) where fd is the
    file descriptor returned by os.open, and name is the filename.
    
    If 'suffix' is specified, the file name will end with that suffix,
    otherwise there will be no suffix.
    
    If 'prefix' is specified, the file name will begin with that prefix,
    otherwise a default prefix is used.
    
    If 'dir' is specified, the file will be created in that directory,
    otherwise a default directory is used.
    
    If 'text' is specified and true, the file is opened in text
    mode.  Else (the default) the file is opened in binary mode.  On
    some operating systems, this makes no difference.
    
    The file is readable and writable only by the creating user ID.
    If the operating system uses permission bits to indicate whether a
    file is executable, the file is executable by no one. The file
    descriptor is not inherited by children of this process.
    
    Caller is responsible for deleting the file when done with it.

"""
    
    tmp = tempfile.mkstemp(**args)
    return (os.fdopen(tmp[0],'w'),tmp[1])



class SparseRowMapper(object):
    """A class to return filled sparse rows
    """

    def __init__(self,colOrder,missingValue):
        """Set up a new instance that knows the colunm order.  
        """
        self.blankRow = []
        self.colIdx={}
        for i,c in enumerate(colOrder):
            self.blankRow.append(copy.copy(missingValue))
            self.colIdx[c]=i

    def __call__(self,rowDict):
        """
        """
        rv = copy.copy(self.blankRow)
        #rv=self.blankRow
        for c,v in rowDict.items():
            rv[self.colIdx[c]] = v
        return rv
        

class MultiDictCluster(dict):
    """A cluster formatted representation of multiple dictionaries.
    Gene and Experiemnt weights are stored in the cluster, and used
    in clustering calculations, but extracted values are not weighted.

    The structure of the dictionary is genecentric.
    """

    def __init__ (self,missingValue=''):

        self.gWeights={}
        self.gOrder=[]
        self.eWeights = {}
        self.eOrder = []
        self.missingValue=missingValue

        self._record=None
        
        self.gNames=AutoDict()
        #self.eWeights=[]
        #self.eNames=[]

    def subSetCluster(self, genes=None, exps=None, expsFirst=True):
        """Return a subset copy of this cluster. If genes is not
        None only that genes will be included in the subset. If
        exps is not None only that experiments will be included in
        the subset.
        `genes` and `exps` can be a callable returning boolean.
        
        Arguments:
        - `genes`: Restrict the genes to this list.
        - `exps`: Restrict the experiments to this list.
        - `expsFirst`: Defines which list is processed first,
        experiments or genes.
        """
        def processGenes(o, newCluster, genes, expsFirst):
            if callable(genes):
                func = genes
            else:
                l = getIterable(genes)
                if len(l) > 0:
                    func = lambda x: x[0] in l
                else:
                    func = lambda x: x
            if expsFirst:
                allowedExps = newCluster.eOrder
            else:
                allowedExps = o.eOrder
            items = [(i, o.gWeights[i],
                      filterDictByKeys(o[i], allowedExps)) for i in o.gOrder]
            for g, w, d in filter(func, items):
                    newCluster[g] = d
                    newCluster.gWeights[g] = w

        def processExps(o, newCluster, exps, expsFirst):
            # TODO: Check if expsFirst is False
            if callable(exps):
                for e, w, d in filter(func, o.iterExpts()):
                    newCluster.eOrder.append(e)
                    newCluster.eWeights[e] = w
                    for g in d.iterkeys():
                        if g not in newCluster:
                            newCluster.gWeights[g] = o.gWeights[g]
                        newCluster[g] = {}
                        newCluster[g][e] = o[g][e]    
            else:
                for ex in exps:
                    d = o.getExpt(ex)
                    if d is not None:
                        newCluster.eOrder.append(d[0])
                        newCluster.eWeights[d[0]] = d[1]
                        for g, v in d[-1].iteritems():
                            newCluster[g] = {}
                            newCluster[g][d[0]] = v                        
        
        newCluster = self.__class__(self.missingValue)
        if expsFirst:
            processExps(self, newCluster, exps, expsFirst)
            processGenes(self, newCluster, genes, expsFirst)
        else:
            processGenes(self, newCluster, genes, expsFirst)
            processExps(self, newCluster, exps, expsFirst)
        return newCluster

    def __setitem__ (self,key,value):
        if key not in self.gWeights:
            self.gAppend(key)
        super(MultiDictCluster,self).__setitem__(key,value)

    def __delitem__(self,key):
        if key not in self.gWeights:
            raise KeyError
        else:
            del self.gWeights[key]
            self.gOrder.remove(key)
            super(MultiDictCluster,self).__delitem__(key)
            
    def renameExps (self,nameMap):
        """rename exps with the nameMap.
        nameMap <-- { oldName1: newName1,
                      oldName2: newName2,
                      ... }
        """
        for i,n in enumerate(self.eOrder):
            if n in nameMap:
                nName = nameMap[n]
                self.eOrder[i]=nName
                self.eWeights[nName]=self.eWeights[n]
                del self.eWeights[n]

                for g in self.gOrder:
                    if n in self[g]:
                        self[g][nName] = self[g][n]
                        del self[g][n]

    def renameGenes (self,nameMap):
        """rename genes with the nameMap.
        nameMap <-- { id1: newName1,
                      id2: newName2,
                      ... }
        """
        self.gNames.update(nameMap)

##         origNames = self.keys()
##         dropNames = []
##         for n in origNames:
##             if n in nameMap:
##                 nName = nameMap[n]
##                 self.gAppend(nName,self.gWeights[n])
##                 self[nName]=self[n]
##                 dropNames.append(n)
##         for n in dropNames:
##             del self[n]

    def normalizeGenes(self, s=1):
        """nomalize to n observationssuch that sum of
        values for gene is s
        """
        for k,v in self.items():
            self[k]=normDict(self[k],s)


    def normalizeExpts(self, s=1):
        """nomalize to n observationssuch that sum of
        values for gene is s
        """
        s=float(s)
        eSums=list(self.dataArray().sum(axis=0))
        nf=dict(zip(self.eOrder,[s/x for x in eSums]))
        
        for g,row in self.items():
            #print row
            for e,v in row.items():
                self[g][e]=self[g][e] * nf[e]

    def logTransform(self, base, zeroValue=0):
        """Transforms the values to its log of base `base`. If the
        value is 0 it will be replaced with `zeroValue`.

        Parameters:
        - `base`: Integer or the letter 'e'. Base of the log to
        transform each value.
        - `zeroValue`: Value to replace the zeros.
        """
        if base == 'e':
            mlog = lambda x: math.log(x)
        else:
            mlog = lambda x: math.log(x, base)

        for g, row in self.items():
            for e, value in row.items():
                if value == 0:
                    self[g][e] = zeroValue
                else:
                    self[g][e] = mlog(value)
        

    def gAppend(self,name,gWeight=1.0):
        """append a gene and weight to cluster.
        data for the gene is not set.
        """
        self.gOrder.append(name)
        self.gWeights[name]=gWeight
        self[name]={}

    def refresh_gOrder(self):
        self.gOrder=self.gWeights.keys()


    def iterExpts(self):
        """Return experiments like:
        (eName,eWeight, { geneID_1:value_1,
                           geneID_2:value_2,
                           ...}
         )
        Missing values are not filled in.  Weights are not applied.
        """
        for e in self.eOrder:
            yield (e, self.eWeights[e],
                   dict([(k,v[e])
                         for k,v in self.items() if e in v ]))


    def getExpt(self,name):
        """Return an experiment like:
        (eName,eWeight, { geneID_1:value_1,
                          geneID_2:value_2,
                           ...})
        Missing values are not filled in.  Weights are not applied.
        name can be an iterable for a group of experiments.
        """
        
        if type(name) in StringTypes:
            names=(name,)
        else:
            names = name
        rv=[]
        for name in names:
            if name in self.eWeights:
                rv.append((name, self.eWeights[name],
                           dict([(k,v[name])
                                 for k,v in self.items() if name in v ])))
            else:
                rv.append(None)
                
        if len(rv) >1:
            return rv
        elif len(rv)==1:
            return rv[0]
        else:
            return None
            
    
    
    def eAppend(self,expt,name,eWeight=1.0,pctNorm=False,sumNorm=False):
        """append a new column ("experiment" in cluster speak)
        from expt (a dict). Expt can be iterable of dicts/experiments.
        'append' is an alias for eAppend.
        """
        
        try:
            expt.keys()
            expt=[expt]
            name=[name]
        except:
            pass

        for e,n in zip(getIterable(expt),getIterable(name)):
            #print e,n
            self.eOrder.append(n)
            #print type(name),name
            self.eWeights[n]=eWeight
            
            if pctNorm and sumNorm:
                raise ValueError, "only one, at most, of pctNorm and sumNorm can be True"
            elif pctNorm:
                nf = 100.0/sum((v for v in e.values() if v is not None))
            elif sumNorm:
                nf=1.0/sum((v for v in e.values() if v is not None))
            else:
                nf=1.0

            for k,v in e.items():
                if v is not None:
                    if k not in self:
                        self[k]={}
                        self.gWeights[k]=1.0

                    self[k][n]=v*nf
                    
    # backward compatability
    append=eAppend
    

    def update (self,other):
        """Add genes in one cluster to another cluster with the same
        experiments.  Other's genes clobber self's in conflicts.
        """
        if (self.eWeights != other.eWeights):
            raise ValueError, "Cluster experiments are not the same"
        else:
            #for g,v in other.items():
                #self[g]=v
                #self.gWeights[g]=other.gWeights[g]
            super(MultiDictCluster,self).update(other)
            self.gWeights.update(other.gWeights)
                
            
    def addNestedDict(self,dod,eCallback=None,**kwArgs):
        """Add dict of dict. keys are experiment names and values
        are experiment gene:value pairs. eWeight and pctNorm passwd
        as keyword arguments will be passed to append.  eCallback is
        preformed on the values before appending.
        """
        for n,e in dod.items():
            if eCallback != None:
                e=eCallback(e)
            self.eAppend(e,n,**kwArgs)

    def gValues(self,g,missingValue=None):
        """returns a list of values for a gene in experiemnt order
        with default values filled in.  If missingValue is not none
        the value overrides self.missingValue.
        """
        rv =[]
        eIdx={}

        if missingValue == None:
            missingValue = self.missingValue

        for i,e in enumerate(self.eOrder):
            if e in self[g]:
                rv.append(self[g][e])
            else:
                rv.append(missingValue)
             
        return rv


    def gPresValues(self,g):
        """list of values for gene in no particular order and
        with on ly the actual values (no missing fill in)
        """
        return self[g].values()


    def gSort(self,keyFcn):
        """
        """
        pass


    def eSort(self,keyFcn):
        """
        """
        pass


    def clusterFmt(self,outFile=sys.stdout,):
        """
        """

        if len(self.eOrder) == 0:
            self.eOrder=sorted(list(self.experiments()))
        if len(self.eWeights) == 0:
            self.eWeights = dict(zip(self.eOrder,[1.0]*len(self.eOrder)))

        #print self.eOrder
        #print self.eWeights.keys()
        
        print >> outFile, '\t'.join(['UID','NAME','GWEIGHT']+
                                    [str(x) for x in self.eOrder])
        print >> outFile, '\t'.join(['EWEIGHT','','']+
                                    [str(self.eWeights[e]) for e in self.eOrder])
        for k in self.gOrder:
            fields= ['"%s"'%k,self.gNames[k],self.gWeights[k]]
            for innerK in self.eOrder:
                try:
                    fields.append(self[k][innerK])
                except KeyError:
                    fields.append(self.missingValue)
            
            print >> outFile, '\t'.join([str(x) for x in fields])

    def setGName(self, k, name):
        self._gNames[k] = name
            
    def kSqlNameMap (self,ksqlclass,lookupFailPrefix='ID: ',
                     objectCallback=lambda x: str(x),
                     dim='g'):
        """Map gNames to kSqlObject string representations (usually Name column).
        keys must be
        dim sould be 'g' for gene operation or 'e' for expt operation
        """
        if dim.startswith('g'):
            nameMap={}
            goCopy = self.gOrder[:]
            for g in goCopy:
                try:
                    newG = objectCallback(ksqlclass(int(g)))
                except:
                    newG="%s%s" % (lookupFailPrefix,g)
                nameMap[g]=newG
            self.renameGenes(nameMap)

        elif dim.startswith('e'):
            nameMap={}
            for e in self.eOrder:
                try:
                    newE=objectCallback(ksqlclass(int(n)))
                except:
                    newE="%s%s" % (lookupFailPrefix,n)
                nameMap[e] =newE

            self.renameExps(nameMap)
                    
        else:
            raise ArgumentError, "dim must be 'g' or 'e'"


    def __add__ (self,other):
        """combine clusters, return copy
        """
        rv=self.__class__()
        #rv=copy.deepcopy(self)
        for n,w,e in self.iterExpts():
            rv.append(e,n,eWeight=w)
        
        for n,w,e in other.iterExpts():
            rv.append(e,n,eWeight=w)
        return rv

        
    def gStats(self,missingValue=0.0):
        """dict of {geneID: (min,max,mean,median,std,stderr,
        Shapiro-Wilk(w,p),normaltest_chisq (D'Agostino and Pearson),...}
        """
        import scipy as S
        import scipy.stats as SS
        
        rv={}
        for k,v in self.items():
            #print k,v
            va = S.array(self.gValues(k,missingValue))

            try:
                normaltest = SS.normaltest(va)
            except:
                normaltest=None
            try:
                shapiro = SS.shapiro(va)
            except:
                shapiro = None
            

            try:
                rv[k] = (va.min(),va.max(),va.mean(),SS.median(va),
                         SS.std(va),SS.stderr(va),normaltest,
                         shapiro
                         )
            except:
                print k,va
                raise
        return rv


    def gTransform(self,fcn):
        """transform values using a callback (fcn).
        fcn is called like so:  fcn(geneID,[exp1:value1, ... ])
        And should return an iterable of exp:value pairs,
        or may return None in which case the gene is removed
        from the cluster
        """
        for k,v in self.items():
            #print k,v
            eIDs,values= unzip(v.items())
            txValues=fcn(k,v)
                         #values)
            if txValues == None:
                del self[k]
            else:
                for e,val in txValues.items():
                    self[k][e] = val

    def transform(self,fcn):
        """transform every value in the cluster with
        fcn.
        """
        for g in self:
            for e,v in self[g].items():
                self[g][e]=fcn(v)
                

    def allValues(self,dtype='d'):
        """returns a 1D array of all the values
        no missing value fillin
        """
        import numpy as N
        c = 0
        for g in self:
            c+=len(self[g])
        rv=N.zeros(c,dtype=dtype)
        i=0
        for g in self:
            for v in self[g].values():
                rv[i]=v
                i+=1
        return rv

    

    def filterGenes(self,fcn,keep=True):
        """The callback is called like fcn(self,geneID), and
        is cast as a bool. The value of keep (True of False),
        controls which genes are kept (others are chucked).
        """
        for k in self.keys()[:]:
            cv = fcn(self,k)
            if bool(cv) != keep:
                del self[k]
        

    def dataArray(self,dtype='d'):
        """Returns a numpy array (with 0's for
        missing values.  The raw values are used to
        populate the array; weights are not applied.
        """
        
        import numpy as N
#        self.refresh_gOrder()

        genesOrdered = sorted(self.gOrder)
        gIdx = {}
        for i,g in enumerate(genesOrdered):
            gIdx[g]=i
            
        eIdx = {}
        for i,e in enumerate(self.eOrder):
            eIdx[e]=i
       
        
        a=N.zeros((len(self.gOrder),len(self.eOrder)),dtype=dtype)
        for g in genesOrdered:
            for e in self.eOrder:
                if e in self[g]:
                    a[gIdx[g],eIdx[e]] = self[g][e]
            
        return a


    def runCluster(self,jobName,
                   gDist='c',gMethod='m',
                   eDist='c',eMethod='m',
                   addScaleGene=None,
                   **kwds):
        """
        Run Pycluster:
        transpose: if equal to 0, genes (rows) are clustered;
                      if equal to 1, microarrays (columns) are clustered.
           dist     : specifies the distance function to be used:
                      dist=='e': Euclidean distance
                      dist=='b': City Block distance
                      dist=='c': Pearson correlation
                      dist=='a': absolute value of the correlation
                      dist=='u': uncentered correlation
                      dist=='x': absolute uncentered correlation
                      dist=='s': Spearman's rank correlation
                      dist=='k': Kendall's tau
           method   : specifies how the distance between two clusters is defined:
                      method=='a': the distance between the arithmetic means of the
                                   two clusters
                      method=='m': the distance between the medians of the two
                                   clusters
                      method=='s': the smallest pairwise distance between members
                      method=='x': the largest pairwise distance between members of
                                   the two clusters
                      method=='v': average of the pairwise distances between
                                   members of the clusters"""
        rec=self.record(clobber=True)
        gTree = None
        eTree = None
        if gDist != None and gMethod != None:
            gTree=rec.treecluster(dist=gDist,method=gMethod)
        if eDist != None and eMethod != None:
            eTree=rec.treecluster(dist=eDist,method=eMethod,
                                   transpose=1)
        rec.save(jobName,gTree,eTree)
        if addScaleGene != None:
            cf=file(jobName + '.cdt','a')
            name,values=addScaleGene
            print >> cf, '\t'.join(['GENE%dX'%len(self),name,name,'1.0']+[str(x) for x in values])
            cf.close()
        
    #def array

    def record(self,clobber=False):
        """return the corresponding Pycluster.Record object.
        This is cached after the first call. Clobber causes
        recalculation .
        """

        if clobber or self._record == None:
            import Pycluster
            buf=StringIO()
            self.clusterFmt(buf)
            buf.seek(0)
            self._record = Pycluster.Record(buf)

        return self._record

    def eDelete(self,e):
        """delete the experiment specified by name
        """
        self.eOrder.remove(e)
        del self.eWeights[e]

        for g in self.keys():
            if e in self[g]:
                del self[g][e]

    @classmethod
    def readCdt(cls,cdtFile):
        """
        """
        rv=cls()
        if type(cdtFile) in StringTypes:
            cdtFile=open(cdtFile)

        lines = cdtFile.__iter__()
        
        rv.eNames=lines.next().strip().split('\t')[4:]
        lines.next()
        rv.eWeights=dict(zip(rv.eNames,
                             [float(x) for x in lines.next().strip().split('\t')[4:]]))
        for l in lines:
            f = l.strip().split('\t')
            gID=f[1]
            rv.gNames[gID]=f[2]
            rv.gOrder.append(gID)
            rv.gWeights[gID] = float(f[3])
            values = [float(x) for x in f[4:]]
            for i,eName in enumerate(rv.eNames):
                if gID not in rv:
                    rv[gID] = {}
                rv[gID][eName]=values[i]
        return rv
        
            
    def setGEV(self,*tripples):
        """Add/set (Gene, Experiement, Value) tripples
        EWIEGHT and GWEIGHT are set to 1.0 if exp or gene are
        new to the cluster.
        """
        #print tripples
        for g,e,v in tripples:
            e=str(e)
            if g not in self:
                self[g]={}
                self.gWeights[g]=1.0
            if e not in self.eWeights:
                self.eOrder.append(e)
                self.eWeights[e]=1.0
                
            self[g][e]=v


    def experiments(self):
        """return set of experiment names
        """
        rv = set([])
        for gene in self.keys():
            for exp in self[gene].keys():
                rv.add(exp)
        return rv


    def minValue(self):
        """return minimum non-missing value
        """
        if len(self) == 0:
            raise ValueError, "No values found in cluster."
        
        rv=sys.maxint
        for rowD in self.itervalues():
            rMin=min(rowD.itervalues())
            if rMin < rv:
                rv=rMin
        return rv   

        
class PickleableObject(object):
    """Generic pickleable class
    """
    def save(self,filename,clobber=False):
        """given an object, e.g. a taxon dictionary,
        pickle it into the specified filename

        Set clobber = True to allow overwriting of the
        output file"""
        f = safeOFW(filename,clobber=clobber,append=False)
        try:
            cPickle.dump(self,f)
        except:
            os.unlink(filename)
            raise
        else:
            f.close()
        

    @classmethod
    def restore(cls,filename):
        """given a filename, return an object containing the
        pickled data in that file"""
        f = open(filename,'rb')
        obj = cPickle.load(f)
        return obj

 
def vec2col(vec):
    """reshape a 1d numpy vector to a single 2d column
    """
    vec.shape=vec.shape[0],1
    return vec

def vec2row(vec):
    """reshape a 1d numpy vector to a single 2d row
    """
    
    vec.shape=1,vec.shape[0]
    return vec

def getIterable(obj):
    if obj is None:
        return []
    try:
        obj.__iter__
    except AttributeError:
        return (obj,)
    return obj        

def partition(func, iterator):
    """Returns a list of elements for which `func` returns True and
    another list of elements that returns False.
    
    Arguments:
    - `func`: Function that returns boolean.
    - `iterator`: Iterator of elements.
    """
    trues, falses = [], []
    for item in iterator:
        if func(item):
            trues.append(item)
        else:
            falses.append(item)
    return trues, falses

def filterDictByKeys(dic, allowedKeys):
    return dict([(k,v) for k,v in dic.iteritems() \
                 if k in allowedKeys])

def lastLine(f):
    """given a path or file returns last line of it
    if path is given file is opend and closed. if file
    is given is repositioned to its postion on calling
    (unless there is an exception)
    """

    if type(f) in StringTypes:
        ipos = None
        f=file(f)
    else:
        ipos=f.tell()

    f.seek(-2,2)
    while f.read(1) != '\n':
        f.seek(-2,1)
    rv=f.readline()
    if ipos==None:
        f.close()
    else:
        f.seek(ipos,0)
    return rv

class FieldConverter(list):
    """class that makes a callabble converter object.
    takes a sequence of fields and performs a position
    specific conversion on each.

    conversion callables are avalible using a list interface.
    default conversion can be set on construction.

    fast construction example:
    FieldConverter((x,y),int),((a,b),datetime.datetime.fromtimestamp),default=str)

    default is str by default

    you could have a strict converter by enumerating all the fields and
    setting the default to a callable that raises the Exception of your choice


    """
    def __init__ (self,*index_fcn_tuples,**kwds):
        if 'default' in kwds:
            self.default=kwds['default']
        else:
            self.default=str

        tmp=[]
        for idxs,fcn in index_fcn_tuples:
            tmp.extend(getIterable(idxs))
        self.extend([self.default]*(max(tmp)+1))
        #print self

        for idxs,fcn in index_fcn_tuples:
            for i in getIterable(idxs):
                self[i]=fcn

    def __call__(self, *fields):
        """fields can be a sequence or indivudual arguments.
        extra fields are returned after conversion with default
        callable
        """
        print fields
        if len(fields) ==1 and len(self) > 1 :
            fields=fields[0]
        rv=[]
        for i,f in enumerate(fields):
            if i < len(self):
                rv.append(self[i](f))
            else:
                rv.append(self.default(f))
        return tuple(rv)
        

class SafeDict(dict):
    """A dictionary with a global default value (None by default).
    Set default at creation time with SafeDict(default=x).
    """
    def __init__(self, *args, **kwargs):
        if 'default' in kwargs:
            self.__default = kwargs['default']
            del kwargs['default']
        else:
            self.__default = None
        super(SafeDict, self).__init__(*args, **kwargs)

    def __getitem__(self, key):
        if key in self:
            return super(SafeDict, self).__getitem__(key)
        else:
            return self.__default


class AutoDict(dict):
    """A safe dictionay, where the key value is returned if
    a lookup fails. e.g.:

    myDict= AutoDict()
    myDict[4] # <-- returns 4
    """
    def __getitem__(self, key):
        if key in self:
            return super(AutoDict, self).__getitem__(key)
        else:
            return key
   



def wrap(text, width):
    """
    A word-wrap function that preserves existing line breaks
    and most spaces in the text. Expects that existing line
    breaks are posix newlines (\n).
    """
    return reduce(lambda line, word, width=width: '%s%s%s' %
                  (line,
                   ' \n'[(len(line)-line.rfind('\n')-1
                          + len(word.split('\n',1)[0]
                                ) >= width)],
                   word),
                  text.split(' ')
                  )
