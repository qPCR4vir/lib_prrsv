#!/usr/local/bin/python 
#
# Common iterator tricks and batching tools
#

import sys

from types import ListType,TupleType,GeneratorType,StringTypes
import os, re, random
import bisect

__verison__ = tuple([int(x) for x in
                     '$Revision: 1.6 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

def getIteratable(obj):
    """
    Return something that has an __iter__ method,
    if obj does not provide that return the
    one-tuple (obj,).

    Arguments:
    - `ojb`:
    """
    try:
        obj.__iter__
    except AttributeError:
        return (obj,)
    else:
        return obj


def fileIterator(*files):
    """returns a generator of open files for reading.  Similar
    to multiFile, which returns a generator of lines from
    multiple files.

    Input (files) can be a single list or comma seperated
    values of ether files paths and of file objects.
    Nested tuples, lists and generators are unfolded.  For
    this reason, the order of files may not be preserved.

    If any 'file' is a path it will be opened and closed
    If any 'file' is a is not a path, seek(0) will be
    performed if the 'file' supports it, after reading to EOF,
    but not before.
    """
    try:
        l = len(files)
    except AttributeError:
        files = [files]
        
    if len(files)==1:
        if type(files[0]) in (TupleType,ListType):
            files = files[0]
    else:
        files = reversed(files)
    files=list(files)
        
    while len(files) > 0:
        f=files.pop(0)
        close=False
        seek0=False
        if hasattr(f,'read'):
            inFile=f
            seek0=True
            
        if type(f) in (TupleType,ListType,GeneratorType):
            files.extend(f)
            continue
        
        if type(f) in StringTypes:
            if f == '-':
                inFile = sys.stdin
            else:
                inFile=file(f,'r')
                close = True

        yield inFile
        # clean up
        if close:
            inFile.close()
        elif seek0:
            try:
                inFile.seek(0)
            except:
                pass


def multiFile(*files):
    """returns a generator of lines from many files.
    Input (files) can be a single list or comma seperated
    values of ether files paths and of file objects.
    Nested tuples, lists and generators are unfolded.  For
    this reason, the order of files may not be preserved.

    If 'file' is a path it will be opened and closed
    If 'file' is a is not a path, seek(0) will be performed,
    after reading to EOF, but not before.
    """
    for inFile in fileIterator(*files):
        for l in inFile:
            yield l


def iterCount(iterable):
    """count number of times each value is seen in iterable.
    return dict((value,count),...)
    """
    rv = {}
    for x in iterable:
        if x in rv:
            rv[x]+=1
        else:
            rv[x]=1
    return rv


def batchList(iterable,batchSize=None,batchCount=None):
    """ return a list of lists representing batches of items from
    iterable. Either batchSize or batchCount should be an integer.
    """
    if batchCount == None and batchSize == None:
        raise ValueError , "only batchSize or batchCount may be specified"
    elif batchCount != None and batchSize != None:
        raise ValueError , "only batchSize or batchCount may be specified"

    rv=[]
    try:
        inList=list(iterable)
    except TypeError:
        # got non sequence type
        return [[iterable]]
        
    if batchSize != None:
        batchSize = int(batchSize)
    else:
        if len(inList)%batchCount==0:
            batchSize=len(inList)/batchCount
        else:
            batchSize=(len(inList)/batchCount) +1

    for start in range(0,len(inList),batchSize):
        rv.append(inList[start:start+batchSize])
    return rv


def batchIndicies(totalSize,batchSize=None,batchCount=None):
    """ return a list of lists if (startIdx, endIdx) representing batches
    from iterable. Either batchSize or batchCount should be an integer.
    """
    rv=[]
    if batchCount == None and batchSize == None:
        raise ValueError , "only batchSize or batchCount may be specified"
    elif batchCount != None and batchSize != None:
        raise ValueError , "only batchSize or batchCount may be specified"

    if batchSize != None:
        batchSize = int(batchSize)
    else:
        if len(inList)%batchCount==0:
            batchSize=totalSize/batchCount
        else:
            batchSize=(totalSize/batchCount) +1

    for start in range(0,totalSize,batchSize):
        if start+batchSize > totalSize:
            rv.append((start,totalSize))
        else:
            rv.append((start,start+batchSize))
    return rv

def randomLineIterator(fname, numResults):
    """Returns a generator of `count` random lines from `fname`
    """
    visitedRanges = {}
    starts = []
    avgLineSize = 2
    lineRe = re.compile(r'\n([^\n]*)\n', re.M)
    lineReStart = re.compile(r'^([^\n]*)\n', re.M)
    lineReEnd = re.compile(r'\n([^\n]*)\s?$', re.M)
    def nextRandomPosition(maxPos):
        while True:
            pos = random.randint(0, maxPos)
            idx = bisect.bisect_right(starts, pos)
            if idx > 0 and visitedRanges[starts[idx-1]] > pos:
                continue
            if idx > 0 and idx < len(starts):
                lo = starts[idx-1]
                hi = visitedRanges[lo]
                if lo <= pos <= hi:
                    continue

            return pos

    def lineInPosition(f, pos):
        raw = ''
        startP, endP = pos, pos
        os.lseek(f, pos, os.SEEK_SET)
        halfSize = avgLineSize / 2
        endFile = False
        while True:
            if startP > 0:
                startP -= halfSize
                frontChunk = halfSize
                if startP < 0:
                    frontChunk = halfSize + startP
                    startP = 0
                os.lseek(f, startP, os.SEEK_SET)
                chunk = os.read(f, frontChunk)
                raw = chunk + raw
            os.lseek(f, endP, os.SEEK_SET)
            chunk = os.read(f, halfSize)
            raw += chunk
            endP += len(chunk)

            endFile = len(chunk) < halfSize
            if startP == 0:
                match = lineReStart.search(raw)
            elif endFile:
                match = lineReEnd.search(raw)
            else:
                match = lineRe.search(raw)
            if not match:
                continue

            line = match.group(1)
            innerStart = raw.find(line)
            sPos = startP + innerStart
            ePos = sPos + len(line)
            if pos > ePos:
                raw = raw[ePos-startP:]
                startP = ePos
                continue
            elif pos < sPos:
                raw = raw[:sPos-startP]
                endP = sPos
                continue
            if sPos > 0:
                sPos -= 1
            ePos += 1
            return (line, sPos, ePos)
        raise Exception

    f = os.open(fname, os.O_RDONLY)
    fsize = os.fstat(f).st_size
    curResults = 0
    while curResults < numResults:
        if len(starts) == 1 and \
               visitedRanges[starts[0]] - starts[0] >= fsize:
            raise StopIteration()

        curP = nextRandomPosition(fsize)
        line, sPos, ePos = lineInPosition(f, curP)
        idx = bisect.bisect_right(starts, sPos)
        lo = None
        if idx > 0:
            lo = starts[idx-1]
            hi = visitedRanges[lo]
            if lo <= sPos <= hi:
                idx -= 1
                visitedRanges[lo] = ePos
            else:
                lo = sPos
                starts.insert(idx, lo)
                visitedRanges[lo] = ePos
        else:
            lo = sPos
            starts.insert(idx, sPos)
            visitedRanges[sPos] = ePos
        if idx < len(starts)-1:
            lo2 = starts[idx + 1]
            hi2 = visitedRanges[lo2]
            if lo2 <= ePos <= hi2:
                visitedRanges[lo] = hi2
                del visitedRanges[lo2]
                del starts[idx + 1]

        curResults += 1
        yield line

    os.close(f)
