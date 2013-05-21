"""Annotate Module
GFF reference: http://www.sequenceontology.org/gff3.shtml
"""
__version__ = tuple([int(x) for x in
                     '$Revision: 1.4 $'.split()[1].split('.')])
__author__ = "Kael Fischer"


import types
import urllib



class GffFile (object):
    """A GFF file container.
    Facilities for parsing GFF files.

    All information about relationships between features is contained
    in this class.

    properties:
    meta = dictionary of gff meta info
    features = list of features (in file order)
               note that features have the .index property which
               is the position in this list
    IDidx = a dictionary of Feature.ID:list position
    children = dict ParentID:[child Features]
    parent = dict Feautre:Parent.ID

    [note that parents are designanetd by ID and children are
    GffFeature objects]

    pos = iterator pointer

    """
    def __init__(self,inFile=None):
        """ Make new instance, parse file (or file at path) if present
        """
        self.meta = {}
        self.features = []
        self.IDidx={}

        self.children = {}
        self.parent = {}

        self.pos = 0

        if inFile != None:
            self._parseFile(inFile)


    def forward(self):
        """return next feature
        """
        self.pos+=1
        return self.currentFeature()


    def forwardType(self):
        """Move forward to the next record which is the sams as the
        current type. or None if none found.
        """
        t = self.currentFeature().type
        self.forward()
        for f in self:
            if f.type == t:
                return f
        else:
            return None


    def back(self):
        """return last feature
        """        
        self.pos-=1
        return self.currentFeature()


    def backType(self):
        """Move backward to the last record which is the sams as the
        current type. or None if none found.
        """
        t = self.currentFeature().type
        self.back()
        for f in self.__back_iter__():
            if f.type == t:
                return f
        else:
            return None

    
    def currentFeature(self):
        """return self.features[self.pos] or None if out of range 
        """
        try:
            return self.features[self.pos]
        except IndexError:
            return None


    def typeGenerator(self,featureType=None,startFrom=None,direction=1):
        """generator over features of a patricular type (str or feature instance),
        starting from position specified (int or feature) or if None from
        current pos.  pos is not effected and multiple generators have independant
        pointers.
        direction can be set to >=0 (forward) or < 0 (backwards)
        """
        if direction <0:
            step = -1
        else:
            step = 1
            
        if startFrom == None:
            startFrom = self.pos

        if type(startFrom) not in (types.IntType, types.LongType):
            startFrom=startFrom.index
        
        if type(featureType) not in types.StringTypes:
            if featureType == None:
                featureType = self.currentFeature().type
            elif isinstance(featureType,GffFeature):
                featureType = featureType.type
            else:
                raise ValueError, "Feature type undetermined"
                    
        while True:
            try:
                if self.features[startFrom].type == featureType:
                    startFrom+=step
                    yield self.features[startFrom]
            except IndexError:
                raise StopIteration
            
            

    def __iter__ (self):
        """iterates over records; changes pos
        """
        
        while True:
            yield self.features[self.pos]
            self.pos+=1
            if self.pos >= len(self.features):
                raise StopIteration


    def __back_iter__ (self):
        """iterates backwards over records; changes pos
        """
        
        while True:
            yield self.features[self.pos]
            self.pos-=1
            if self.pos < 0:
                raise StopIteration


    def __getitem__(self,ID):
        """lookup feature by ID
        """
        if type(ID) == types.IntType:
            return self.features[ID]
        else:
            return self.features[self.IDidx[ID]]


    def seek(self,pos):
        """move pos to a particular rec (int, or feature)
        """
        if type(pos) == types.IntType:
            self.pos=pos
        else:
            self.pos = pos.index
        return self.features[self.pos]


    def _parseFile(self,inFile):
        """dirty file parser
        """

        if type(inFile) in types.StringTypes:
            inFile=file(inFile)

        i=0
        for l in inFile:
            l=l.strip()
            if len(l) == 0:
                continue
            if not l.startswith('#'):
                # parse a line
                newF = GffFeature(l,i)
                self.features.append(newF)
                if newF.ID != None:
                    self.IDidx[newF.ID] = i

                if 'Parent' in newF.attributes:
                    p = newF.attributes['Parent']  # this is an ID
                    c = newF                       # child may not have
                                                   # an ID  
                    self.parent[c]=p
                    if p not in self.children:
                        self.children[p] = []
                    self.children[p].append(c)

                i+=1
            else:
                if l.startswith("###FASTA"):
                    pass
                    #do fasta stuff
                elif l.startswith("###"):
                    #close feature refs
                    pass
                elif l.startswith('##'):
                    k,v = [urllib.unquote(x)
                           for x in l[2:].split('\t') ]
                    self.meta[k]=v
                else:
                    # more that 3 #'s
                    continue


    def findByAttribute(self,criteria = {},**kwargs):
        """Find next (from pos) feature that has attributes specified.
        """
        criteria = dict(criteria.items() + kwargs.items())
        self.forward()
        for f in self:
            if f.attrMatch(criteria):
                return f


    def countByAttribute(self,criteria = {},**kwargs):
        """Global count of features with attributes specified.
        """
        criteria = dict(criteria.items() + kwargs.items())

        rv=0
        for f in self.features:
            if f.attrMatch(criteria):
                rv+=1
        return rv


    def generatorByAttribute(self,criteria = {},startFrom=None,
                             direction=1,**kwargs):
        """Generator of features with attributes specified, 
        starting from position specified (int or feature) or if None from
        current pos.
        
        Direction can be set to >=0 (forward, default) or < 0 (backwards)

        Instance's pos attr is not effected.
        """
        criteria = dict(criteria.items() + kwargs.items())

        if direction <0:
            step = -1
        else:
            step = 1
            
        if startFrom == None:
            startFrom = self.pos

        if type(startFrom) not in (types.IntType, types.LongType):
            startFrom=startFrom.index
                
        while True:
            try:
                startFrom+=step
                f=self.features[startFrom]
                if f.attrMatch(criteria):
                    yield f
            except IndexError:
                raise StopIteration
        

class GffFeature (object):
    """Feature Class
    """
    def __init__(self,line=None,index=None):

        self.seqid = None
        self.source = None
        self.type = None
        self.start = None
        self.end = None
        self.score = None
        self.strand = None
        self.phase = None
        self.attributes = {'ID':None}
        self.index = index
        if line != None:
            self.parseLine(line)
        

 
    def parseLine(self,line):
        """dirty parser
        """
        f = line.strip().split('\t')
        self.seqid = urllib.unquote(f[0])
        self.source = urllib.unquote(f[1])
        self.type = urllib.unquote(f[2])
        self.start = int(f[3])
        self.end = int(f[4])
        if f[5] == '.':
            self.score = f[5]
        else:
            self.score=float(f[5])
        if f[6] in ('+-.?'):
            self.strand = f[6]
        else:
            raise ValueError, "'%s' given for 'strand'.  Must be:+,-,., or ?" % f[6]
        if f[7] == '.':
            self.phase = '.'
        else:
            x=int(f[7])
            if x<0 or x>3:
                raise ValueError, "'%s' given for 'phase'.  Must be: .,0,1,2 or 3" % f[7]
            self.phase = x

        for attr in f[8].split(';'):
            k,v = attr.split('=')
            if ',' not in v:
                v=urllib.unquote(v)
            else:
                v=tuple([urllib.unquote(x) for x in v.split(',')])
            self.attributes[urllib.unquote(k)] = v
        
        self.ID = self.attributes['ID']


    def __str__(self):
        """returns ID or name attribute or "No ID/Name"
        """
        if self.ID == None:
            if self.attributes['Name'] == None:
                return "No ID/Name"
            return self.attributes['Name']
        return self.ID


    def __repr__(self):
        """returns a valid GFF feature line coresponding to
        the instance
        """
        def squashTupple(mbt):
            if type(mbt) != types.TupleType:
                return urllib.quote(mbt)
            else:
                return ','.join([urllib.quote(vp) for vp in mbt])

        
        attrStr = ';'.join(['%s=%s' % (urllib.quote(k),
                                       squashTupple(v))
                            for k,v in self.attributes.items() if v != None])
        return '\t'.join([urllib.quote(str(x)) for x in
                          (self.seqid,self.source,self.type,self.start,self.end,
                           self.score,self.strand,self.phase)]+[attrStr])


    def __eq__(self,feat):
        """Returns True if all elements of features are equal 
        """
        if not isinstance(feat,GffFeature):
            return False
        
        return  (self.seqid == feat.seqid and
                 self.source == feat.source and
                 self.type == feat.type and
                 self.start == feat.start and
                 self.end == feat.end and
                 self.score == feat.score and
                 self.strand == feat.strand and
                 self.phase == feat.phase and
                 self.attributes == feat.attributes )  


    def attrMatch(self,criteria = {},**kwargs):
        """returns True if features attributes match the given criteria
        """
        criteria = dict(criteria.items() + kwargs.items())

        for k,v in criteria.items():
            if k not in self.attributes:
                return False
            if v != self.attributes[k]:
                return False
            pass
        return True


    def noquoteStr(self):
        """return __repr_ with url sytle quoteing unquoted.
        """
        return urllib.unquote(self.__repr__())

    def __getitem__(self,item):
        """get item interface to instance's attributes

        i.e.:
        x['gene'] = x.attributes['gene']

        """
        return self.attributes[item]
