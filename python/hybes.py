#!/usr/local/bin/python
# prrsv hybe analysis
# Kael Fischer, et al
# 2013

import re
from copy import copy

import numpy as N
from numpy import ma as MA
from scipy import stats

from utils import MultiDictCluster

class IntensityArray(object):
    """numpy holder for intensities in a cluster file.
    I assume the data are un normalized.
    attributes are:
    - array (numpy masked array of floats, starting mask is True where
    the intensity is None - i.e flagged during gridding)
    - uidIdx (a dictionary that maps row #'s to oligo IDs)
    - uidList (an ordered list of UIDs as they appear in the data)
    - arrayIdx (a dictionary that maps column #'s to arrayIDs,
    and chip#s - integers)
    - oligoLongName (a dictionary that maps short uids to the
    oligo Names in the cluster)
    """

    def __init__ (self, clusterPath, intensityCutoff=None):
        """cluster path specifies the cluster to read in.
        """
        self.clusterPath=clusterPath
        f=file(clusterPath)
        flds=f.next().strip().split('\t')
        if flds[2]=='GWEIGHT':
            colOne=3
        else:
            colOne=2
            
        self.arrayList=flds[colOne:]
        self.arrayIdx=dict(((x,i) for (i,x) in enumerate(self.arrayList)))
        self.arrayIdx.update(dict(((int(x.split('_')[0]),i) for (i,x)
                              in enumerate(self.arrayList))))
        self.uidList=[]
        self.uidIdx={}
        self.oligoLongName={}
        self.annotatedOligos={}
        
        i=0
        data=[]

        for l in f:
            flds=l.strip().split('\t')
            if flds[0] == 'EWEIGHT':
                continue
            uid=flds[0].strip('"')
            oName=flds[1].strip('"')
            self.uidList.append(uid)
            self.uidIdx[uid]=i
            self.oligoLongName[uid]=oName
            if re.match('\d+',oName):
                self.annotate(uid,'PRRSV')
            else:
                self.annotate(uid,'NON-VIRAL')
            annFields=oName.split('::')
            if len(annFields)>1:
                for af in annFields[1:]:
                    self.annotate(uid,af)
            row=[]
            for v in flds[colOne:]:
                try:
                    row.append(float(v))
                except ValueError:
                    row.append(None)
            while len(row) < len(self.arrayIdx)/2:
                row.append(None)
            data.append(row)
            i+=1
            
        self.array=MA.masked_array(data,mask=N.equal(data,None),
                                   dtype='float',
                                   fill_value=1e-20)

        self.sumRawIntensity=self.array.sum(0)
        self.operations=[]

        self.norm='raw'
        self.otherTransforms='none'
        
        if intensityCutoff is not None:
            self.maskByIntenisty(intensityCutoff)

    def annotate(self,uid,annotation):
        """
        """
        try:
           self.annotatedOligos[annotation].add(uid)
        except KeyError:
            self.annotatedOligos[annotation]=set([uid])

    def knownAnnotations(self):
        """return a list of the known annotation groups
        """
        return self.annotatedOligos.keys()


    def maskByIntenisty(self,minIntensity):
        """mask any values less than minIntensity
        """
        self.array.mask = self.array.mask | (self.array.filled(fill_value=1e+20) < minIntensity)
        self.operations.append('Mask Intensity < %s' % minIntensity) 

    def maskAnnotatedOligos(self,annotation):
        """Mask oligos with annotation.
        """
        try:
            for row in (self.uidIdx[x] for x in self.annotatedOligos[annotation]):
                self.array.mask[row,:]=True
            self.operations.append('Mask %s Oligos' % annotation)
        except KeyError:
            raise KeyError, ("%s: unknown. The known annotations are:\n%s"%
                             (annotation,sorted(self.knownAnnotations())))

    def maskNonViral(self):
        """Mask NON-VIRAL oligos.
        """
        self.maskAnnotatedOligos('NON-VIRAL')

    def maskArrays(self,arrays):
        """Mask out arrays specified by a list (or sequence)
        of integers. 
        """
        arrays=list(arrays)
        try:
            for col in (self.arrayIdx[x] for x in arrays):
                self.array.mask[:,col]=True
            self.operations.append('Mask Arrays: %s' % arrays)
        except KeyError:
            raise KeyError, ("%s: unknown." % col)

    def log2(self):
        """Log base 2 transform the intensities.
        """
        self.array=N.log2(self.array)
        self.opperations.append('log2 transform')

       
    def sumNormalize(self):
        """divide intensities by sum of unmasked values on array.
        """
        self.norm='sum'
        self.array=self.array/(self.array.sum(0))
        self.operations.append('Sum Normalized Arrays')

    def pctNorm(self):
        """divide intensities by sum of unmasked values on an array,
        scale to 100%.
        """
        self.norm='percent'
        self.data=100.0*self.array/(self.array.sum(0))
        self.operations.append('Percent Normalized Arrays')

    def normalizeArrays(self,normFactors,description):
        """Arbitrary normalzation.  Divide arrays by normFactors,
        which shoudl be an array the length of the number of columns
        in the main data array.
        Provided description will be added to self.operations as:
        "XXXX Normalized Arrays".
        """
        self.norm=description
        self.array=self.array/(normFactors)
        self.operations.append('%s Normalized Arrays' % description)
       
        
    def arrayGroupCompare(self,grpA,grpB,minN=3, maxP=1):
        """given 2 groups of arrays (lists of integers),
        return a table of uid,A_mean,A_std,A-n,
        B_mean,B_std,B_n, t-statistic,2-tailed p-value
        Where A or B has fewer than N unmasked observarions,
        that oligo is not reported.
        """
        rv=[]
        
        subA=self.array[:,[self.arrayIdx[i] for i in grpA]]
        subB=self.array[:,[self.arrayIdx[i] for i in grpB]]

        for i in xrange(subA.shape[0]):
            a=subA[i,:].compressed()
            if len(a)<minN:
                continue
            b=subB[i,:].compressed()
            if len(b) <minN:
                continue
            t,p=stats.ttest_ind(a,b)
            if p> maxP:
                continue

            rv.append((self.oligoLongName[self.uidList[i]],
                       a.mean(),
                       a.std(),
                       len(a),
                       b.mean(),
                       b.std(),
                       len(b),
                       t,
                       p))
            rv.sort(key=lambda x: x[-1])
            
        return rv


    def cluster(self):
        """The corresponding MultiDictCluster 
        """
        rv=MultiDictCluster()
        rv.gNames=copy(self.oligoLongName)
        for i in xrange(len(self.uidList)):
            if not N.equal(self.array.mask,False)[i,:].any():
                continue
            rv[self.uidList[i]]={}
            row=self.array[i,:]
            for j in xrange(len(self.arrayList)):
                if row.mask[j]:
                    continue
                rv[self.uidList[i]][self.arrayList[j]]=row[j]
        return rv

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
            
