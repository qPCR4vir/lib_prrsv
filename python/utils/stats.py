"""Module for tracking and reporting statistcs from observations
"""

__verison__ = tuple([int(x) for x in
                     '$Revision: 1.3 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

import random
import numpy as N
import matplotlib as M
import pylab as P
import utils

def returnX(x):
    """returns the input argument
    
    Arguments:
    - `x`:
    """
    return x


class Tracker (utils.PickleableObject):
    """Observation tracker and statistic reporter
    """

    def __init__(self,callback=returnX,names=[]):
        """
        """
        self.callback=callback
        self.observations=[]        
        self.array=None
        self.names=names

    def append(self,items):
        """add observations
        
        Arguments:
        - `self`:
        - `items`: an element or iterable of elements
        """
        try:
            iter(items)
        except TypeError:
            items = (items,)

        self.observations.extend(
            (self.callback(i) for i in items))
        self.array=None


    def nArray(self):
        """return a numpy array of the observations
        
        Arguments:
        - `self`:
        """
        if self.array == None:
            self.array = N.array(self.observations)
        return self.array


    def dropCache(self):
        """remove the cached numpy array.  It is regenerated as needed.
        
        Arguments:
        - `self`:
        """
        self.array=None

    def mean(self,**kwds):
        """return the mean of the observations
        
        Arguments:
        - `self`:
        """
        return N.mean(self.nArray(),**kwds)

    def median(self,**kwds):
        """return the median of the observations
        
        Arguments:
        - `self`:
        """
        return N.median(self.nArray(),**kwds)
    
    def std(self,**kwds):
        """standard deviation of the observations
        
        Arguments:
        - `self`:
        """
        return N.std(self.nArray(),**kwds)

    def stderr(self,**kwds):
        """not implemented
        """
        pass

    def histogram(self,**kwds):
        """histogram of observations
        
        Arguments:
        - `self`:
        - `**kwds`: keywords arguments to numpy.histogram
        """
        return N.histogram(self.nArray(),**kwds)
        
    def xy(self,normed=False):
        """a histogram with integral bins (width=1)
        
        Arguments:
        - `self`:
        - `*kwds`:
        """
        return N.histogram(self.nArray(),normed=normed,
                           bins=range(int(N.floor(self.min())),
                                      int(N.ceil(self.max()+1))),
                           new=True)
    

    def max(self,**kwds):
        """returns max of observations
        
        Arguments:
        - `self`:
        """
        return self.nArray().max(**kwds)

    def min(self,**kwds):
        """returns min of observations
        
        Arguments:
        - `self`:
        """
        return self.nArray().min(**kwds)


    def boxplot(self,**kwds):
        """
        """
        P.boxplot(self.nArray(),**kwds)
    

    def scatter(self,xIdx=0,yIdx=1,**kwds):
        """
        
        Arguments:
        - `self`:
        - `xIdx`:
        - `yIdx`:
        - `**kwds`:
        """
        
        P.scatter(self.nArray()[:,xIdx],self.nArray()[:,yIdx])
        if xIdx < len(self.names):
            P.xlabel(self.names[xIdx])
        if yIdx < len(self.names):
            P.ylabel(self.names[yIdx])        
    
    
def randomData(low,high,size=1):
    """generator of random float tuples of length size.
    If size=1, the datum is returned as a single float.
    each datum is in the interval [low,high]
    """
    if size==1:
        while True:
            yield random.uniform(low,high)
    else:
        while True:
            yield tuple([random.uniform(low,high)
                         for x in xrange(size)])
        
