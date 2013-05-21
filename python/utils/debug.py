#!/usr/local/bin/python 
#
# utils/debug.py
#
# Debugging helper functions
#

__verison__ = tuple([int(x) for x in
                     '$Revision: 1.1 $'.split()[1].split('.')])
__author__ = "Kael Fischer"

from __init__ import *
import types

def dictTypes (thing):
    """Return a list of things __dict__ keys and
    type(value)
    """
    return [ (k,type(v)) for k,v in thing.__dict__.items()]

def itemsOfType (inDict,matchType):
    """return a generator of key, value tuples from
    dictionary where the values are in matchType
    (which may be a sequence of types or objects).
    """
    myTypes = []
    for t in getIterable(matchType):
        if type(t) != types.TypeType:
            tt=type(t)
        else:
            tt = t
        if t not in myTypes:
            myTypes.append(tt)

    for k,v in inDict.items():
        if type(v) in myTypes:
            yield (k,v)


def moduleWalk(top,history=None):
    """go through module (top) name space and recurse through
    modules in namespace.  returns a generator that yields
    each module found once.
    """
    if history is None:
        history=[]
    else:
        # note that using and passing a list for history
        # allows recursive calls to know about prior history
        # and the caller to know what the children did  
        history=history

    myModules=[]

    if type(top) == types.DictType:
        myModules=[v for k,v in itemsOfType(top,types.ModuleType)
                   if v not in history]
    elif type(top) == types.ModuleType:
        myModules=[top]
    else:
        myModules=list(getIterable(myModules))

    while len(myModules) > 0:
        m=myModules.pop()
        if type(m) == types.ModuleType:
            if m not in history:
                history.append(m)
                yield m
                for sm in (x for x in moduleWalk(m.__dict__,history)):
                    yield sm

    
def safeVersion(mod):
    """if mod.__version__ is defined return it
    else return None
    """
    try:
        return mod.__version__
    except AttributeError:
        return None
    
def versionWalk(top):
    """returns generator of module,module.__version__
    when it is defined.  Walks module tree in top module
    namespace.  top can be a sequence (e.g. globals) of
    objects, non-modules in the sequence are ignored.
    """
    for m in moduleWalk(top):
        v=safeVersion(m)
        if v is not None:
            yield (m,v)
            
