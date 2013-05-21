# snmp utilities
# and class that provides simple nagios
# plugin interface.
#
__version__ = [int(x) for x in "$Revision: 1.2 $".split()[1].split('.')]
__author__ = "Kael Fischer <kael.fischer@gmail.com>"

import datetime
import re
from utils import iterCount, dict2optStr
from commands import getstatusoutput

WALK_PATH = "/usr/local/bin/snmpwalk"
GET_PATH = "/usr/local/bin/snmpget"

#see snmpcmd(1)
SNMP_OPTIONS = {
    #'-O':'q',   # output options
    '-v':'1',   # snmp version
    '-m':'ALL'  # gobble all MIBs
    }

COMMUNITY = "RORO"

def oidSplit(oid):
    """return (head,tail) where tail is everything after the last '.'
    and head is every thing before it. tail may be empty string.
    """
    if oid.count('.') == 0:
        return (oid,'')
    else:
        oSplit = oid.split('.')
        return ('.'.join(oSplit[:-1]),oSplit[-1])

    

def snmpConvert(snmpResult,trimMIB=False):
    """Given snmp data in the form: TYPE FORMATTED_STUFF,
    returns a converted value as follows:

    snmp type      python type
     INTEGERx         int  (note some integer data are like 
                            state(idx) - these are returned as strings
                            e.g.: direct(3),up(1),softwareLoopback(24),
                            etc..)
     COUNTERx         int
     GAUGEx           (int, unit string)     
     IpAddress        str
     Network Address  str
     OID              str
     STRING           str
     Timeticks        timedelta
     <all others>     str

     If a '=' is present in the input it is assumed that the input is a
     complete line for the output of get or walk.  In that case a tuple is
     returned like (oid, converted_value).  Otherwise only the converted
     value is returned. That could be a problem if there is an '=' in the
     value and the line in not a complete line.  So just sent the complete
     line, OK!
    """

    oid = None

    if snmpResult.count(' = ') > 0:
        oid,snmpResult = snmpResult.split(' = ',1)
        if trimMIB:
            oid = oid.split('::')[-1]
    try:
        sSplit = snmpResult.strip().split(': ',1)
        if len(sSplit) > 1:
            sUnit,sValue = snmpResult.strip().split(': ',1)
        else:
            sUnit=snmpResult.strip(snmpResult.strip(": "))
            sValue =''
            
        if sUnit.startswith('INTEGER') or sUnit.startswith('COUNTER'):
            try:
                value=int(sValue.strip())
            except:
                value = sValue
        elif sUnit.startswith("GAUGE"):
            vSplit = sValue.split(' ',1)
            if len(vSplit) == 2:
                value = (int(vSplit[0].strip()),vSplit[1])
            else:
                value = (int(vSplit[0]),'')
        elif sUnit.startswith('Timeticks'):
            ticks = re.search('\((\d+)\)',sValue).group(1)
            value = datetime.timedelta(milliseconds=10*int(ticks))

        else:
            value=sValue.strip("\'\"")

        if oid == None: 
            return value
        else:
            return (oid,value) 
    except:
        print  oid,snmpResult
        raise

def walk (ip,oid,trimMIB=True):
    """snmp walk.  Returns a list of (oid,stringValues).
    By default, returned oids have <MIB>:: trimmed off.
    """
    s, o = getstatusoutput("%s %s -c %s %s %s " % (
        WALK_PATH,dict2optStr(SNMP_OPTIONS),
        COMMUNITY,ip,oid))
    if s != 0:
        raise IOError , "query of %s at %s failed" % (oid,ip)
    else:
        return [snmpConvert(l,trimMIB=trimMIB) for l in o.split('\n') if l.count(' = ') >0]
    

def get (ip,oid):
    """snmp get.  Returns a the value requested as a string.
    """
    s, o = getstatusoutput("%s %s -c %s %s %s " % (
        GET_PATH,dict2optStr(SNMP_OPTIONS),
        COMMUNITY,ip,oid))
    if s != 0:
        raise IOError , "query of %s at %s failed" % (oid,ip)
    else:
        return snmpConvert(o.strip().split(' = ')[1])


def walkList(ip,oid):
    """returns a list of string values 
    oid must be evaluate to an instance list.
    e.g. disk.0, disk.1, etc...
    """
    tups = walk(ip,oid)
    rootOid,startIdx = oidSplit(tups[0][0])
    startIdx =  int(startIdx)

    for n,tup in enumerate (tups):
        oid, val = tup
        oHead,oTail = oidSplit(oid)
        myRoot = oHead 
        myN = int(oTail)
        if myRoot != rootOid:
            raise ValueError , "returned OIDs are not a single list e.g (%s , %s)" %(tups[0][0],oid)
        elif myN != n+startIdx:
            raise ValueError , "returned oids not a consecutive list: %s" % oid

    return [val for oid,val in tups]
        
    
def countList(ip,oid):
    """Go through a list (specified by oid) and
    return a dict of { value:count,...} for allt values in list.
    """
    return iterCount(walkList(ip,oid))


def walkTable(ip,tableOid,trimMIB=True):
    """return a dict of OID: [walkList(OID)], walking under tableOid.
    Returned oids are trimmed to remove 'MIB::' prefix, if any. 
    """
    rv = {}
    tups = walk(ip,tableOid,trimMIB=False)

    for o,v in tups:
        oHead = oidSplit(o)[0]
        if oHead not in rv:
            rv[oHead] = []

    for oid in rv.keys():
        rv[oid]=walkList(ip,oid)
        if trimMIB:
            trimOid = oid.split('::')[-1]
            rv[trimOid] = rv[oid]
            del rv[oid]

    return rv

        
    
