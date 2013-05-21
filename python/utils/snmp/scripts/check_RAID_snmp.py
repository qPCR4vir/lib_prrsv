#!/usr/local/bin/python -u
"""
 Check status of RAID arrays using snmp

 Both the proware and areca controlers are handled

 On Proware these are checked:
 fan, power....

 On Acera these are checked:
 
 

"""
__version__ = [int(x) for x in "$Revision: 1.1 $".split()[1].split('.')]
__author__ = "Kael Fischer <kael.fischer@gmail.com>"


import sys
import time,datetime

from utils import time2datetime,snmp,nagios,iterCount

#
# oids
#
pwOID = {'serNum': '.1.3.6.1.4.1.14752.1.1.1.2',
         'event':'.1.3.6.1.4.1.14752.1.1.6',
         'defaultMsg':'.1.3.6.1.4.1.14752.1.1.1.1.0',
         'slotTable':'.1.3.6.1.4.1.14752.1.1.4'
         }

aceraOID= {'serNum': '.1.3.6.1.4.1.19736.1.1.1.3',
           'defaultMsg':'.1.3.6.1.4.1.19736.1.1.1.2.0',
           'hddStates':'.1.3.6.1.4.1.19736.1.1.3.1.1.9',
           'hddIdx': '.1.3.6.1.4.1.19736.1.1.3.1.1.1',
           'hwStatus': '.1.3.6.1.4.1.19736.1.1.2',
           'raidNames': '.1.3.6.1.4.1.19736.1.1.4.1.1.2',
           'raidStates': '.1.3.6.1.4.1.19736.1.1.4.1.1.4',
           'volNames':'.1.3.6.1.4.1.19736.1.1.5.1.1.2',
           'volStates':'1.3.6.1.4.1.19736.1.1.5.1.1.5',
           12000:'.1.3.6.1.4.1.19736.1.1.2.2.0',
           5000 :'.1.3.6.1.4.1.19736.1.1.2.3.0',
           3300 :'.1.3.6.1.4.1.19736.1.1.2.4.0',
           2500 :'.1.3.6.1.4.1.19736.1.1.2.5.0',
           1300 :'.1.3.6.1.4.1.19736.1.1.2.6.0',
           1250 :'.1.3.6.1.4.1.19736.1.1.2.7.0',
           'pwrTbl':'.1.3.6.1.4.1.19736.1.1.2.12.1.2',
           'fanTbl':'.1.3.6.1.4.1.19736.1.1.2.13.1.2',
           'dskTbl':'.1.3.6.1.4.1.19736.1.1.2.14',
           'boardTemp': '.1.3.6.1.4.1.19736.1.1.2.1.0'
           }


WARN_TEMP=40
WARN_MVOLT_PCT_TOL=0.03
WARN_LOW_FAN_RPM=2800
CRIT_LOW_FAN_RPM=500
ACREA_MV=(12000,5000,3300,2500,1300,1250)

def platform(ip):
    """return 'ACREA' or 'PROWARE' or 'UNRECOGNIZED' based on reply to
    snmp request for serial number.
    """
    if len(snmp.walk(ip,pwOID['serNum'])) >0:
        return 'PROWARE'
    elif len(snmp.walk(ip,aceraOID['serNum'])) >0:
        return 'ACREA'
    else:
        return 'UNRECOGNIZED'


def pwQuery(ip):
    """run the proware query
    """

    status = nagios.OK
    message = snmp.get(ip,pwOID['defaultMsg']).strip()
    errMessages = []
    
    # check proware events
    events = snmp.walk(ip,pwOID['event'])
    for name,sStat in events:
        if sStat != 'Good(1)':
            status=nagios.CRITICAL
            name = name.rstrip('.0').replace('event','')
            errMessages.append('%s-%s'%(name,sStat))

    # disk status
    sTable = snmp.walkTable('192.168.0.112',pwOID['slotTable'])
    dStatCount = iterCount(sTable['slotStatus'])

    #hot spare check
    if 'Free(2)' not in dStatCount:
        dStatCount['Free(2)'] = 0
        if dStatCount['Free(2)'] != 1:
            status = nagios.CRITICAL
            errMessages.append('%s Hot Spares'%(dStatCount['Free(2)']))

    # indivdual slot check
    slotStat = dict(zip(sTable['slotDiskIndex'],zip(sTable['slotStatus'],sTable['slotBadBlockNumber'])))
    sIdx = slotStat.keys()
    sIdx.sort()

    for s in sIdx:
        sStat,sBB = slotStat[s]
        if sStat not in ('Arranged(1)','Free(2)'):
            status=nagios.CRITICAL
            errMessages.append('Slot %s Status-%s'%(s,sStat))
        if sBB > 0:
            errMessages.append('Slot %s-%s Bad Blocks'%(s,sBB))
            if status == nagios.OK:
                status=nagios.WARN
    
    if len(errMessages) > 0:
        message = ';'.join(errMessages)
    return status,message


def aceraQuery(ip):
    """run the proware query
    """
    status = nagios.OK
    message = snmp.get(ip,aceraOID['defaultMsg']).strip()
    errMessages = []

    # check voltages
    for mv in ACREA_MV:
        sMV=snmp.get(ip,aceraOID[mv])
        if float(abs(sMV-mv))/float(mv) >= WARN_MVOLT_PCT_TOL:
            status=nagios.WARN
            errMessages.append('%s mV line is %s mV'%(mv,sMV))
    # temperatures
    sBT=snmp.get(ip,aceraOID['boardTemp'])
    if sBT >= WARN_TEMP:
            status=nagios.WARN
            errMessages.append('Controller Temp-%s C'%(sBT))
    dskTbl=snmp.walkTable(ip,aceraOID['dskTbl'])
    dskDict=dict(zip(dskTbl['hwHddIndex'],dskTbl['hwHddTemp']))
    for dsk,tmp in dskDict.items():
        if tmp >= WARN_TEMP:
            status=nagios.WARN
            errMessages.append('Disk %s Temp-%s C'%(dsk,tmp)) 
    # fan RPM
    minRPM=min(snmp.walkList(ip,aceraOID['fanTbl']))
    if minRPM <= WARN_LOW_FAN_RPM:
            status=nagios.WARN
            errMessages.append('Fan RPM Low-%s'%(minRPM))
    # CRITICAL tests follow
    if minRPM <= CRIT_LOW_FAN_RPM:
            status=nagios.CRITICAL
            errMessages.append('Fan RPM Low-%s'%(minRPM))
    # power supplies
    pwrCt=snmp.countList(ip,aceraOID['pwrTbl'])
    for stat,n in pwrCt.items():
        if stat != 'Ok(1)':
            status=nagios.CRITICAL
            errMessages.append('Power Supply Status-%s'%(stat))
    # disk roles
    dskStates = snmp.walkList(ip,aceraOID['hddStates'])
    dskIdx=snmp.walkList(ip,aceraOID['hddIdx'])

    dsCt = iterCount(dskStates)
    if 'Hot Spare' not in dsCt:
        dsCt['Hot Spare'] = 0
        if dsCt['Hot Spare'] != 1:
            status = nagios.CRITICAL
            errMessages.append('%s Hot Spares'%(dsCt['Hot Spare']))
    for n,s in enumerate(dskStates):
        if s not in ('RaidSet Member','Hot Spare'):
            status = nagios.CRITICAL
            errMessages.append('Slot %s Status-%s'%(dskIdx[n],s))
    # raid+volume status
    rNames = snmp.walkList(ip,aceraOID['raidNames'])+snmp.walkList(ip,aceraOID['volNames'])
    rStates = snmp.walkList(ip,aceraOID['raidStates'])+snmp.walkList(ip,aceraOID['volStates'])
    for n,s in enumerate(rStates):
        if s != "Normal":
            status = nagios.CRITICAL
            errMessages.append('%s Status-%s'%(rNames[n],s))
    
    if len(errMessages) > 0:
        message = ';'.join(errMessages)
    return status,message




if __name__ == '__main__':
    ip = nagios.getopt(sys.argv[1:],__doc__)
    try:
        p = platform(ip)
    except IOError,e:
        status = nagios.CRITICAL
        message='snmp error: '+ str(e)
    else:
        if p == 'PROWARE':
            status,message = pwQuery(ip)
        elif  p == 'ACREA':
            status,message = aceraQuery(ip)
        else:
            status = nagios.CRITICAL
            message="unknown hardware at %s" %ip

    nagios.print4nagios(status,message)
    sys.exit(status)
    

