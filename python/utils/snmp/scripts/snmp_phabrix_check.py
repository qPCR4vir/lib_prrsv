#!/usr/local/bin/python -u
"""
 example of using low level snmp to the PHABRIX
 lock file status, using pstatus.py

 This requires a directive like:
 exec pstatus /var/phabrix/bin/pstatus.py ''

 plock files older that 1/2 hour produce warn states
 plock files older than 1 hour produce critical states
"""
__version__ = [int(x) for x in "$Revision: 1.1 $".split()[1].split('.')]
__author__ = "Kael Fischer <kael.fischer@gmail.com>"


import sys
import time,datetime

from utils import time2datetime,snmp,nagios

extTableOid = '.1.3.6.1.4.1.2021.8'
warnTime = datetime.timedelta(minutes=30)
deadTime = datetime.timedelta(minutes=60)

def query (ip):
    status = None
    message = ''

    try:
        extTable = snmp.walkTable(ip,extTableOid)
    except IOError ,eStr:
        return nagios.CRITICAL,'snmp error: '+ str(eStr)

    if len(extTable)==0:
        return nagios.CRITICAL,"snmp error: no extention table"
        
    pDict= dict(zip(extTable['extNames'],extTable['extOutput']))

    if 'pstatus' not in pDict:
        return nagios.CRITICAL,"snmp error: pstatus not in extention table"
        
    pOut = pDict['pstatus'].strip()
    
    cFields = pOut.split('Time:')
    if len(cFields) == 1:
        #no plock file, no plock history
        return nagios.CRITICAL,pOut

    try:
        timeStr = cFields[-1].strip()
        pTime = time2datetime(time.strptime(timeStr))
    except:
        return nagios.WARN, "(time parsing error) " + pOut

    if datetime.datetime.now() - pTime > deadTime:
        return nagios.CRITICAL,pOut
    elif datetime.datetime.now() - pTime > warnTime:
        return nagios.WARN,pOut
    else:
        return nagios.OK,pOut
    
    return (nagios.OK,extTable['extOutput'][0])


if __name__ == '__main__':
    ip = nagios.getopt(sys.argv[1:],__doc__)
    status,message = query(ip)
    nagios.print4nagios(status,message)
    sys.exit(status)
    

