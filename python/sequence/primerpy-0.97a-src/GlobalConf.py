#!/usr/bin/env python

"""GlobalConf.py
    Read or Write (from GUI input) Global configuration file for OligoBench/PrimerPy
    """

import sys

#for POSIX systems
conf_file = '_primerpy_posix.conf'
default_conf = 'UserDataDir=result\nPrimer3_exe_Path=primer3/primer3_core -format_output\nLogFile=primerpy.log'
if 'win32' in sys.platform:
    conf_file = '_primerpy_win32.conf'
    default_conf = 'UserDataDir=result\nPrimer3_exe_Path=primer3\\primer3_core -format_output\nLogFile=primerpy.log'

class GlobalConf:
    def __init__(self):
        #conf check is carried out in primerpy/MyMainFrame
        self.ConfDict = self.readconf()
        
    def readconf(self, input_file=conf_file):
        # read from conf.txt
        confdict = {}
        w = open(input_file).readlines()
        for line in w:
            try:
                a = line.rstrip().split('=')
                confdict[a[0]] = a[1]
            except (IndexError, KeyError):
                pass
        return confdict
        
    def writeconf(self, user_defined, write_file=conf_file):
        #write over primerpyconf by user defined parameters from GUI
        w = open(write_file, 'w')
        s = self.dict_write(user_defined)
        w.write(s)
        w.close()
        
    def reset(self, write_file=conf_file):
        #reset conf.txt to default
        w = open(write_file, 'w')
        w.write(default_conf)
        w.close()
        
    def dict_write(self, dict):
        b=""
        for key, value in dict.items():
            b += key+'='+str(value)+'\n'
        return b
        
    #----------------------------
    #handling record file, where run_num and file name of last_run are kept
    def get_record(self, recordfile='_primerpy.record'):
        #get run_num, last_run_result from record file
        #records for run_num, last_run_result
        w = open(recordfile).readlines()
        record = {'run_num': int(w[1].rstrip()),
                    'last_run_result': w[2].rstrip()}
        return record
    
    def flush_record(self, record, recordfile='_primerpy.record'):
        #records for run_num, last_run_result
        w = open(recordfile, 'w')
        w.write('#records for run_num, last_run_result\n' + str(record['run_num'])
                +'\n'+ record['last_run_result']+'\n')
        w.close()
 
