#!/usr/bin/env python
"""
Primer3Wrapper.py
* wrapper utility to run primer3, and catch std i/o/error
* real-time PCR primer design strategy incorporated
"""

import time
from subprocess import *

run_time = time.asctime().replace(" ", "_").replace(":", "")

class InputDict:
    """class InputDict takes initial parameters and those from GUI input"""
    def __init__(self, gui_input_dict):
        self.input_dict={"PRIMER_PICK_ANYWAY" : 1,
                      "PRIMER_PRODUCT_MAX_TM" : 90,
                      "PRIMER_PRODUCT_MIN_TM" : 75,
                      "PRIMER_PRODUCT_SIZE_RANGE" : "90-150 80-300",
                      "PRIMER_MIN_SIZE" : 18,
                      "PRIMER_MAX_SIZE" : 22,
                      "PRIMER_MIN_TM" : 59,
                      "PRIMER_MAX_TM" : 61,
                      "PRIMER_MIN_GC" : 45,
                      "PRIMER_MAX_GC" : 55,
                      "PRIMER_SELF_ANY" : 3,
                      "PRIMER_SELF_END" : 0,
                      "PRIMER_NUM_RETURN" : 3}
        self.input_dict.update(gui_input_dict)

class WrapperOutput:
    """data interface as ThinWrapper output.
    output will go to result file or/and log file"""
    def __init__(self):
        self.log = ''
        self.result = ''
    def receive(self, l='', r=''):
        self.log += l
        self.result = r
    def warn(self, WarnMessage):
        self.result = WarnMessage + '\n\n' + self.result
    def get_results(self):
        return self.log, self.result

#------------------------------------------------------------------------------------------
class ThinWrapper:
    """
    ThinWrapper class takes care of running primer3_core,
       feeding input and collect output; channel to a WrapperOutput instance
       """
    def __init__(self):
        self.TmpInputFile = "tmpPrimer3coreInputFile"
        self.FailureMessage = 'PrimerPy failed to get acceptable primers from your sequence :('
        self.WarningMessage = "@@@ No acceptable primers were found for real-time PCR,      @@@\n@@@ the provided primers may not work, use at your own risk. @@@"
    
    def boulder_write(self, dict):
        "write text stream in BoulderIO format from input dictionary"
        b=""
        for key, value in dict.items():
            b += key+'='+str(value)+'\n'
        b+='=\n' #trailing "=" in BoulderIO file
        return b

    def write_input(self, input_str):
        myhandle=open(self.TmpInputFile, 'w')
        myhandle.write(input_str)
        myhandle.close()

    #-----------------------------------------------------------------
    #    core function of using primer3
    #    read in GUI parameters, update default parameters
    #-----------------------------------------------------------------
    def call_primer3_core(self, input_dict, output_object, pcmd, run_NUM):
        # execution unit, called in design_logic_loop
        # output_object is an instance of WrapperOutput, to store info to log, result
        # create input temporary file
        input_str = self.boulder_write(input_dict)
        self.write_input(input_str)
        #pcmd: exe command str, from design_logic_loop
        p=Popen(pcmd, shell=True, stdin=PIPE, stdout=PIPE,
                stderr=PIPE)#, close_fds=True)
        (primer3_in, primer3_out, primer3_err) = (p.stdin, p.stdout, p.stderr)
        #older version was (primer3_in, primer3_out, primer3_err)=os.popen3(pcmd)
        primer3_out_text = primer3_out.read()
        if primer3_out_text.rstrip() == "":
            output_object.receive(l='Primer3 runtime error occured', r=primer3_err.read())
            return 1, output_object
        elif "NO PRIMERS FOUND" not in primer3_out_text[:500]:
            # primer3 run success
            outheader = run_time + '_run' + str(run_NUM) + '\n\n'
            if input_dict.has_key('Users_Note'):
                outheader = input_dict['Users_Note'] + '\n' + outheader
            output_object.receive(l='Done.\n', r=outheader+primer3_out_text)
            return 1, output_object
        else:
            output_object.receive(l='Revised...')
            return 0, output_object
        
    #-------------------------------------------------------------------
    #   Design function 
    #
    #-------------------------------------------------------------------
    def design_logic_loop(self, Primer3_exe_Path, input_dict, run_NUM=1):
        # input_dict contain updated parameters (both default and gui)
        pcmd = Primer3_exe_Path + " < " + self.TmpInputFile
        run_output = WrapperOutput()
        run_output.receive(l='\n'+input_dict["PRIMER_SEQUENCE_ID"]+'\n')
        result, run_output = self.call_primer3_core(input_dict, run_output, pcmd, run_NUM)
        if result != 1:
            input_dict["PRIMER_SELF_ANY"]=4
            result, run_output = self.call_primer3_core(input_dict, run_output, pcmd, run_NUM)
            if result != 1:
                input_dict["PRIMER_MIN_TM"]=58
                input_dict["PRIMER_MAX_TM"]=62
                result, run_output = self.call_primer3_core(input_dict, run_output, pcmd, run_NUM)
                if result != 1:
                    input_dict["PRIMER_MIN_GC"]=40
                    input_dict["PRIMER_MAX_GC"]=60
                    result, run_output = self.call_primer3_core(input_dict, run_output, pcmd, run_NUM)
                    if result != 1:
                        input_dict["PRIMER_MIN_GC"]=30
                        input_dict["PRIMER_MAX_GC"]=70
                        result, run_output = self.call_primer3_core(input_dict, run_output, pcmd, run_NUM)
                        if result != 1:
                            input_dict["PRIMER_SELF_END"]=1
                            result, run_output = self.call_primer3_core(input_dict, run_output, pcmd, run_NUM)
                            if result != 1:
                                run_output = self.noluck(input_dict, run_output, pcmd, run_NUM)
        return run_output

    def noluck(self, input_dict, output_object, pcmd, run_NUM):
        # design at loose condition, will give warning
        loose_dict={"PRIMER_PICK_ANYWAY" : 1,
                      "PRIMER_PRODUCT_MAX_TM" : 93,
                      "PRIMER_PRODUCT_MIN_TM" : 75,
                      "PRIMER_PRODUCT_SIZE_RANGE" : "80-300 80-700",
                      "PRIMER_MIN_SIZE" : 17,
                      "PRIMER_MAX_SIZE" : 25,
                      "PRIMER_MIN_TM" : 56,
                      "PRIMER_MAX_TM" : 64,
                      "PRIMER_MIN_GC" : 30,
                      "PRIMER_MAX_GC" : 70,
                      "PRIMER_SELF_ANY" : 8,
                      "PRIMER_SELF_END" : 8,
                      "PRIMER_NUM_RETURN" : 3}
        input_dict.update(loose_dict)
        result, output_object = self.call_primer3_core(input_dict, output_object, pcmd, run_NUM)
        if result != 1:
            output_object.receive(l='Failed.', r=self.FailureMessage)
        else:
            output_object.warn(self.WarningMessage)
        return output_object



if __name__ == '__main__':
    #obsolete, 9-30-07
    from GlobalConf import *
    global_config = GlobalConf().readconf()
    Primer3_exe_Path = global_config['Primer3_exe_Path']
    test_sequence='CCAGTAGCATATGCTTGTCTCAAAGATTAAGCCATGCAAGTGTAAGTACACACGGGCTGTACAGTGAAACTGCGAATGGCTCATTAAATCAGTTATGGTTCCTTTGATCGCTCCACCGTTACTTGGATAACTGTGGCAATTCTAGAGCTAATACATGCAAACGAGCGCTGACCCCCGGGGATGCGTGCATTTATCAGATCCAAAACCCATGCGGGACGGGCCCTTCCGGGGGCCCGCCCCGGCCGCTTTGGTG'
    test_in = GuiInputs(template=test_sequence)
    test_dict = InputDict(test_in).returndict()
    test = ThinWrapper()
    test_result = test.design_logic_loop(Primer3_exe_Path, test_dict, run_NUM=13)
    print test_result.log, test_result.result
    
