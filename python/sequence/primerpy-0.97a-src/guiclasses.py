# coding=utf-8
"""guiclasses.py,
classes for PrimerPy graphic interface


"""

import wx, os
import wx.html
import wx.lib.flatnotebook as fnb
from Primer3Wrapper import *
from dimercheck import *

from GlobalConf import *
#get global parameters and run records
GC = GlobalConf()
record = GC.get_record()
#run_num, last_run_result
RUN_NUM = record['run_num']
LAST_RUN_RESULT = record['last_run_result']
#path & files
global_config = GC.ConfDict
Primer3_PATH = global_config['Primer3_exe_Path']
UserDataDir = global_config['UserDataDir']
LogFile = os.path.join(UserDataDir, global_config['LogFile'])
#   Check configurations in primerpy/MyMainFrame.
#   If primer3 or usr path incorrect, launch configure_dialog

class TemplatePanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        vbox = wx.BoxSizer(wx.VERTICAL)
        text1 = wx.StaticText(self, label='Paste your template sequence here:')
        vbox.Add(text1)
        self.template = wx.TextCtrl(self, size=(270,90), style=wx.TE_MULTILINE)
        self.template.SetInsertionPoint(0)
        vbox.Add(self.template)
        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.paste_button = wx.Button(self, 11, 'Paste')
        self.Bind(wx.EVT_BUTTON, self.paste, self.paste_button)
        hbox1.Add(self.paste_button)
        self.clear_button = wx.Button(self, 12, 'Clear All')
        self.Bind(wx.EVT_BUTTON, self.clear_all, self.clear_button)
        hbox1.Add(self.clear_button)
        hbox1.Add((15,-1))
        self.design_button = wx.Button(self, 13, 'Design')
        self.Bind(wx.EVT_BUTTON, self.design, self.design_button)
        hbox1.Add(self.design_button)
        vbox.Add(hbox1)
        self.SetSizer(vbox)
        
    def ProcessTemplate(self, text):
        "process input template, return PRIMER_SEQUENCE_ID, SEQUENCE"
        global run_time, RUN_NUM
        if text[0]==">":
            firstline=max(text.find('\n'), text.find('\r'))
            t = text[firstline+1: ].replace("\n", "").replace("\r", "")
            return text[ :firstline], t
        else:
            #this is executed before RUN_NUM+=1
            return run_time+'_run'+str(RUN_NUM+1), text.replace("\n", "").replace("\r", "")
    def get(self):
        "check for empty template, return seq dict"
        text = self.template.GetValue()
        text=text.replace(" ", "")
        if text == "":
            raise Exception('No template')
        else:
            a = self.ProcessTemplate(text)
            return {'PRIMER_SEQUENCE_ID' : a[0], 'SEQUENCE' : a[1]}
    def paste(self, event):
        self.template.Paste()
    def clear_all(self, event):
        self.template.Clear()
        self.GetParent().clear2()
    def design(self, event):
        self.GetParent().design()

class OptionPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        sizer = wx.GridSizer(rows=4, cols=2)
        sizer.Add(wx.StaticText(self, label="Predefine left primer:"))
        sizer.Add(wx.StaticText(self, label="Users' Note:"))
        self.LeftPrimer = wx.TextCtrl(self, size=(160,-1))
        sizer.Add(self.LeftPrimer)
        self.UserNote = wx.TextCtrl(self, size=(160,-1))
        sizer.Add(self.UserNote)
        sizer.Add(wx.StaticText(self, label="Predefine right primer:"))
        sizer.Add(wx.StaticText(self, label="Desired Amplicon size:"))
        self.RightPrimer = wx.TextCtrl(self, size=(160,-1))
        sizer.Add(self.RightPrimer)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.AmpliconSize = wx.TextCtrl(self)
        hbox.Add(self.AmpliconSize)
        hbox.Add(wx.StaticText(self, label="(80~300)"))
        sizer.Add(hbox)
        self.SetSizer(sizer)
        self.Show(False)
    def get(self):
        "return a dictionary of paramters"
        options_dict={
            'PRIMER_LEFT_INPUT' : self.LeftPrimer.GetValue().replace(" ", ""),
            'PRIMER_RIGHT_INPUT' : self.RightPrimer.GetValue().replace(" ", ""),
            'PRIMER_PRODUCT_OPT_SIZE' : self.AmpliconSize.GetValue().replace(" ", ""),
            'Users_Note' : self.UserNote.GetValue()
            }
        return options_dict
    def clear3(self):
        self.LeftPrimer.Clear()
        self.RightPrimer.Clear()
        self.UserNote.Clear()
        self.AmpliconSize.Clear()


class InputPage(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        icon = wx.StaticBitmap(self, -1, wx.Bitmap('eleboa.jpg'), pos=(30,60), size=(120,80))
        #-------------Template Panel
        self.tp = TemplatePanel(self)
        self.tp.SetPosition((180, 30))
        self.tp.SetSize((300, 150))
        #-------------
        self.more_button = wx.ToggleButton(self, -1, 'More Options >>>', pos=(180,180))
        #-------------Option Panel
        self.mo = OptionPanel(self)
        self.mo.SetPosition((160, 220))
        self.mo.SetSize((360, 88))
        self.Bind(wx.EVT_TOGGLEBUTTON, self.ShowOptions, self.more_button)
        
        self.Center()
        self.Show()

    def ShowOptions(self, event):
        self.mo.Show(self.more_button.GetValue())
    def clear2(self):
        # continued from self.tp.clear_all
        if self.more_button.GetValue():
            self.mo.clear3()
    def GetValues(self):
        gui_input_dict = self.tp.get()
        if self.more_button.GetValue():
            x = self.mo.get()
            for akey in x:
                if x[akey] != '':
                    gui_input_dict[akey] = x[akey]
        return gui_input_dict

    def design(self):
        global GC, RUN_NUM, UserDataDir, Primer3_PATH, LAST_RUN_RESULT, run_time
        try:
            gui_input_dict = self.GetValues()
            #except "No template"
            RUN_NUM += 1
            LAST_RUN_RESULT = run_time + 'run' + str(RUN_NUM) + '.html'
            LAST_RUN_RESULT = os.path.join(UserDataDir, LAST_RUN_RESULT)
            #
            input_dict = InputDict(gui_input_dict).input_dict
            run = ThinWrapper()
            run_result = run.design_logic_loop(Primer3_PATH, input_dict, RUN_NUM)
            self.write_result_log(run_result, LAST_RUN_RESULT)
            record = {'run_num' : RUN_NUM, 'last_run_result' : LAST_RUN_RESULT}
            GC.flush_record(record)
            self.GetParent().show_result()
            self.GetParent().GetParent().statusbar.SetStatusText("Done...")
        except:
            self.GetParent().GetParent().statusbar.SetStatusText("No Template!")

    def write_result_log(self, run_result, LAST_RUN_RESULT):
        global LogFile
        de = SequenceUtil()
        out = open(LAST_RUN_RESULT, 'w')
        #reformat the output
        html_text = de.reformat(run_result.result)
        out.write(html_text)
        out.close()
        out = open(LogFile, 'a')
        out.write(run_result.log)
        out.write('result saved as: ' +LAST_RUN_RESULT+ '\n')
        out.close()

class CheckDimerPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        wx.StaticText(self, label='Evaluate the quality of any primers:', pos=(50, 20))
        wx.StaticText(self, label='Sense primer(s)', pos=(50, 50))
        self.sense = wx.TextCtrl(self, size=(200,50), pos=(50, 70),  style=wx.TE_MULTILINE)
        wx.StaticText(self, label='Antisense primer(s)', pos=(270, 50))
        self.antisense = wx.TextCtrl(self, size=(200,50), pos=(270, 70),  style=wx.TE_MULTILINE)
        self.dimer_button = wx.Button(self, -1, 'Start Evaluation', pos=(270, 140))
        self.Bind(wx.EVT_BUTTON, self.dimer_evaluate, self.dimer_button)
        wx.StaticText(self, label="Notes:\n"
                                  "one primer per line, any number of primers.\n"
                                  "degenerate primers are interpreted into all possibilities.\n"
                                  "non-base characters are skipped.\n",
                                    pos=(50, 250))

    def dimer_evaluate(self, event):
        #get primers
        sense_list = list(set(self.sense.GetValue().splitlines()))
        anti_list = list(set(self.antisense.GetValue().splitlines()))
        sl, al, batch = [], [], []
        de = SequenceUtil()
        for item in sense_list:
            sl += de.deambiguity(item)
        for item in anti_list:
            al += de.deambiguity(item)
        if len(sl) > 0 and len(al) > 0:
            for x in sl:
                for y in al:
                    batch.append((x,y))
            html_text = de.BatchCheck(batch)
        elif len(sl) > 0:
            html_text = de.BatchCheck_single(sl)
        elif len(al) > 0:
            html_text = de.BatchCheck_single(al)
        else:
            html_text = ''
        if html_text != '':
            batchcheck_result = os.path.join(UserDataDir, 'DIMER_CHECK_RESULTS.html')
            de.html_write(html_text, batchcheck_result)
            self.GetParent().show_dimers(batchcheck_result)


class PrimerNoteBook(fnb.FlatNotebook):
    def __init__(self, parent):
        fnb.FlatNotebook.__init__(self, parent, -1, style=fnb.FNB_VC8|fnb.FNB_NO_NAV_BUTTONS)
        #self.page_changed = fnb.EVT_FLATNOTEBOOK_PAGE_CHANGED
        self.show_design()

    def show_design(self):
        #if not hasattr(self, 'sheet_design'):
        self.sheet_design = InputPage(self)
        self.AddPage(self.sheet_design, "design")
        self.sheet_design.SetFocus()
    def show_dimercheck(self):
        self.sheet_dimercheck = CheckDimerPanel(self)
        self.AddPage(self.sheet_dimercheck, "Check dimers")

    def show_result(self):
        global LAST_RUN_RESULT
        self.sheet_result = wx.html.HtmlWindow(self)
        self.AddPage(self.sheet_result, "result")
        self.sheet_result.LoadFile(LAST_RUN_RESULT)
        #self.sheet_result.SelectionToText()
    def show_dimers(self, f):
        self.sheet_dimers = wx.html.HtmlWindow(self)
        self.AddPage(self.sheet_dimers, "dimers check")
        self.sheet_dimers.LoadFile(f)
    def show_log(self):
        global LogFile
        #if not hasattr(self, 'sheet_log'):
        self.sheet_log = wx.html.HtmlWindow(self)
        self.AddPage(self.sheet_log, "log")
        self.sheet_log.LoadFile(LogFile)
    def show_manual(self):
        #if not hasattr(self, 'sheet_manual'):
        self.sheet_manual = wx.html.HtmlWindow(self)
        self.AddPage(self.sheet_manual, "manual")
        self.sheet_manual.LoadFile('primerpy_tutorial.htm')


class configure_dialog(wx.Dialog):
    def __init__(self, parent):
        wx.Dialog.__init__(self, parent, -1, "Configure PrimerPy", size=(300, 300))
        wx.StaticBox(self, -1, 'Only modify the defaults if you have to', (5, 5), size=(290, 235))
        panel = wx.Panel(self)
        topgrid = wx.GridSizer(rows=3, cols=2)
        topgrid.Add(wx.StaticText(panel, -1, "Directory of results:"))
        self.datadir = wx.TextCtrl(panel, -1, "result", size=(120,-1))
        topgrid.Add(self.datadir)
        topgrid.Add(wx.StaticText(panel, -1, "Directory to primer3:"))
        self.primer3dir = wx.TextCtrl(panel, -1, "primer3", size=(120,-1))
        topgrid.Add(self.primer3dir)
        topgrid.Add(wx.StaticText(panel, -1, "Name your log file:"))
        self.logname = wx.TextCtrl(panel, -1, "primerpy.log", size=(120,-1))
        topgrid.Add(self.logname)
        panel.SetSizer(topgrid)
        panel.SetPosition((10, 35))
        panel.SetSize((280, 90))
        wx.StaticText(self, -1,
                                "Explanations:\n"
                                "  PrimerPy stores each run as a simple\n  html file in the Directory of results,\n  where the log file also resides.\n  Directory to primer3 is where you\n  find primer3_core.\n",
                                pos = (10, 140),
                                style=wx.ALIGN_LEFT)
        
        self.restore_button = wx.Button(self, -1, "Restore to Default", (15, 245))
        self.Bind(wx.EVT_BUTTON, self.reset, self.restore_button)
        self.conf_button = wx.Button(self, -1, "Save Configuration", (160, 245))
        self.Bind(wx.EVT_BUTTON, self.configure, self.conf_button)

    def input_dict(self):
        #collect user inputs, return a dict
        mydict = {}
        mydict['UserDataDir'] = self.datadir.GetValue().rstrip()
        mydict['Primer3_exe_Path'] = os.path.join(self.primer3dir.GetValue().rstrip(),
                                                  'primer3_core') + ' -format_output'
        mydict['LogFile'] = self.logname.GetValue().rstrip()
        return mydict

    def reset(self, event):
        global GC
        GC.reset()
        self.Close(True)
        self.popup()
    def configure(self, event):
        input_dict = self.input_dict()
        if not os.path.isdir(input_dict['UserDataDir']):
            os.mkdir(input_dict['UserDataDir'])
        global GC
        GC.writeconf(input_dict)
        self.Close(True)
        self.popup()
    def popup(self):
        p = wx.MessageDialog(self, 'Please restart PrimerPy\nto use the changed configuration', 'need restart', wx.OK|wx.ICON_INFORMATION)
        p.ShowModal()
        p.Destroy()


class AboutPrimerPy(wx.Dialog):
    def __init__(self, parent):
        wx.Dialog.__init__(self, parent, -1, "About PrimerPy", size=(300, 300) )
        icon = wx.StaticBitmap(self, -1, wx.Bitmap('eleboa.jpg') )
        info = wx.StaticText(self, -1,
            "PrimerPy is a primer design tool\n"
            "optimized for SyBr Green QPCR, \nfree under GNU/GPL.\n"
            "Copyright 2007, Shuzhao Li\n"
            "http://www.bioinformatics.org/primerpy\n\n"
            "This version of PrimerPy uses Primer3 (primer3.sourceforge.net) as backend engine.\n", 
            style=wx.ALIGN_CENTER
            )
        button = wx.Button(self, wx.ID_OK, "Close")
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add((-1, 5))
        sizer.Add(icon, 0, wx.ALIGN_CENTER)
        sizer.Add((-1, 15))
        sizer.Add(info, 0, wx.EXPAND|wx.ALL)
        sizer.Add((-1, 10))
        sizer.Add(button, 0, wx.ALIGN_CENTER)
        self.SetSizer(sizer)

class AboutBoa(wx.Dialog):
    def __init__(self, parent):
        wx.Dialog.__init__(self, parent, -1, "About the Boa", size=(300, 300) )
        title = wx.StaticText(self, -1, "Elephant in a snake")
        pic = wx.StaticBitmap(self, -1, wx.Bitmap('boa0.jpg'), size=(100, 66) )
        caption = wx.StaticText(self, -1, "from The Little Prince\n"
            u"by Antoine de Saint Exupéry\n")
        button = wx.Button(self, wx.ID_OK, "Close")
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add((-1, 35))
        sizer.Add(title, 0, wx.ALIGN_CENTER)
        sizer.Add((-1, 25))
        sizer.Add(pic, 0, wx.ALIGN_CENTER)
        sizer.Add((-1, 25))
        sizer.Add(caption, 0, wx.ALIGN_CENTER)
        sizer.Add((-1, 15))
        sizer.Add(button, 0, wx.ALIGN_CENTER)
        self.SetSizer(sizer)

class PrintPage(wx.html.HtmlEasyPrinting):
    def __init__(self):
        wx.html.HtmlEasyPrinting.__init__(self)
    def Print(self, htmlfile):
        self.SetHeader(htmlfile)
        self.PrintFile(htmlfile)



