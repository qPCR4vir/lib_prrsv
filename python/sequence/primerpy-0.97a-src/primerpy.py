#!/usr/bin/env python
"""primerpy.py, main GUI frame of primerpy-tortilla

PrimerPy version 0.97a (tortilla)
Copyright 2007, 2008, Shuzhao Li
http://www.bioinformatics.org/primerpy

PrimerPy is released under GNU/GPL v2.0 or newer;
provided as is, with no guanrantee of any form implied.

The current release uses primer3 as a backend engine.
Primer3 is copyrighted by
Whitehead Institute for Biomedical Research, Steve Rozen
(http://jura.wi.mit.edu/rozen), and Helen Skaletsky.
Primer3 is under GNU/GPL and the new BSD license.

"""

from guiclasses import *

class MyMainFrame(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(600, 450))
        self.CreateMenu()
        # PrimerNoteBook, main organizer
        self.nb = PrimerNoteBook(self)
        self.StatusBar()
        self.Centre()
        self.Show()
        self.SanityCheck()

    def CreateMenu(self):
        menu0 = wx.Menu()
        item = menu0.Append(-1, "Design")
        self.Bind(wx.EVT_MENU, self.OnDesign, item)
        item = menu0.Append(-1, "Check dimers")
        self.Bind(wx.EVT_MENU, self.OnDimers, item)
        menu0.AppendSeparator()
        #
        item = menu0.Append(-1, "Print")
        #self.Bind(self.nb.page_changed, ?)
        self.Bind(wx.EVT_MENU, self.Print, item)
        #
        menu0.AppendSeparator()
        item = menu0.Append(-1, "Configure")
        self.Bind(wx.EVT_MENU, self.OnConfigure, item)
        item = menu0.Append(-1, "Exit")
        self.Bind(wx.EVT_MENU, self.OnExit, item)
        # display
        menu1 = wx.Menu()
        item = menu1.Append(-1, "Last Run")
        self.Bind(wx.EVT_MENU, self.OnResult, item)
        # original primer3 output?
        item = menu1.Append(-1, "Log")
        self.Bind(wx.EVT_MENU, self.OnLog, item)
        menu1.AppendSeparator()
        Manual = menu1.Append(-1, "Manual")
        self.Bind(wx.EVT_MENU, self.OnManual, Manual)
        # about
        menu2 = wx.Menu()
        item = menu2.Append(-1, "About PrimerPy")
        self.Bind(wx.EVT_MENU, self.OnAbout, item)
        item = menu2.Append(-1, "About the boa")
        self.Bind(wx.EVT_MENU, self.OnBoa, item)

        menubar = wx.MenuBar()
        menubar.Append(menu0, "&PrimerPy")
        menubar.Append(menu1, "&Show")
        menubar.Append(menu2, "&About")
        self.SetMenuBar(menubar)

    def OnConfigure(self, *event):
        abox = configure_dialog(self)
        abox.ShowModal()
        abox.Destroy()
    def OnExit(self, event):
        self.Close(True)

    def OnDesign(self, event):
        self.nb.show_design()
    def OnDimers(self, event):
        self.nb.show_dimercheck()
    def Print(self, event):
        "print the page in selection if html"
        page = self.nb.GetSelection()
        try:
            htmlfile = self.nb.GetPage(page).OpenedPage
            p = PrintPage()
            p.Print(htmlfile)
        except AttributeError:
            self.statusbar.SetStatusText("Can only print text results.")

    def OnResult(self, event):
        self.nb.show_result()
    def OnLog(self, event):
        self.nb.show_log()
    def OnManual(self, event):
        self.nb.show_manual()

    def OnAbout(self, event):
        abox = AboutPrimerPy(self)
        abox.ShowModal()
        abox.Destroy()
    def OnBoa(self, event):
        abox = AboutBoa(self)
        abox.ShowModal()
        abox.Destroy()

    def StatusBar(self):
        self.statusbar = self.CreateStatusBar()
    def SanityCheck(self):
        global UserDataDir, Primer3_PATH
        p = Primer3_PATH.split(' ')[0]
        if 'win32' in sys.platform:
            p += '.exe'
        if not os.path.isfile(p) or not os.path.isdir(UserDataDir):
            wx.FutureCall(5000, self.OnConfigure())


#----------------------------------------------------------------------

if __name__ == '__main__':

    class MyApp(wx.App):
        def OnInit(self):
            frame = MyMainFrame(None, -1, 'PrimerPy, QPCR primer design')
            frame.Show(True)
            self.SetTopWindow(frame)
            return True

    app = MyApp(False)
    app.MainLoop()

