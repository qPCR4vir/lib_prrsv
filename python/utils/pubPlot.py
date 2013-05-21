#!/usr/local/bin/python

import wx
import wxmpl
import matplotlib

class TopPanel(wxmpl.PlotPanel):
    pass
    #def __init__(self, parent, id):
        #wx.Panel.__init__(self, parent, id, style=wx.BORDER_SUNKEN)
    #    pass

class BottomPanel(wx.Panel):

    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id)
        self.prompt = wx.StaticText(self, -1, 'Search Term:', (40, 60))
        #self.term = wx.TextCtrl(self,-1,"")
        self.goBtn = wx.Button(self,-1,"Update Results")


        vbox = wx.BoxSizer(wx.VERTICAL)
        hbox = wx.BoxSizer(wx.HORIZONTAL)

        hbox.Add(self.prompt,0,wx.ALIGN_RIGHT)
        #hbox.Add(self.term,0,wx.EXPAND)
        
        vbox.Add(self.goBtn)
        vbox.Add(hbox)


class SearchHistory(wx.Frame):
    def __init__(self, parent, id, title):
        wx.Frame.__init__(self, parent, id, title, size=(600, 600))

        self.topPanel = TopPanel(self,-1)
        self.bottomPanel = BottomPanel(self, -1)


        vbox = wx.BoxSizer(wx.VERTICAL)
        vbox.Add(self.topPanel, 1, wx.EXPAND | wx.ALL, 5)
        vbox.Add(self.bottomPanel, 1, wx.EXPAND | wx.ALL, 5)

        self.menu = MenuBar()
        self.SetMenuBar(self.menu)
        self.Bind(wx.EVT_MENU,self.OnQuit,self.menu.quit)
         
        self.SetSizer(vbox) 
        self.Show(True)

    def OnQuit(self,event):
        """
        
        Arguments:
        - `self`:
        """
        self.Close()
 

class MenuBar(wx.MenuBar):
    """
    """

    def __init__ (self):
        wx.MenuBar.__init__(self)
        self.file=wx.Menu()
        self.quit = self.file.Append(-1,'Quit','1 guess!')
        self.Append(self.file,'&File')



def main():
    global app,frame
    app = wx.App()
    frame=SearchHistory(None, -1, 'PubMed History')
    app.MainLoop()    

if __name__ == '__main__':
    main()
