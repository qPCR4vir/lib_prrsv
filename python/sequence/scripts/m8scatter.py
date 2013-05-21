#!/usr/bin/env python


import math

import wx
import matplotlib
import numpy

import wxmpl

from sequence import blastNoSQL


class MainFrame(wx.Frame):
    """Main frame of the app
    """

    def __init__(self, *args, **kwargs):
        kwargs['style'] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwargs)

        self._filenameXAxis = None
        self._filenameYAxis = None
        self._xAxis = None
        self._yAxis = None
        
        self.__doLayout()

    def __createButton(self, label, handler):
        """Creates a button. Returns the button widget.
        
        Arguments:
        - `label`: Label for the button
        - `handler`: Event handler for the button
        """
        button = wx.Button(self, -1, label)
        self.Bind(wx.EVT_BUTTON, handler, button)
        return button

    def __createButtonAndTextCtrl(self, sizer, label, handler):
        """Creates a StaticText and a Button in the same row. Returns
        the text and button widgets.
        
        Arguments:
        - `sizer`: Sizer that will contain the horizontal sizer with
        both controls
        - `label`: Label for the button
        - `handler`: Event handler for the button
        """
        sizerContainer =  wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(sizerContainer, 0, wx.BOTTOM, 4)
        button = self.__createButton(label, handler)
        sizerContainer.Add(button, 0, wx.ALIGN_LEFT)
        text = wx.StaticText(self, -1)
        sizerContainer.Add(text, 1,
                           wx.EXPAND | wx.ALIGN_LEFT | wx.LEFT, 4)
        return text, button

    def __createMenuItem(self, parentMenu, label, handler=None):
        """Creates a menu item. Return the menu item widget. If label
        is 'SEPARATOR' then a separator is appended to `parentMenu`.
        
        Arguments:
        - `parentMenu`: Menu widget to attach to this menu item
        - `label`: Label for the menu item
        - `handler`: Event handler for the menu item
        """
        if label == 'SEPARATOR':
            parentMenu.AppendSeparator()
            return None
        
        menuItem = wx.MenuItem(parentMenu, -1, label, '',
                               wx.ITEM_NORMAL)
        parentMenu.AppendItem(menuItem)
        if handler is not None:
            self.Bind(wx.EVT_MENU, handler, menuItem)
        return menuItem

    def __doLayout(self):
        """Creates all the widgets.
        """
        menuBar = wx.MenuBar()
        fileMenu = wx.Menu()
        itemsInfo = [('Save as PNG...', self._onMenuSavePng),
                     ('Save as PS...', self._onMenuSavePs),
                     ('Print...', self._onMenuPrint),
                     ('SEPARATOR', None),
                     ('Quit', self._onMenuQuit),
                     ]
        for label, handler in itemsInfo:
            self.__createMenuItem(fileMenu, label, handler)
            
        menuBar.Append(fileMenu, 'File')
        self.SetMenuBar(menuBar)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        vSizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(vSizer, 1, wx.EXPAND, wx.ALL, 4)
        t, b = self.__createButtonAndTextCtrl(vSizer,
                                              'Select file (X-axis)',
                                              self._onSelectFileX)
        self._textFileX = t
        t, b = self.__createButtonAndTextCtrl(vSizer, 
                                              'Select file (Y-axis)',
                                              self._onSelectFileY)
        self._textFileY = t

        hSizer = wx.BoxSizer(wx.HORIZONTAL)
        vSizer.Add(hSizer, 0, wx.ALIGN_LEFT | wx.BOTTOM, 4)
        drawButton = self.__createButton('Draw', self._onDrawButton)
        hSizer.Add(drawButton, 0, wx.ALIGN_LEFT | wx.RIGHT, 4)
        swapButton = self.__createButton('Swap', self._onSwapButton)
        hSizer.Add(swapButton, 0, wx.ALIGN_LEFT)

        self.plotPanel = wxmpl.PlotPanel(self, -1, size=(6.0, 6.0),
                                         zoom=False, selection=False,
                                         location=False,
                                         crosshairs=False,
                                         cursor=False)
        vSizer.Add(self.plotPanel, 1, wx.EXPAND, wx.ALL, 4)
        
        self.SetSizer(sizer)
        sizer.Fit(self)

    def _onSelectFileX(self, event):
        """
        
        Arguments:
        - `event`:
        """
        self._filenameXAxis = self.openM8File()
        self._textFileX.SetLabel(self._filenameXAxis)
        self._xAxis = None
        
    def _onSelectFileY(self, event):
        """
        
        Arguments:
        - `event`:
        """
        self._filenameYAxis = self.openM8File()
        self._textFileY.SetLabel(self._filenameYAxis)
        self._yAxis = None

    def _onDrawButton(self, event):
        """
        
        Arguments:
        - `event`:
        """
        self.plotFigure()

    def _onSwapButton(self, event):
        """
        
        Arguments:
        - `event`:
        """
        if self._yAxis is not None and self._xAxis is not None:
            self._filenameXAxis, self._filenameYAxis = \
                                self._filenameYAxis, self._filenameXAxis 
            self._xAxis, self._yAxis = self._yAxis, self._xAxis
            self._textFileY.SetLabel(self._filenameYAxis)
            self._textFileX.SetLabel(self._filenameXAxis)
            self.plotFigure()
            
    def _onMenuPrint(self, event):
        """Handler of menu item 'Print'
        
        Arguments:
        - `event`:
        """
        printer = wxmpl.FigurePrinter(self.plotPanel)
        printer.printFigure(self.plotPanel.get_figure())

    def _onMenuQuit(self, event):
        """Handler of menu item 'Quit'
        
        Arguments:
        - `event`:
        """
        self.Close()

    def _onMenuSavePng(self, event):
        """Handler of menu item 'Save as PNG'.
        
        Arguments:
        - `event`:
        """
        self.saveFigure('png')

    def _onMenuSavePs(self, event):
        """Handler of menu item 'Save as PostScript'.
        
        Arguments:
        - `event`:
        """
        self.saveFigure('ps')

    def saveFigure(self, fmt='png'):
        """Shows a file selection dialog and saves the current figure
        in the selected file.
        
        Arguments:
        - `fmt`: Format to use. Any format supported by matplotlib
        (png, pdf, ps, eps and svg)
        """
        dlg = wx.FileDialog(self, style=wx.SAVE | wx.OVERWRITE_PROMPT)
        if dlg.ShowModal() == wx.ID_OK:
            figure = self.plotPanel.get_figure()
            filename = dlg.GetPath()
            if not filename.endswith(fmt):
                filename = '%s.%s' % (filename, fmt)
            figure.savefig(filename, format=fmt)
        dlg.Destroy()

    def openM8File(self):
        """Shows a file selection dialog and returns a generator of
        tuples with the content of the selected file.
        """
        dlg = wx.FileDialog(self, style=wx.OPEN)
        filename = None
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
        dlg.Destroy()
        return filename

    def plotFigure(self):
        """Plots the figure using the m8generator
        """
        if self._filenameXAxis is None or \
               self._filenameYAxis is None:
            return
        figure = self.plotPanel.get_figure()
        figure.clear()
        axis = figure.add_subplot(111)

        if self._xAxis is None or self._yAxis is None:
            querys = {}
            lastIndex = -1

            generator = blastNoSQL.m8tupleGenerator(self._filenameXAxis)
            eValues = []
            for row in generator:
                if not querys.has_key(row[0]):
                    lastIndex += 1
                    querys[row[0]] = lastIndex
                eValues.append(row[-2])
            maxX = max(eValues)
            minX = min(eValues)

            generator = blastNoSQL.m8tupleGenerator(self._filenameYAxis)
            eValues = []
            for row in generator:
                if not querys.has_key(row[0]):
                    lastIndex += 1
                    querys[row[0]] = lastIndex
                eValues.append(row[-2])
            maxY = max(eValues)
            minY = min(eValues)

            xs = numpy.empty((1,))
            xs[0] = maxX
            xs = xs.repeat(len(querys))
            ys = numpy.empty((1,))
            ys[0] = maxY
            ys = ys.repeat(len(querys))

            generator = blastNoSQL.m8tupleGenerator(self._filenameXAxis)
            for row in generator:
                value = row[-2]
                index = querys[row[0]]
                xs[index] = min(value, xs[index])

            generator = blastNoSQL.m8tupleGenerator(self._filenameYAxis)
            for row in generator:
                value = row[-2]
                index = querys[row[0]]
                ys[index] = min(value, ys[index])

            xs = -1 * numpy.log10(xs)
            ys = -1 * numpy.log10(ys)
            self._xAxis = xs
            self._yAxis = ys
            
        collection = axis.hexbin(self._xAxis, self._yAxis, bins='log',
                                 mincnt=1)
        xMin, xMax = axis.get_xlim()
        yMin, yMax = axis.get_ylim()
        globalMin = min(xMin, yMin)
        globalMax = max(xMax, yMax)
        axis.set_xlim(globalMin, globalMax)
        axis.set_ylim(globalMin, globalMax)
        figure.colorbar(collection)
        self.plotPanel.draw()


if __name__ == '__main__':
    app = wx.App(False)
    mainFrame = MainFrame(None)
    mainFrame.Show(True)
    app.MainLoop()
