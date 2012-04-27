import logging
import wx
import log, cmd, easy, persist, events, autobind

from log import gui_handler

# FIXME
# would be good for exceptions etc
# Need to make an HTML formatter then...
# or use the VListBox from the same Demo

htmlformatter = logging.Formatter(
    '''
    <font size=2>
        <b><font color='blue'>%(levelname)s: %(name)s</font></b>
        <font color='red'>%(asctime)s</font>
        <br>
        %(filename)s, %(lineno)s: %(message)s
    </font>
    ''')
    # <pre>
    # blarg is 
    # oeueu
    # </pre>

class LogHtmlListBox(wx.HtmlListBox):
    def __init__(self, parent, **kwargs):
        wx.HtmlListBox.__init__(self, parent, **kwargs)
        # gui_handler.add_watch(self)
        autobind.bind(self)
        self.Update()
        events.subscribe('newx.log.message', self.Update)
    
    def Update(self, evt=None):
        global gui_handler
        self.SetItemCount(len(gui_handler))
        self.ScrollToLine(self.GetLineCount())
        self.Refresh()

    def OnGetItem(self, n):
        global gui_handler
        rec = gui_handler[n]
        return htmlformatter.format(rec)

    def OnDrawSeparator(self, dc, rect, n):
        # Place some simple separation between each of the items.

        dc.SetBrush(wx.NullBrush)
        dc.SetPen(wx.GREY_PEN)
        dc.DrawLine(rect.left, rect.bottom, rect.right, rect.bottom)

    def OnDblClicked(self, evt):
        # We are going to extract the file name and line number and send them
        # out over the events system. The idea is that you can pick this up, 
        # and then send the user to the offending line in an editor.

        n = evt.GetInt()
        rec = gui_handler[n]

        # This has been modified from the traceback module in the standard 
        # library.
        if rec.exc_info:
            etype, value, tb = rec.exc_info
            print dir(value), value
            # If it is a Syntax Error get it from there, instead of the error 
            # that the line occured on. This will be when someone is using eval or 
            # execfile, or some late import.
            if etype is SyntaxError:
                try:
                    msg, (filename, lineno, offset, line) = value
                except:
                    pass
            else:
                # Get the path and line for the exception. We trace down the 
                # stack to the offending line
                while tb is not None:
                    f = tb.tb_frame
                    lineno = tb.tb_lineno
                    co = f.f_code
                    filename = co.co_filename
                    tb = tb.tb_next
        else:
            filename = rec.pathname
            lineno = rec.lineno

        # Tell anyone who is interested...
        events.publish('newx.log.source', (filename, lineno))

class LogPanel(easy.Panel):

    __wxlayout__ = '''
    P : 1
    1 : loglist
    '''

    def __init__(self, parent):
        easy.Panel.__init__(self, parent)
        self.loglist = LogHtmlListBox(self, style=wx.SUNKEN_BORDER)
        self.Layout()

    # def OnDrawSeparator(dc, rect, n):
        # print 'sep'
        # pass


class LogListCtrl(wx.ListCtrl):
    '''For showing incoming log messages'''

    COLUMNS = 'name levelname message'.split()

    def __init__(self, parent):
        global gui_handler
        wx.ListCtrl.__init__(self, parent, wx.NewId(), size=(300,200),
                style = wx.NO_BORDER | wx.LC_REPORT | wx.LC_VIRTUAL)

        self.SetItemCount(len(gui_handler))

        self.lookup = {}
        for c, nm in enumerate(self.COLUMNS):
            self.InsertColumn(c, nm, wx.LIST_FORMAT_LEFT)
            self.lookup[c] = nm

        # gui_handler.add_watch(self)
        events.subscribe('newx.log.message', self.Update)

    def Update(self, evt=None):
        global gui_handler
        items = len(gui_handler)
        self.SetItemCount(items)
            # setattr(self, nm, c)
        self.Focus(items-1)

    #---------------------------------------------------
    # These methods are callbacks for implementing the
    # "virtualness" of the list...
    def OnGetItemText(self, item, col):
        global gui_handler
        return getattr(gui_handler[item], self.lookup[col])

    def OnGetItemImage(self, item):
        return 0

    def OnGetItemAttr(self, item):
        return None


class LogFrame(wx.MiniFrame):
    '''Holds a LogListCtrl'''

    def __init__(self, parent):
        # size = persist.make('logframe_size', (500, 500), self.GetSizeTuple)
        # pos = persist.make('logframe_pos', (0, 0), self.GetPositionTuple)
        wx.MiniFrame.__init__(self, parent, wx.NewId(), "Log", 
            size=(200, 200))
        self.loglist = LogHtmlListBox(self)
        sizer = wx.BoxSizer()
        sizer.Add(self.loglist, 1, wx.GROW |wx.ALL, 0)
        self.SetAutoLayout(1)
        self.SetSizer(sizer)
        sizer.SetSizeHints(self)
        self.Layout()
        # sizer.SetSizeHints(self)
        #
    #FIXME Shouldnt' really close from button
    # EVT_CLOSE(func) 
    # def OnClose(self, event):
        # print event


class ShowLogCmd(cmd.ToggleCmd):
    '''Allows toggling of the Frame Gui'''

    def __init__(self, text):
        cmd.ToggleCmd.__init__(self, text)
        self.logwin = None
    
    def get_on(self):
        # if not visible assume that not split
        return self.logwin and self.logwin.IsShown()

    def set_on(self, v):
        self.show_log(v)

    on = property(get_on, set_on)

    def show_log(self, v):
        if v:
            if not self.logwin:
                tw = wx.GetApp().GetTopWindow()
                self.logwin = LogFrame(tw)
            self.logwin.Show(True)
        else:
            if self.logwin:
                self.logwin.Show(False)

cmd.add(SHOW_LOG=ShowLogCmd('Show &Log|Show the Log Window'))
