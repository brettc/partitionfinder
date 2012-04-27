import wx
from newx import cmd, toolbar, menu
from newx import autobind, layout, split, persist, events
from newx.log import *
from newx.idfree import *

from wx.lib.evtmgr import eventManager as evt

cmd.add(FILE_EXIT=cmd.Command('E&xit', help='Exit the application'))

# our standard mainframe
# TODO -- get the size from persist
class MainFrame(wx.Frame):

    def __init__(self, title="Please Give a Title"):

        # size = persist.make('mainframe_size', (500, 500), self, wx.EVT_MOVE, self.GetSizeTuple)
        # pos = persist.make('mainframe_pos', (0, 0), self, wx.EVT_SIZE, self.GetPositionTuple)

        wx.Frame.__init__(self, None, -1, title, size=(500, 500),
                        style=wx.DEFAULT_FRAME_STYLE|wx.NO_FULL_REPAINT_ON_RESIZE)

        # need to first register commands
        # self.Setup()

        if hasattr(self, '__wxtoolbar__'):
            self.toolbar = toolbar.make_toolbar(self, self.__wxtoolbar__)
            self.SetToolBar(self.toolbar)
            self.toolbar.Realize()

        if hasattr(self, '__wxmenu__'):
            self.menubar = menu.make_bar(self.__wxmenu__)
            self.SetMenuBar(self.menubar)

        # statusbar
        self.CreateStatusBar()

        # finished, now layout
        # if hasattr(self, '__wxlayout__'):
            # fmt = self.__wxlayout__
            # grid = layout.parse_format(fmt)
            # self.sizer = grid.make_sizer(self) 
            # self.SetSizerAndFit(self.sizer)

        # lastly bind all commands
        cmd.bind(self)
        # and autobind win funcs and children
        # autobind.bind(self)
        # self.Bind(wx.EVT_CLOSE, self.OnClose)

        # server.subscribe(topic='mainframe', listener=self.notify)
        # evt.Register(self.OnSize, wx.EVT_SIZE, self)
        # evt.Register(self.OnMove, wx.EVT_MOVE, self)

        autobind.bind(self)
        autobind.bindcontrols(self)

        # if hasattr(self, '__wxlayout__'):
            # gbs = self.__wxlayout__.make_sizer(self)
            # self.SetAutoLayout(1)
            # self.SetSizerAndFit(gbs)
            # self.Layout()

    # def Layout(self):
        # g = layout.parse_format(self.__wxlayout__)
        # gbs = g.make_sizer(self)
        # self.SetSizerAndFit(gbs)
        # autobind.bindcontrols(self)

    def OnSize(self, event):
        event.Skip()
        # persist.set('mainframe_size', self.GetSizeTuple())

    def OnMove(self, event):
        event.Skip()
        # persist.set('mainframe_pos', self.GetPositionTuple())

    # def notify(self, message):
        # # print 'got', message
        # if message.topic[1] == 'size':
            # self.SetSize(message.data)
        # elif message.topic[1] == 'pos':
            # self.SetPosition(message.data)


    def OnFileExit(self, event):
        self.Close()

    def OnClose(self, event):
        # persist.update()
        # persist.save('xxxx')
        self.Destroy()

class Notebook(wx.Notebook):

    def __init__(self, parent):
        wx.Notebook.__init__(self, parent, wx.NewId(), style=wx.NB_BOTTOM)
        self.AddPages()

        # Tell the pages they have been changed -- this is much more useful
        # than having it in the Notebook. IMHO.
        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self._PageChanged)

    def AddPages(self):
        raise NotImplementedError
        
    def _PageChanged(self, event):
        frompage = event.GetOldSelection()
        topage = event.GetSelection()
        fromname, toname = '<nopage>', '<nopage>'

        if frompage >= 0:
            pg = self.GetPage(frompage)
            fromname = self.GetPageText(frompage)
            if hasattr(pg, 'OnPageClosed'):
                pg.OnPageClosed()
        if topage >= 0:
            pg = self.GetPage(topage)
            toname = self.GetPageText(topage)
            if hasattr(pg, 'OnPageOpened'):
                pg.OnPageOpened()

        debug('changing pages in notebook from %s to %s', fromname, toname)
        event.Skip()
        events.publish("newx.page")


class Panel(wx.Panel):

    def __init__(self, parent):
        wx.Panel.__init__(self, parent, wx.NewId())

    def Layout(self):
        if hasattr(self, '__wxlayout__'):
            gbs = self.__wxlayout__.make_sizer(self)
            self.SetAutoLayout(1)
            self.SetSizerAndFit(gbs)
            # self.Layout()
        autobind.bindcontrols(self)
    
class Dialog(wx.Dialog):

    def __init__(self, parent, *args, **kwargs):
        if not kwargs.has_key('style'):
            kwargs['style'] = wx.DEFAULT_DIALOG_STYLE
        wx.Dialog.__init__(self, parent, wx.NewId(), *args, **kwargs)
        self.Layout()

    def Layout(self):
        if hasattr(self, '__wxlayout__'):
            gbs = self.__wxlayout__.make_sizer(self)
            self.SetAutoLayout(1)
            self.SetSizerAndFit(gbs)
            # self.Layout()
        # g = layout.parse_format(self.__wxlayout__)
        # gbs = g.make_sizer(self)
        # self.SetSizerAndFit(gbs)
        # magically bind controls
        autobind.bindcontrols(self)

class ImageDisplayWindow(wx.ScrolledWindow):

    def __init__(self, parent, **kw):
        more_style = 0
        if kw.has_key('style'):
            more_style = kw['style']

        wx.ScrolledWindow.__init__(self, 
            parent, style=wx.VSCROLL|wx.HSCROLL|more_style)
        self.SetScrollbars(1, 1, 50, 50)
        self.bmp = None
        self.SetBackgroundColour(wx.WHITE)
        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnEraseBackground)

    def Update(self, msg):
        # msg.data should be a filename
        #
        try:
            self.bmp = wx.Image(msg.data).ConvertToBitmap()
        except:
            error('could not load %s', msg.data)
            self.bmp = None

        if self.bmp:
            self.SetVirtualSize((self.bmp.GetWidth(), self.bmp.GetHeight()))
        self.Refresh()

    def OnPaint(self, pdc):
        pdc = wx.PaintDC(self)
        self.PrepareDC(pdc)
        if self.bmp:
            pdc.DrawBitmap(self.bmp, 0, 0)

    def OnEraseBackground(self, event):
        dc = event.GetDC()
        dc.Clear()
  

class App(wx.App):
    def __init__(self, create_frame_fun, *args, **kwargs):
        # persist.load('xxxx')
        self.create_frame_fun = create_frame_fun
        self.winargs = args
        self.winkwargs = kwargs
        wx.App.__init__(self, 0)

    def OnInit(self):
        # pbly NOT needed, though possible
        # cmd.bind(self)
        # print self.create_frame_fun
        # print '-----'
        self.frame = self.create_frame_fun(*self.winargs, **self.winkwargs)
        self.SetTopWindow(self.frame)
        # persist.broadcast()
        self.frame.Show(True)
        return True

    def OnFileExit(self, event):
        self.frame.Close()

def run(mainwin_class, *args, **kwargs):
    app = App(mainwin_class, *args, **kwargs)
    app.MainLoop()


# smarter about what they are doing
# maybe something like this
