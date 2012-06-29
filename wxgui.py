#
# Many ideas are taken from here:
# http://wiki.wxpython.org/Optimizing%20for%20Mac%20OS%20X
#
#
# TODO
# * Set a flag for breaking out of the analysis cleanly (make it part of the
# progress monitor?)
# * Sort out the logging later
# * Editing, on another tab
# * Multiple windows, or not... 
# * Use the dv.PyDataViewIndexListModel Example for the logger
# * How to resize vertically!

import logging
import thread
log = logging.getLogger("gui")
import wx
import wx.stc as stc

from partfinder import config, analysis_method, reporter

myformatter = logging.Formatter(
    """%(levelname)s: %(name)s, %(asctime)s,%(filename)s, %(lineno)s: %(message)s\n"""
)

# Make Publishing easier
from wx.lib.pubsub import Publisher
publisher = Publisher()

def publish(path, data=None):
    assert isinstance(path, str)
    topic = tuple(path.split('.'))
    publisher.sendMessage(topic, data)

def subscribe(path, listener):
    assert isinstance(path, str)
    topic = tuple(path.split('.'))
    publisher.subscribe(listener, topic)

class Handler(logging.Handler):
    def emit(self, record):
        self.output.AddText(myformatter.format(record))

class BasePanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        self.sizer = self.make_contents(parent)
        self.SetSizer(self.sizer)
        self.sizer.Fit(self)
        self.Fit()

class DetailsPane(wx.Panel):
    def __init__(self, parent):
        pass

class MainPanel(BasePanel):
    def make_contents(self, parent):
        box = wx.BoxSizer()

        self.path = wx.StaticText(parent, -1, "?")
        box.Add(self.path, wx.EXPAND|wx.ALL, 5)

        self.button = wx.Button(parent, -1, "...")
        box.Add(self.button, wx.EXPAND|wx.ALL, 5)
        parent.Bind(wx.EVT_BUTTON, lambda x: publish('button.press', x), self.button)

        subscribe('data.change', self.DataUpdate)
        self.SetMinSize((400, 400))
        return box

    def DataUpdate(self, msg):
        self.path.SetLabel(msg.data.base_path)


class MainFrame(wx.Frame):
    def __init__(self, title = "Partition Finder"):
        wx.Frame.__init__(self, None , -1, title)

        self.CreateMenu()
        # self.CreateContent()
        self.CreateStatusBar()

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.panel = MainPanel(self)
        # self.Bind(wx.EVT_BUTTON, self.OnButton, self.panel.button)
        self.sizer.Add(self.panel, 1, wx.EXPAND)
        self.SetSizer(self.sizer)

        self.Fit()
        self.SetMinSize(self.GetSize())

        self.cfg = config.Configuration()
        subscribe('button.press', self.GetPath)

        self.running = None

    def CreateMenu(self):

        MenuBar = wx.MenuBar()
        FileMenu = wx.Menu()
        item = FileMenu.Append(wx.ID_EXIT, text = "&Exit")
        self.Bind(wx.EVT_MENU, self.OnQuit, item)
        item = FileMenu.Append(wx.ID_ANY, text = "&Open")
        self.Bind(wx.EVT_MENU, self.OnOpen, item)
        item = FileMenu.Append(wx.ID_PREFERENCES, text = "&Preferences")
        self.Bind(wx.EVT_MENU, self.OnPrefs, item)
        MenuBar.Append(FileMenu, "&File")
        
        HelpMenu = wx.Menu()
        item = HelpMenu.Append(wx.ID_HELP, "Test &Help",
                                "Help for this simple test")
        self.Bind(wx.EVT_MENU, self.OnHelp, item)
        ## this gets put in the App menu on OS-X
        item = HelpMenu.Append(wx.ID_ABOUT, "&About",
                                "More information About this program")
        self.Bind(wx.EVT_MENU, self.OnAbout, item)
        MenuBar.Append(HelpMenu, "&Help")

        self.SetMenuBar(MenuBar)

    # def OnPressed(self, evt):
        # if self.running is None:
            # self.running = PFThread()
            # self.running.Start()

    def GetPath(self, evt):
        # In this case we include a "New directory" button. 
        dlg = wx.DirDialog(self, "Choose a directory:",
                        style=wx.DD_DEFAULT_STYLE
                        #| wx.DD_DIR_MUST_EXIST
                        #| wx.DD_CHANGE_DIR
                        )

        # If the user selects OK, then we process the dialog's data.
        # This is done by getting the path data from the dialog - BEFORE
        # we destroy it. 
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            try:
                self.cfg.load_base_path(path)
            except:
                pass
            publish('data.change', self.cfg)

        # Only destroy a dialog after you're done with it.
        dlg.Destroy()

        
    def OnQuit(self,Event):
        self.Destroy()
        
    def OnAbout(self, event):
        dlg = wx.MessageDialog(self, "This is a small program to test\n"
                                     "the use of menus on Mac, etc.\n",
                                "About Me", wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

    def OnHelp(self, event):
        dlg = wx.MessageDialog(self, "This would be help\n"
                                     "If there was any\n",
                                "Test Help", wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

    def OnOpen(self, event):
        dlg = wx.MessageDialog(self, "This would be an open Dialog\n"
                                     "If there was anything to open\n",
                                "Open File", wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

    def OnPrefs(self, event):
        dlg = wx.MessageDialog(self, "This would be an preferences Dialog\n"
                                     "If there were any preferences to set.\n",
                                "Preferences", wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()
        
class App(wx.App):
    def __init__(self, *args, **kwargs):
        wx.App.__init__(self, *args, **kwargs)
        # This catches events when the app is asked to activate by some other
        # process
        self.Bind(wx.EVT_ACTIVATE_APP, self.OnActivate)

    def OnInit(self):
        frame = MainFrame()
        frame.Show()

        import sys
        for f in  sys.argv[1:]:
            self.OpenFileMessage(f)

        return True

    def BringWindowToFront(self):
        try: # it's possible for this event to come when the frame is closed
            self.GetTopWindow().Raise()
        except:
            pass
        
    def OnActivate(self, event):
        # if this is an activate event, rather than something else, like iconize.
        if event.GetActive():
            self.BringWindowToFront()
        event.Skip()
    
    def OpenFileMessage(self, filename):
        dlg = wx.MessageDialog(None,
                               "This app was just asked to open:\n%s\n"%filename,
                               "File Dropped",
                               wx.OK|wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()

    # def MacOpenFile(self, filename):
        # """Called for files droped on dock icon, or opened via finders context menu"""
        # print filename
        # print "%s dropped on app"%(filename) #code to load filename goes here.
        # self.OpenFileMessage(filename)
        
    def MacReopenApp(self):
        """Called when the doc icon is clicked, and ???"""
        self.BringWindowToFront()

    def MacNewFile(self):
        pass
    
    def MacPrintFile(self, file_path):
        pass


class PFThread(object):
    def __init__(self):
        pass

    def Start(self):
        self.keepGoing = self.running = True
        thread.start_new_thread(self.Run, ())

    def Stop(self):
        self.keepGoing = False

    def IsRunning(self):
        return self.running

    def Run(self):
        # logging.getLogger("").addHandler(gui_logging.gui_handler) # add at base level
        # logging.getLogger("").setLevel(logging.INFO)
        cfg = config.Configuration()
        cfg.load_base_path('example')
        method = analysis_method.choose_method(cfg.search)
        rpt = reporter.TextReporter(cfg)
        anal = method(cfg, rpt, False, False)
        anal.analyse()
 
app = App(False)
app.MainLoop()
