#
# Many ideas are taken from here:
# http://wiki.wxpython.org/Optimizing%20for%20Mac%20OS%20X
#

import logging
import thread
log = logging.getLogger("gui")
import wx
import wx.lib.sized_controls as sc

import gui_logging 
from gui_logging import LogHtmlListBox, LogListCtrl

class MainFrame(sc.SizedFrame):
    """ This window displays a button """
    def __init__(self, title = "Partition Finder"):
        sc.SizedFrame.__init__(self, None , -1, title)

        self.CreateMenu()
        self.CreateContent()
        # self.CreateStatusBar()

        self.Fit()
        self.SetMinSize(self.GetSize())
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


    def CreateContent(self):

        pane = self.GetContentsPane()
        # pane.SetSizerType("grid", {"cols":3}) # 3-column grid layout
        # 
        # # row 1
        # wx.TextCtrl(pane, -1).SetSizerProps(halign="left")
        # wx.TextCtrl(pane, -1).SetSizerProps(halign="center")
        # wx.TextCtrl(pane, -1).SetSizerProps(halign="right")
        # 
        # # row 2
        # wx.TextCtrl(pane, -1).SetSizerProps(valign="center")
        # wx.TextCtrl(pane, -1).SetSizerProps(expand=True, proportion=1)
        # wx.TextCtrl(pane, -1).SetSizerProps(valign="center")
        # 
        # # row 3
        # wx.TextCtrl(pane, -1).SetSizerProps(halign="left")
        # wx.TextCtrl(pane, -1).SetSizerProps(halign="center")
        # wx.TextCtrl(pane, -1).SetSizerProps(halign="right")

        pane.SetSizerType("grid", options={'cols':2})
        
        # row 1
        wx.StaticText(pane, -1, "Dir")
        self.button = wx.Button(pane, -1, "...")
        self.Bind(wx.EVT_BUTTON, self.OnButton, self.button)
        # textCtrl = wx.TextCtrl(pane, -1, "Your name here")
        self.button.SetSizerProps(expand=True)
        
        # row 2
        wx.StaticText(pane, -1, "Email")
        self.pressme = wx.Button(pane, -1, "push me")
        self.Bind(wx.EVT_BUTTON, self.OnPressed, self.pressme)
        
        # row 3
        wx.StaticText(pane, -1, "Gender")
        wx.Choice(pane, -1, choices=["male", "female"])
        
        # row 4
        wx.StaticText(pane, -1, "State")
        wx.TextCtrl(pane, -1, size=(60, -1)) # two chars for state
        
        # row 5
        # radioPane = sc.SizedPanel(pane, -1)
        wx.StaticText(pane, -1, "Log")
        html = LogHtmlListBox(pane, size=(60, 100))
        html.SetSizerProps(expand=True)
        # loglist = LogListCtrl(pane)
        # loglist.SetSizerProps(expand=True)
        
        # here's how to add a 'nested sizer' using sized_controls
        # radioPane = sc.SizedPanel(pane, -1)
        # radioPane.SetSizerType("horizontal")
        # radioPane.SetSizerProps(expand=True)
        
        # make these children of the radioPane to have them use
        # the horizontal layout
        # wx.RadioButton(radioPane, -1, "Mr.")
        # wx.RadioButton(radioPane, -1, "Mrs.")
        # wx.RadioButton(radioPane, -1, "Dr.")
        # end row 5

    def OnPressed(self, evt):
        if self.running is None:
            self.running = PFThread()
            self.running.Start()

    def OnButton(self, evt):
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
            self.button.SetLabel(dlg.GetPath())

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

        iconFile = "Resources/pf512.ico"
        icon = wx.Icon(iconFile, wx.BITMAP_TYPE_ICO)
        tbicon = wx.TaskBarIcon()
        tbicon.SetIcon(icon, "I am an Icon")
        
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
        log.warning("FUCK!")
        self.keepGoing = self.running = True
        thread.start_new_thread(self.Run, ())

    def Stop(self):
        self.keepGoing = False

    def IsRunning(self):
        return self.running

    def Run(self):
        from partfinder import config, analysis_method, reporter
        logging.getLogger("").addHandler(gui_logging.gui_handler) # add at base level
        logging.getLogger("").setLevel(logging.INFO)
        cfg = config.Configuration()
        cfg.load_base_path('example')
        method = analysis_method.choose_method(cfg.search)
        rpt = reporter.TextReporter(cfg)
        anal = method(cfg, rpt, False, False)
        anal.analyse()
 
app = App(True)
app.MainLoop()
