import wx
import cmd

# TODO SetMinimumPaneSize based on the size of the windows that are being
# inserted?

# Window styles
# wx.SP_3D  Draws a 3D effect border and sash.  
# wx.SP_3DSASH  Draws a 3D effect sash.  
# wx.SP_3DBORDER  Synonym for wx.SP_BORDER.  
# wx.SP_BORDER  Draws a standard border.  
# wx.SP_NOBORDER  No border (default).  
# wx.SP_NO_XP_THEME  Under Windows XP, switches off the attempt to draw the
# splitter using Windows XP theming, so the borders and sash will take on the
# pre-XP look.  
# wx.SP_PERMIT_UNSPLIT  Always allow to unsplit, even with the minimum pane
# size other than zero.  
# wx.SP_LIVE_UPDATE  Don't draw XOR line but resize the child windows
# immediately.  

# souped up version
class Splitter(wx.SplitterWindow):
    def __init__(self, parent, style):
        wx.SplitterWindow.__init__(self, parent, wx.NewId(), style=style)

        self.Bind(wx.EVT_SPLITTER_SASH_POS_CHANGED, self.OnSashPosChanged) 
        self.Bind(wx.EVT_SPLITTER_DCLICK, self.OnDoubleClickSash) 
        self.sashpos = 100

    def OnDoubleClickSash(self, event):
        # DAMN How to stop this fucking thing!
        # event.Skip()
        pass

    def OnSashPosChanged(self, event):
        self.sashpos = event.GetSashPosition()
        # print 'sash now at ', event.GetSashPosition()
        #
        #
    def split(self):
        if self.splitmode == wx.SPLIT_VERTICAL:
            splitfun = self.SplitVertically
        else:
            splitfun = self.SplitHorizontally
        # don't need self as it is already bound
        splitfun(*self.wins)

    # eventManager.Register(printEvent,        wx.EVT_TOGGLEBUTTON, button)
    # eventManager.Register(enableFrameEvents, wx.EVT_TOGGLEBUTTON, button)


def make_split(parent, create1, create2, splitmode=wx.SPLIT_VERTICAL,
        style=wx.SP_LIVE_UPDATE|wx.SP_3D, splitcls=Splitter):
    '''Make a split window in the parent

    Creates in the parent, using the result of calling the 2 create functions
    with parent as the first param ie. win1 = create1(parent)

    name gives us a unique identifier for the general purposes
    '''
    splitter = splitcls(parent, style)
    splitter.wins = create1(splitter), create2(splitter)
    splitter.splitmode = splitmode
    splitter.split()

    return splitter

def make_split_cmd(splitter, win, **kw):
    for cmdname, txt in kw.items():
        scmd = SplitCmd(splitter, win, txt)
        cmd.add(**{cmdname : scmd})
        break

class SplitCmd(cmd.ToggleCmd):
    '''A command class for handling split and unsplit events automatically'''

    def __init__(self, splitter, num, text):
        assert num == 0 or num == 1

        cmd.ToggleCmd.__init__(self, text)
        self.splitter = splitter
        self.num = num
        self.win = splitter.wins[num]

    def get_on(self):
        # if not visible assume that not split
        return self.win.IsShown()

    def set_on(self, v):
        self.do_split(v)
    on = property(get_on, set_on)

    def update(self, event):
        # disable this if we are the only one visible
        on = self.on
        event.Check(on)
        otherwin = self.splitter.wins[(self.num+1) % 2]
        if on:
            if otherwin.IsShown():
                self.enabled = True
            else:
                self.enabled = False

        event.Enable(self.enabled)


    def do_split(self, v):
        if v:
            if self.win.IsShown():
                return # assume split
            # split
            if self.splitter.IsSplit():
                # replace 
                self.splitter.ReplaceWindow(self.win, self.num)
            else:
                # split it -- get current win first
                curwin = self.splitter.GetWindow1()
                if self.num == 0:
                    wins = (self.win, curwin, self.splitter.sashpos)
                else:
                    wins = (curwin, self.win, self.splitter.sashpos)
                self.win.Show(True)
                self.splitter.split()
        else:
            self.splitter.Unsplit(self.win)
        

    
