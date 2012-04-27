import wx

def getname(win):
    # TODO think of something better..?
    if isinstance(win, wx.Window):
        tit = win.GetTitle()
        if not tit:
            return win.__class__.__name__
    elif isinstance(win, wx.App):
        return win.GetAppName()
