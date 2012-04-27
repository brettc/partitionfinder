import wx
def Notebook(parent, **kwargs):
    return wx.Notebook(parent, id=wx.NewId(), **kwargs)

def TextCtrl(parent, **kwargs):
    return wx.TextCtrl(parent, id=wx.NewId(), **kwargs)
