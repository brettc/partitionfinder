import wx
from wx.py import shell, editwindow
from newx import easy, layout

class ShellPanel(wx.Panel):

    # we might add some buttons
    #
    # __wxlayout__ = layoout.Grid(rows=(
        # layout.Row
                                # '''
    # P : 1
    # 1 : shell
    # '''

    MYFACES = editwindow.FACES.copy()
    # reset to my liking
    MYFACES['size'] = 9

    def __init__(self, parent, **kw):
        '''Send keywords to update the defaults'''
        wx.Panel.__init__(self, parent, wx.NewId())

        myfaces = self.MYFACES.copy()
        if kw:
            myfaces.update(kw)

        # this should be fixed in wx -- locals shouldn't be used!
        self.shell = shell.Shell(
            self, 
            wx.NewId(), 
            wx.DefaultPosition, 
            wx.DefaultSize,
            wx.CLIP_CHILDREN, 
            '', 
            {})

        self.shell.setStyles(myfaces)
        self.locals = self.shell.interp.locals

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.shell, 1, wx.EXPAND | wx.ALL, 3) 

        self.SetAutoLayout(1)
        self.SetSizerAndFit(sizer)
        self.Layout()
        # server.subscribe(topic=('world','new'), listener=self.update_locals)
        # server.subscribe(topic=('organism', 'new'), listener=self.update_locals)
        # self.Layout()

    # def update_locals(self, msg):
        # self.locals[msg.topic[0]] = msg.data

    # def OnPageClosed(self):
        # print 'leaving Shell page'
