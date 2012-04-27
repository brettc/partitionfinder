import logging
log = logging.getLogger('newx.menu')
import string
import wx
import cmd

class InvalidMenuItem(Exception):
    pass

class Menu(object):
    '''A newx Menu object

    Simplifies the creation of the menus (a lot)
    '''

    def __init__(self, text, *args):
        self.wxid = wx.NewId()
        self.text = text
        self.entries = args
               
    def construct(self, parent):
        '''Construct a menu item

        This relies on the newx.cmd framework being in place
        '''
        wxmenu = wx.Menu()
        if isinstance(parent, wx.MenuBar):
            parent.Append(wxmenu, self.text)
        else:
            parent.AppendMenu(self.wxid, self.text, wxmenu)

        for entry in self.entries:
            if isinstance(entry, cmd.Command):
                if isinstance(entry, cmd.ToggleCmd):
                    flag = wx.ITEM_CHECK
                else:
                    flag = wx.ITEM_NORMAL

                if entry.keyboard:
                    txt = '\t'.join((entry.text, entry.keyboard))
                else:
                    txt = entry.text
                wxmenu.Append(entry.wxid, txt, entry.help, flag)

                # Record here that the menu is attached to the command.
                # Any further updates to the command will automatically update
                # the gui item
                entry.add_gui_item(wxmenu)

            elif isinstance(entry, Menu):
                entry.construct(self)
            elif entry is cmd.SEPARATOR:
                wxmenu.AppendSeparator()
            else:
                raise InvalidMenuItem

        return wxmenu


def make_bar(menus):
    bar = wx.MenuBar()
    for menu in menus:
        assert isinstance(menu, Menu)
        menu.construct(bar)

    return bar

def make_popup(menu):
    assert isinstance(menu, Menu)
    return menu_layout.construct()


def _test():
    import cmd
    cmd.add(
        OPEN=cmd.Command('Open up',help='Whatever'),
        CLOSE='Close',
        DO_SOMETHING='Do Something',
        TOGGLE=cmd.ToggleCmd('Do Something'),
        EXIT='Exit Out',
        )


    class MainWindow(wx.Frame):
        def __init__(self, parent, id, title):
            wx.Frame.__init__(self, parent, -1, title, size = (500, 500),
                            style=wx.DEFAULT_FRAME_STYLE|wx.NO_FULL_REPAINT_ON_RESIZE)

            menu = (Menu('&File', 
                         cmd.ids.OPEN, 
                         cmd.ids.CLOSE, 
                         cmd.ids.EXIT
                         ),
                    Menu('&Stuff', 
                         cmd.ids.DO_SOMETHING,
                         cmd.ids.TOGGLE
                         ),
                    )

            self.SetMenuBar(make_bar(menu))
            self.Show(True)
            cmd.bind(self)

        def OnToggle(self, event):
            print 'xxx'

        def OnDoSomething(self, event):
            dlg = wx.MessageDialog (self, 'wooooh', 'Error!', wx.OK)
            # cmd.ids.OPEN.text = 'WAzza'
            dlg.ShowModal()

        def OnExit(self, event):
            self.Close()

    
    class App(wx.App):
        def OnInit(self):
            frame = MainWindow(None, -1, "Testing")
            # dc = wx.WindowDC(frame)
            # print dc.GetHDC()
            self.SetTopWindow(frame)
            return True

    app = App(0)
    app.MainLoop()


if __name__ == '__main__':
    _test()



