import wx
import autobind

# Extra pixels that each column takes up. May have to change for different
# platforms
EXTRA_COLUMN_WIDTH = 12
BETWEEN_COLUMN_WIDTH = 2

class AutoVirtualListCtrl(wx.ListCtrl):
    '''Automatic creation of columns and data based on special methods

    NOTE: You need member called self.items, this is a sequence, which
    '''

    def __init__(self, parent, autowidth=True, **kw):
        more_style = kw.pop('style', 0)
        wx.ListCtrl.__init__(self, parent, wx.NewId(),
            style = more_style | wx.LC_REPORT | wx.LC_VIRTUAL)

        fnt = self.GetFont()
        dc = wx.WindowDC(self)
        dc.SetFont(fnt)

        self.get_functions = []
        for colnum, coldef, in enumerate(self.__wxcolumns__):
            coldef.calc_width(dc)
            self.InsertColumn(colnum, coldef.heading, coldef.align,
                              width=coldef.width)
            self.get_functions.append(coldef.getfun)

        # Resize ourselves based on all of the columns.
        if autowidth:
            totalw = wx.SystemSettings.GetMetric(wx.SYS_VSCROLL_X)
            for i in range(self.GetColumnCount()):
                totalw += self.GetColumnWidth(i) + BETWEEN_COLUMN_WIDTH
            sz = self.GetSize()
            sz.width = totalw 
            # debug('Resizing %s to size %s', util.getname(self), totalw)
            self.SetSize(sz)
            self.SetSizeHints(totalw, -1)

        autobind.bind(self)


    def RedrawAll(self):
        self.SetItemCount(len(self.items))
        # self.RefreshItems(0, len(self.items)-1)
        self.Refresh()

    #---------------------------------------------------
    # These methods are callbacks for implementing the
    # "virtualness" of the list...
    def OnGetItemText(self, item_num, col_num):
        item = self.items[item_num]
        return self.get_functions[col_num](self, item)

    # TODO: Add images to our column definitions etc.
    # def OnGetItemImage(self, item):
        # return 0

    # def OnGetItemAttr(self, item):
        # return None
        #

LIST_ALIGN = { 
    'l' : wx.LIST_FORMAT_LEFT,
    'r' : wx.LIST_FORMAT_RIGHT,
    'c' : wx.LIST_FORMAT_CENTER,
    }

class Column(object):
    def __init__(self, getfun, heading, align='l', width=100):
        self.heading = heading
        self.width = width
        self.getfun = getfun

        if align not in LIST_ALIGN.values():
            try:
                align = LIST_ALIGN[align[0:1].lower()]
            except KeyError:
                align = LIST_ALIGN['l']

        self.align = align

    def calc_width(self, dc):
        # print self.heading, 
        if type(self.width) is not int:
            # measure it as text
            text = str(self.width)
            self.width, h = dc.GetTextExtent(text) 
            self.width += EXTRA_COLUMN_WIDTH
            # debug('measuring the label "%s", has size %s', text, width)

def _test():

    class TestListCtrl(AutoVirtualListCtrl):

        items = range(10000)
        bigstring = 'This is a long string of course'

        def get_test(self, num):
            return self.bigstring

        def get_num(self, num):
            return str(num)

        def get_x(self, num):
            return num * 2

        __wxcolumns__ = [
            Column(get_test, 'Test', width=bigstring),
            Column(get_num, 'Number', width='Number', align='r'),
            Column(get_x, 'Blarg', width='xxxxx'),
            ]

    class MainWindow(wx.Frame):
        def __init__(self, parent, id, title):
            wx.Frame.__init__(self, parent, -1, title, size = (500, 500),
                            style=wx.DEFAULT_FRAME_STYLE|wx.NO_FULL_REPAINT_ON_RESIZE)

            self.vlist = TestListCtrl(self)
            self.vlist.RedrawAll()
            self.Show(True)

    
    class App(wx.App):
        def OnInit(self):
            frame = MainWindow(None, -1, "Testing")
            self.SetTopWindow(frame)
            return True

    app = App(0)
    app.MainLoop()


if __name__ == '__main__':
    _test()



