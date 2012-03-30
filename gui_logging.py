import logging
import wx
import wx.lib.evtmgr as events
from wx.lib.pubsub import Publisher

publisher = Publisher()

def publish(path, data=None):
    assert isinstance(path, str)
    topic = tuple(path.split('.'))
    # print topic, data
    publisher.sendMessage(topic, data)

def subscribe(path, listener):
    assert isinstance(path, str)
    topic = tuple(path.split('.'))
    # print listener, topic
    publisher.subscribe(listener, topic)

def unsubscribe(path, listener):
    assert isinstance(path, str)
    topic = tuple(path.split('.'))
    publisher.unsubscribe(listener, topic)

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
    
# import wx
# class wxLogger(wx.Log):
    # def OnLog(*args, **kwargs):
        # print kwargs

# wxLogger.SetActiveTarget(wxLogger())

class LogHtmlListBox(wx.HtmlListBox):
    def __init__(self, parent, **kwargs):
        wx.HtmlListBox.__init__(self, parent, **kwargs)
        # gui_handler.add_watch(self)
        # autobind.bind(self)
        self.Update()
        subscribe('newx.log.message', self.Update)
    
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
        dc.SetPen(wx.LIGHT_GREY_PEN)
        dc.DrawLine(rect.left, rect.bottom, rect.right, rect.bottom)

    # def OnDblClicked(self, evt):
        # We are going to extract the file name and line number and send them
        # out over the events system. The idea is that you can pick this up, 
        # and then send the user to the offending line in an editor.

        # n = evt.GetInt()
        # rec = gui_handler[n]
# import events

class LogListCtrl(wx.ListCtrl):
    '''For showing incoming log messages'''

    COLUMNS = 'name levelname'.split()

    def __init__(self, parent):
        global gui_handler
        wx.ListCtrl.__init__(self, parent, wx.NewId(), size=(300,200),
                style = wx.LC_REPORT | wx.LC_VIRTUAL)

        self.SetItemCount(len(gui_handler))

        self.lookup = {}
        for c, nm in enumerate(self.COLUMNS):
            self.InsertColumn(c, nm, wx.LIST_FORMAT_LEFT)
            self.lookup[c] = nm

        # gui_handler.add_watch(self)
        subscribe('newx.log.message', self.Update)

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

class KeepAroundHandler(logging.Handler):
    """
    A handler class which sends keeps around the last N records for active
    display in a GUI
    """

    def __init__(self, limit=400, remove=50):
        '''When there are 'limit' records, it will flush 'remove' number of
        them'''
        logging.Handler.__init__(self)
        self.limit = limit
        self.remove = remove
        self.base = 0
        self.current = 0
        self.records = {}


    def emit(self, record):
        '''Keep it around, flush 'remove' records if 'limit' is reached'''

        self.records[self.current] = record
        self.current += 1
        if len(self.records) == self.limit:
            newbase = self.base + self.remove
            for j in range(self.base, newbase):
                del self.records[j]

            self.base = newbase
            publish('newx.log.flush')

        # Hook it up to our event system
        publish('newx.log.message', record)

    def __len__(self):
        return self.current - self.base

    def __getitem__(self, i):
        i += self.base
        try:
            rec = self.records[i]
        except KeyError:
            raise IndexError # Change to this ..
        return rec

gui_handler = KeepAroundHandler()
logging.getLogger().addHandler(gui_handler) # add at base level

