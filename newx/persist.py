
# from wx.evtmgr eventManager
import wx
from wx.lib.pubsub import Publisher as server
from wx.lib.evtmgr import eventManager
from weakref import WeakKeyDictionary, WeakValueDictionary

# TODO 
'''
Register a class
1. We detect what kind of class it is, must have a name...
2. We register the right kind of events for it, to record any changes
3. We write these somewhere

persist.register('name', window)

We do:
    Sizes, Position for Frame

Splits for splitters

Columns for ListCtrl

'''

_windows = WeakKeyDictionary()
_byname = WeakValueDictionary()
_data = {}

def OnColumnSize(event):
    assert isinstance(event, wx.ListEvent)
    # get the name
    win = event.GetEventObject()
    name = _windows[win]
    col = event.GetColumn()
    sz = win.GetColumnWidth(col)
    print 'saving col', name, col, sz
    print 'data is', event.GetData()

    print _data
    _data[('colsize', name, col)] = sz

def UpdateColumnSize(message):
    name = message.topic[1]
    win = _byname[name]
    assert isinstance(win, wx.ListCtrl)
    print 'got colw', message
    win.SetColumnWidth(int(message.topic[2]), message.data)

server.subscribe(topic='colsize', listener=UpdateColumnSize)

def register(name, win):
    # save the name
    _windows[win] = name
    _byname[name] = win
    print 'registering', name

    if isinstance(win, wx.ListCtrl):
        eventManager.Register(OnColumnSize, wx.EVT_LIST_COL_END_DRAG, win)


#
# publish persist.broadcast messages after startup ?
# before initial show?
#
        # a = ('sports')
        # b = ('sports','baseball')
        # a.matches(b) --> 1
        # b.matches(a) --> 0


# class Entry:
    # def __init__(self, name, dflt, parent, event, getv):
        # self.name = name
        # self.getv = getv
        # self.data = _data.setdefault(name, dflt)
        # parent.Bind(event, self.update)

    # def __call__(self):
        # return self.data

    # def update(self, event):
        # print 'updating'
        # self.data = self.getv()
        # _data[self.name] = self.data


_data = {}

def load(fname):
    global _data
    d = dict()
    try:
        execfile(fname, d)
    except:
        return

    _data =  d.setdefault('values', {})
    # print d['values']
    # print 'loaded', _data

def broadcast():
    # now broadcast
    for k, v in _data.copy().iteritems():
        server.sendMessage(topic=k, data=v)

import os

def save(fname):
    if os.path.exists(fname):
        os.remove(fname)

    out = open(fname, 'w')
    out.write('values = {\n')
    for name, value in _data.items():
        out.write('%s : %s\n' % (`name`, `value`))
    out.write('}\n')
    out.close()

_entries = {}

def set(name, value):
    pass
    # _data[name] = value


# def update():
    # for e in _entries.values():
        # e.update()
        # print e.data


    # Everyone's interested in politics...
    # for x in lList:
        # Publisher().subscribe(topic='politics', listener=x.notify)  # also tests singleton

    # # But only the first four are interested in trivia.
    # for x in lList[:4]:
        # server.subscribe(topic='trivia',   listener=x.notify)

    # # This one subscribes to everything.
    # everythingListener = SimpleListener(999)
    # server.subscribe(topic=(), listener=everythingListener.notify)

    # # Now send out two messages, testing topic matching.
    # server.sendMessage(topic='trivia',               data='What is the capitol of Oregon?')
    # server.sendMessage(topic=('politics','germany'), data='The Greens have picked up another seat in the Bundestag.')


