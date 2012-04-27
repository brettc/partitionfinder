import logging
import events

# import wx
# class wxLogger(wx.Log):
    # def OnLog(*args, **kwargs):
        # print kwargs

# wxLogger.SetActiveTarget(wxLogger())

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
            events.publish('newx.log.flush')

        # Hook it up to our event system
        events.publish('newx.log.message', record)

    def __len__(self):
        return self.current - self.base

    def __getitem__(self, i):
        i += self.base
        try:
            rec = self.records[i]
        except KeyError:
            raise IndexError # Change to this ..
        return rec

import sys
# FIXME output to stdout for now
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

gui_handler = KeepAroundHandler()
logging.getLogger().addHandler(gui_handler) # add at base level

if __name__  == '__main__':
    # A wee test

    def callback(info):
        print info.topic
        print info.data

    events.subscribe('newx', callback)
    log = logging.getLogger('newx.log.test')
    log.warning('eek')


