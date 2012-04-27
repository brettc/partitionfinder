'''Translation of options from strings to wx types

This simplifies the use of the wx options. I can never remember them, and
combining them makes things to BIG. We want concise.
'''

import logging
log = logging.getLogger('newx.optionset')
import wx

class _options(object):
    def __init__(self, name, optset):
        self.name = name
        self.xlate = optset
        self.valid = set(optset.values())

    def parse(self, opts):
        ret = 0
        for c in opts.lower():
            try:
                ret |= self.xlate[c]
            except KeyError:
                log.error('Unknown option translation %s in class %s', 
                          c, name)

        return ret

wx_alignments = {
    'c' : wx.ALIGN_CENTRE,
    'l' : wx.ALIGN_LEFT,
    't' : wx.ALIGN_TOP,
    'r' : wx.ALIGN_RIGHT,
    'b' : wx.ALIGN_BOTTOM,
    'v' : wx.ALIGN_CENTRE_VERTICAL,
    'x' : wx.EXPAND,
    's' : wx.SHAPED,
    }

wx_border = {
    't' : wx.TOP,
    'b' : wx.BOTTOM,
    'l' : wx.LEFT,
    'r' : wx.RIGHT,
    'a' : wx.ALL,
    }

wx_style = {
    's' : wx.SUNKEN_BORDER,
    'n' : wx.NO_BORDER,
    }

alignment = _options('alignment', wx_alignments)
border = _options('border', wx_border)
style = _options('style', wx_style)

