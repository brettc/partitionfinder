import wx
import cmd
import logging
import policy
log = logging.getLogger('newx.toolbar')

BITMAPS = {}

def add_bitmaps(**kw):
    for n, bm in kw.items():
        if not isinstance(n, str):
            error('%s is not a string' % s)
            continue
        if not isinstance(bm, wx.Bitmap):
            error('%s is not a bitmap' % s)
            continue

        BITMAPS[n] = bm

# def add_error(toolbar, txt):
    # bmp = make_bitmap(txt, bg=wx.RED, fg=wx.BLACK)
    # # add it anyway
    # toolbar.AddSimpleTool(-1, bmp, "ERROR, no command for %s" % txt, '')
    # error('No command defined for toolbar item %s' % txt)

def make_bitmap(text, fg=wx.WHITE, bg=wx.BLACK):
    """Make a bitmap with bg color and text in it """
    bw, bh = policy.layout.bitmap_size

    dc = wx.MemoryDC()
    dc.SetFont(wx.SMALL_FONT)
    w, h = dc.GetTextExtent(text)

    bmp = wx.EmptyBitmap(bw, bh)
    dc.SelectObject(bmp)

    dc.SetBackground(wx.Brush(bg))
    dc.Clear()

    dc.SetTextForeground(fg)
    # centre it
    x = bw/2 - w/2
    y = bh/2 - h/2
    # dc.DrawText(text, (x, y))

    dc.SelectObject(wx.NullBitmap)

    return bmp

def make_toolbar(win, cmds):
    toolbar = win.CreateToolBar(wx.TB_HORIZONTAL|wx.NO_BORDER|wx.TB_FLAT)
    toolbar.SetToolBitmapSize(policy.layout.bitmap_size)
    sz = toolbar.GetToolBitmapSize()

    for com in cmds:
        if com is cmd.SEPARATOR:
            toolbar.AddSeparator()
            continue

        # FIXME For now, no bitmaps 
        # Need to find a nice way of loading these ...
        # bmp = BITMAPS.get(com.bmname, None)
        # bmp = make_bitmap(com.idname)
        bmp = None

        log.debug('Adding %s to toolbar', com.text)
        if bmp:
            if com.istoggle:
                toolbar.AddCheckLabelTool(com.wxid, com.text, bmp)
            else:
                toolbar.AddSimpleTool(com.wxid, bmp, com.shorthelp, com.longhelp)
        else:
            # warning('no bitmap name "%s", making text button', com.bmname)

            if isinstance(com, cmd.ToggleCmd):
                but = wx.ToggleButton(toolbar, com.wxid, com.text)
            else:
                but = wx.Button(toolbar, com.wxid, com.text)

            toolbar.AddControl(but)

        com.add_gui_item(toolbar)

    return toolbar

