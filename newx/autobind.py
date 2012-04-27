import wx
import logging
log = logging.getLogger('newx.autobind')

def get_window_name(w):
    return "TODO"

# TODO -- add mappings for specific classes, and ONLY traverse them if they
# are appropriate (isinstance..)
#
# FIXME THIS is actually a better way of doing it
#
# Faster way of doing this -- assume that all On... functions are the ones
# that we want. First pass, get ALL children (wx.Window) then,
#
# Iterate through all of them, adding seeing if they match
# children , or are just window based...
#
# Should probably kept a check in the window of what is bound and what is not 
# ie self.boundfun = dict[self.fun] = 1 
#
# EG. wx.ListCtrl
# EVT_LIST_BEGIN_DRAG(id, func)  Begin dragging with the left mouse button.  
# EVT_LIST_BEGIN_RDRAG(id, func)  Begin dragging with the right mouse button.  
# EVT_LIST_BEGIN_LABEL_EDIT(id, func)  Begin editing a label. This can be prevented by calling Veto().  
# EVT_LIST_END_LABEL_EDIT(id, func)  Finish editing a label. This can be prevented by calling Veto().  
# EVT_LIST_DELETE_ITEM(id, func)  Delete an item.  
# EVT_LIST_DELETE_ALL_ITEMS(id, func)  Delete all items.  
# EVT_LIST_ITEM_SELECTED(id, func)  The item has been selected.  
# EVT_LIST_ITEM_DESELECTED(id, func)  The item has been deselected.  
# EVT_LIST_ITEM_ACTIVATED(id, func)  The item has been activated (ENTER or double click).  
# EVT_LIST_ITEM_FOCUSED(id, func)  The currently focused item has changed.  
# EVT_LIST_ITEM_MIDDLE_CLICK(id, func)  The middle mouse button has been clicked on an item.  
# EVT_LIST_ITEM_RIGHT_CLICK(id, func)  The right mouse button has been clicked on an item.  
# EVT_LIST_KEY_DOWN(id, func)  A key has been pressed.  
# EVT_LIST_INSERT_ITEM(id, func)  An item has been inserted.  
# EVT_LIST_COL_CLICK(id, func)  A column (m_col) has been left-clicked.  
# EVT_LIST_COL_RIGHT_CLICK(id, func)  A column (m_col) has been right-clicked.  
# EVT_LIST_COL_BEGIN_DRAG(id, func)  The user started resizing a column - can be vetoed.  
# EVT_LIST_COL_DRAGGING(id, func)  The divider between columns is being dragged.  
# EVT_LIST_COL_END_DRAG(id, func)  A column has been resized by the user.  
# EVT_LIST_CACHE_HINT(id, func)  Prepare cache f 

_mapped = {}

def addmapping(**kw):
    global _mapped
    for fn, evt in kw.iteritems():
        _mapped[evt] = fn

addmapping(
    Paint=wx.EVT_PAINT,
    LeftDown=wx.EVT_LEFT_DOWN,
    LeftUp=wx.EVT_LEFT_UP,
    RightDown=wx.EVT_RIGHT_DOWN,
    RightUp=wx.EVT_RIGHT_UP,
    Size=wx.EVT_SIZE,
    Motion=wx.EVT_MOTION,
    EraseBackground=wx.EVT_ERASE_BACKGROUND,
    KeyPressed=wx.EVT_KEY_DOWN,
    Idle=wx.EVT_IDLE,
    Timer=wx.EVT_TIMER,
    SetCursor=wx.EVT_SET_CURSOR,
)

_childmap = []
def addchildmapping(cls, **kw):
    _childmap.append((cls, kw))

addchildmapping(wx.Button,
    Pressed=wx.EVT_BUTTON
    )
addchildmapping(wx.TextCtrl,
    Updated=wx.EVT_TEXT
    )
addchildmapping(wx.ListCtrl,
    Selected=wx.EVT_LIST_ITEM_SELECTED,
    Deselected=wx.EVT_LIST_ITEM_DESELECTED,
    Focused=wx.EVT_LIST_ITEM_FOCUSED,
    )
addchildmapping(wx.ListBox,
    Clicked=wx.EVT_LISTBOX,
    DblClicked=wx.EVT_LISTBOX_DCLICK
    )
addchildmapping(wx.VListBox,
    Clicked=wx.EVT_LISTBOX,
    DblClicked=wx.EVT_LISTBOX_DCLICK
    )
addchildmapping(wx.Notebook,
    PageChanged=wx.EVT_NOTEBOOK_PAGE_CHANGED,
    PageChanging=wx.EVT_NOTEBOOK_PAGE_CHANGING
    )

def bind(window, toponly=True):
    for evt, name in _mapped.iteritems():
        name = 'On' + name
        # by default bind only top level, otherwise we rebind parent funcs
        if toponly and name not in window.__class__.__dict__:
            continue
        fn = getattr(window, name, None)
        if fn:
            log.debug('binding %s', name)
            # log.debug('binding %s in window %s', 
                      # name, 
                      # get_window_name(window))
            window.Bind(evt, fn)

    # Do specific windows too
    for cls, mapping in _childmap:
        if isinstance(window, cls):
            for nm, bind in mapping.items():
                possname = 'On' + nm
                fn = getattr(window, possname, None)
                if fn:
                    log.debug('binding %s', name)
                    # log.debug('binding %s in window %s', 
                              # possname,
                              # get_window_name(window))
                    window.Bind(bind, fn)

def bindchild(parent, child, prefix, toponly=True):
    for evt, name in _mapped.iteritems():
        # As above
        name = 'On' + prefix + name
        if toponly and name not in parent.__class__.__dict__:
            continue
        fn = getattr(parent, name, None)
        if fn:
            child.Bind(evt, fn)


# TODO make into a dictionary!! will need to look for __mro__ or something
def findchildmapping(inst):
    for cls, mapping in _childmap:
        if isinstance(inst, cls):
            return mapping
    return None

def bindcontrols(parent):
    # not the same as the others
    ctrls = {}

    for k, v in parent.__dict__.iteritems():
        if isinstance(v, wx.Control):
            # print 'found ', v
            # see if there are any functions the same
            lookup = findchildmapping(v)
            if lookup:
                for nm, bind in lookup.items():
                    possname = 'On' + k + nm
                    fn = getattr(parent, possname, None)
                    if fn:
                        log.debug('binding %s', possname)
                        # log.debug('binding %s in window %s',
                                  # possname,
                                  # get_window_name(parent))
                        parent.Bind(bind, fn, v)

