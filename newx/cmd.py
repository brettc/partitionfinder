'''Commands as actual objects

Doing this simplifies much of the crap usually done with commands in a
windowing system. For example, you no longer need UPDATE_UI events as you can
modify the command directly in the registry and it will automatically update
the associated gui items.

'''
import logging
log = logging.getLogger('newx.cmd')
import policy
import wx
import string
from weakref import WeakKeyDictionary

# This is used in menus and toolbars...
SEPARATOR = '__a_placeholder_object__'
__ = SEPARATOR

class Command(object):
    '''A Command Object

    This captures everything about a particular command, include its text,
    associated bitmaps, and which gui objects are attached to it.
    '''

    def __init__(self, text, help='', longhelp='', keyboard='', bitmap=None):

        self.wxid = wx.NewId()
        self._enabled = True
        self._handlers = WeakKeyDictionary()
        self._gui_items = WeakKeyDictionary()

        self._text = text
        self._help = help
        self._longhelp = longhelp
        self._bitmap = bitmap

        # Is this modifiable?
        # It ain't implemented yet...
        self.keyboard = keyboard

    def set_text(self, text):
        self._text = text
        self.update()
    text = property(lambda x: x._text, set_text)

    def set_help(self, help):
        self._help = help
        self.update()
    help = property(lambda x: x._help, set_help)

    def set_enable(self, state=True):
        self._enabled = state
        if state: statetext = "En" 
        else: statetext = "Dis"
        log.debug("%sabling Command %s", statetext, self.text)
        self.update()
    enabled = property(lambda x: x._enabled, set_enable)

    def update(self):
        '''Update attached menus and toolbars'''
        for item in self._gui_items:
            if isinstance(item, wx.Menu):
                log.debug('Updating menu item for %s' % self.text)
                item.SetLabel(self.wxid, self.text)
                item.SetHelpString(self.wxid, self.help)
                item.Enable(self.wxid, self._enabled)

            elif isinstance(item, wx.ToolBar):
                # TODO Need to handle resetting of TextButton !
                log.debug('Updating toolbar item for %s' % self.text)
                ctrl = item.FindControl(self.wxid)
                if ctrl:
                    if isinstance(ctrl, wx.Button):
                        ctrl.SetLabel(self.text)
                        ctrl.Enable(self._enabled)
                    else:
                        log.error('cannot find gui button on toolbar for %s', self.cmd_name)
                else:
                    item.SetToolShortHelp(self.wxid, self.help)
                    item.SetToolLongHelp(self.wxid, self.longhelp)
                    item.EnableTool(self.wxid, self._enabled)

    def handle(self, event):
        '''Default handler for commands

        If any commands are unbound, they'll end up here
        This is an error, as stuff *should* be bound to something, either by
        binding to some window function or overriding this function
        '''
        log.error('Unhandled %s Command', self.cmd_name)

    def bind_button(self, handler, func):
        handler.Bind(wx.EVT_BUTTON, func, id=self.wxid) 

    def bind_menu(self, handler, func):
        handler.Bind(wx.EVT_MENU, func, id=self.wxid)

    def bind(self, handler):
        '''This binds this command into a handler (usually a window).
        '''

        # TODO Hmm. what happens if we do this twice?
        func = getattr(handler, self.function_name, None)
        if func:
            #FIXME - won't work with wxAPP -- as we need a name for window
            try:
                title = handler.GetTitle()
                log.debug('binding command %s to window %s', self.cmd_name,
                        title)
            except:
                pass

            # Not sure about this?
            self._handlers[handler] = 1
        elif not self._handlers:
            log.warning('binding command %s to default command handler',
                      self.cmd_name)
            # TODO should have default binding happening at the App level??
            # bind it to the handle function in the command class
            func = self.handle

        # Add both menu and button handlers, as we might be adding a text
        # button to the toolbar
        self.bind_menu(handler, func)
        self.bind_button(handler, func)

    def add_gui_item(self, item):
        # add a weakref into here
        self._gui_items[item] = None

class ToggleCmd(Command):
    '''A command that is Toggled

    This binds using a different event
    '''

    # TODO needs to update gui items

    def __init__(self, text, **kw):
        Command.__init__(self, text, **kw)
        self.state_on = False

    def bind_button(self, handler, func):
        handler.Bind(wx.EVT_TOGGLEBUTTON, func, id=self.wxid) 

    # defaults just use internal state
    def get_on(self):
        return self.state_on

    def set_on(self, v):
        if v: 
            self.state_on = True
        else: 
            self.state_on = False

    on = property(get_on, set_on)

    def toggle(self):
        self.on = not self.on

    def handle(self, event):
        self.toggle()

class CommandRegistry(object):
    '''The namespace for all commands ids'''

    def __getattr__(self, cmd_id):
        '''
        Any attempt to access a command id that does not exist will end up
        here. That is fine. We just warn about it, then insert it.
        '''
        log.warning('No Command name for %s, inserting one...' % cmd_id)

        # Call the command below to add it
        return add_single_cmd(cmd_id)

# Registry of all commands indexed available mnemonic name: FILE_OPEN etc
ids = CommandRegistry()

# TODO rename to this?
registry = ids

# Table of all commands indexed by the unique wxid
command_table = {}

def from_event(event):
    '''Return the associated command from an recieved event'''
    return command_table[event.GetId()]

def add(**kw):
    '''Add some commands to the registry
    '''
    # TODO toggle is done -- but how to handle groups??
    for cmd_name, stuff in kw.iteritems():
        add_single_cmd(cmd_name, stuff)

def add_single_cmd(cmd_name, stuff=None):
    if stuff is None:
        stuff = cmd_name

    # need to bypass __getattr__ here, otherwise we recurse..
    if hasattr(ids.__dict__, cmd_name):
        # TODO Wondering about stuff that has already been connected to
        # gui_items here, should we check ...?
        log.warning('%s is already a command; replacing it...' % cmd_name)

    if isinstance(stuff, Command):
        command = stuff
    else:
        # Check for keyboard stuff
        keyboard = ''
        if isinstance(stuff, str):
            if '\t' in stuff:
                stuff = stuff.split('\t')
                keyboard = stuff[1]
                stuff = stuff[0]

        if policy.naming.is_toggle_command(cmd_name):
            command = ToggleCmd(stuff, keyboard=keyboard)
        else:
            command = Command(stuff, keyboard=keyboard)

    # Add to our list of commands
    log.debug('Adding command %s' % cmd_name)
    setattr(ids, cmd_name, command)

    # Embed some stuff into the class
    command.cmd_name = cmd_name
    command.function_name = policy.naming.cmd_to_function_name(cmd_name)
    command.bitmap_name = policy.naming.cmd_to_bitmap_name(cmd_name)

    # lastly, add to global command table indexed by id
    command_table[command.wxid] = command
    return command

def bind(handler):
    '''Bind the commands into a handler implicitly by name'''
    for cid, com in command_table.iteritems():
        com.bind(handler)

# Accelerator stuff -- how should we handle this?
# if 0:
    # # This is another way to set Accelerators, in addition to
    # # using the '\t<key>' syntax in the menu items.
    # aTable = wx.AcceleratorTable([(wx.ACCEL_ALT,  ord('X'), exitID),
                                    # (wx.ACCEL_CTRL, ord('H'), helpID),
                                    # (wx.ACCEL_CTRL, ord('F'), findID),
                                    # (wx.ACCEL_NORMAL, WXK_F3, findnextID)
                                    # ])
    # self.SetAcceleratorTable(aTable)
    #

