'''Policies -- this controls preferences for such things as ...

    Naming: 

    GUI positioning

    NOTE: We shouldn't really have any imports from newx here, as this file
    should be able to be loaded and modified before any other stuff 
'''

import string

class DefaultNaming(object):
    def cmd_to_function_name(self, cmd_name):
        '''Given FILE_OPEN make OnFileOpen'''
        cmd_name = string.capwords(cmd_name, '_').replace('_', '')
        return 'On' + cmd_name

    def cmd_to_bitmap_name(self, cmd_name):
        '''Given FILE_OPEN make file_open'''
        return cmd_name.lower()

    def is_toggle_command(self, name):
        return name.lower().startswith('toggle')

class DefaultLayout:
    # TODO: Maybe do some sort of translation thingy from Font size?
    bitmap_size = 21, 18
    margin = 4
    vmargin = margin
    hmargin = margin

    # Around the edge of the controls in a layout
    border_margin = 6

    # Default values for Layouts
    align = 'x'
    border = 'a'
    gap = 4
    style = 0

naming = DefaultNaming()
layout = DefaultLayout()
