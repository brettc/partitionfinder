'''Simplified Layouts for wx
'''

# TODO Make some deeper namespaces like pyx. It has a nice system.

import logging
log = logging.getLogger('newx.layout')

import string
import re
import wx

import policy
import optionset

# class LayoutError(RuntimeError): pass
# class AlignmentParseError(RuntimeError): pass
# class BorderParseError(RuntimeError): pass

class DisplayErrorCtrl(wx.StaticText):
    '''Used for showing errors in layout interactively'''

    def __init__(self, parent, cell):
        txt = 'ERROR: no child "' + cell.name + '"'
        style = wx.SIMPLE_BORDER|wx.ALIGN_CENTRE
        wx.StaticText.__init__(self, parent, wx.NewId(), txt, style=style)


# TODO need to make these special class that are directives so that they can
# be at the beginning of list


span_right, span_down = range(2)
# synonyms
colspan = span_right
rowspan = span_down

class _Cells(list):
    def __init__(self, cell=None):
        list.__init__(self)
        if cell:
            self.append(cell)

    def __or__(self, cell):
        self.append(cell)
        return self

    def append(self, cell):
        # This is where all cells get preprocessed
        if cell is span_right:
            cell = SpanRight()
        elif cell is span_down:
            cell = SpanDown()
        elif isinstance(cell, str):
            cell = Cell(cell)

        assert isinstance(cell, _CellBase)
        list.append(self, cell)

class Row(object):
    def __init__(self, cells, expanding=None):
        # Handle single cell case, or multiple cells
        self.expanding = expanding

        if isinstance(cells, _Cells):
            self.cells = cells
        else:
            self.cells = _Cells(cells)

    def __iter__(self):
        return iter(self.cells)

    def __len__(self):
        return len(self.cells)

    def __call__(self, name, **kwargs):
        # allow simple prototyping
        newrow = self.__class__(name, **self.__dict__)
        newrow.__dict__.update(kwargs)
        return newrow


class _CellBase(object):
    ROW, COL = 0, 1 # for tuples

    def __or__(self, other):
        cells = _Cells(self)
        return cells | other

    # def __add__(self, directive):
        # self._directive = directive
        # return self

    def __rshift__(self, expanding):
        # This actually goes on the row, we put it there later
        assert type(expanding) is int
        self._expanding = expanding
        return self

    def _bind(self, cells, i, j):
        self.pos = i, j
        self.span = [1, 1]
        cells[self.pos] = self

    col = property(lambda x: x.pos[x.COL])
    row = property(lambda x: x.pos[x.ROW])
    colspan = property(lambda x: x.span[x.COL])
    rowspan = property(lambda x: x.span[x.ROW])

# class _Directives(list):
    # def __init__(self, directive=None):
        # list.__init__(self)
        # if directive:
            # self.append(directive)

    # def __or__(self, directive):
        # self.append(directive)
        # return self

    # def append(self, directive):
        # assert isinstance(directive, _Directive)
        # list.append(self, directive)

# class _Directive(object):
    # def __or__(self, other):
        # directives = _Directives(self)
        # return directives | other

# class ColumnWidth(_Directive):
    # def __init__(self, expanding=1):
        # self.expanding = expanding

# class FixedWidth(ColumnWidth):
    # def __init__(self):
        # ColumnWidth.__init__(self, 0)

# class ColDef(object):
    # def __init__(self, directives):
        # if isinstance(directives, _Directives):
            # self.directives = directives
        # else:
            # self.directives = _Directives(directives)

class Span(_CellBase):
    def to_string(self):
        return ''

    def _find_span(self, cells, findtype):
        found = None
        lookat = list(self.pos)
        for i in range(self.pos[findtype]-1, -1, -1):
            try:
                lookat[findtype] = i
                found = cells[tuple(lookat)]
                break
            except:
                continue

        # we found something, could be cell or span
        if isinstance(found, Cell):
            found.span_to(self, findtype)
            self.cell = found
        elif isinstance(found, Span):
            found.cell.span_to(self, findtype)
        else:
            # TODO really need to check more than this
            raise BadLayoutError('No cell found for span')


class SpanRight(Span):
    def _bind(self, cells, i, j):
        Span._bind(self, cells, i, j)
        self._find_span(cells, self.COL)

class SpanDown(Span):
    def _bind(self, cells, i, j):
        Span._bind(self, cells, i, j)
        self._find_span(cells, self.ROW)


class Cell(_CellBase):
    '''Represents what is contained inside the cell'''

    def __init__(self, name=None, **kwargs):
        # FIXME -- need name AND label, cos name maps to proper name in python
        if name:
            self.name = name
        self.cls = kwargs.pop('cls', None)
        self.gap = kwargs.pop('gap', policy.layout.margin)

        # TODO: get policy defaults
        align = kwargs.pop('align', policy.layout.align)
        border = kwargs.pop('border', policy.layout.border)
        style = kwargs.pop('style', 0)

        if isinstance(align, str):
            align = optionset.alignment.parse(align)
                    
        if isinstance(border, str):
            border = optionset.border.parse(border)

        self.align = align
        self.border = border
        self.style = style

        # Anything left, we store for now.
        # TODO This is not really very cool
        self.__dict__.update(kwargs)

    def span_to(self, cell, findtype):
        self.span[findtype] = cell.pos[findtype] - self.pos[findtype] + 1

    def create_content(self, parent):
        # Blank names only allowed in prototypes
        # print self.__dict__
        assert hasattr(self, 'name')

        # NOTE -- override this in derived classes
        # TODO add some special names HLINE, VLINE
        content = getattr(parent, self.name, None)
        if not content and self.cls:
            # Try some different argument combinations
            ctrlid = wx.NewId()
            for args in [], [ctrlid]:
                try:
                    content = self.cls(parent, *args)
                    break
                except:
                    content = None

        return content


    def construct(self, gbs, parent):
        '''Construct the cell, using parent window, and place inside sizer
        '''
        content = self.create_content(parent)
        if not content:
            log.error('Cannot find child window named %s, inserting DUMMY', self.name)
            content = DisplayErrorCtrl(parent, self)

            # FIXME Is this right?
            setattr(parent, self.name, content)

        flags = self.border | self.align
        if self.style != 0:
            content.SetWindowStyle(self.style | ctrl.GetWindowStyle())

        # print content
        gbs.Add(content, self.pos, self.span, flags, self.gap)

        # TODO is this the best thing to do
        sz = content.GetSize()
        gbs.SetItemMinSize(content, sz.width, sz.height)

    def __call__(self, name, **kwargs):
        # allow simple prototyping
        newcell = self.__class__(name, **self.__dict__)
        newcell.__dict__.update(kwargs)
        return newcell


class Attributes(object):
    def __init__(self, **kw):
        pass

class Layout(object):
    def __init__(self, columndef, rows=[]):
        self._by_pos = {}
        self.rows = []
        self.columndef = columndef
        self.extend(rows)

    def __iter__(self):
        return iter(self.rows)

    def append(self, row):
        if not isinstance(row, tuple):
            row = row,
        arg_expanding = None
        arg_cells = _Cells()
        for bit in row:
            if isinstance(bit, int):
                arg_expanding = bit
            elif isinstance(bit, _Cells):
                arg_cells = bit
            elif isinstance(bit, Cell):
                arg_cells = _Cells(bit)
            else:
                # TODO generate error
                pass
        self.rows.append(Row(arg_cells, arg_expanding))

    def extend(self, rows):
        for r in rows:
            self.append(r)

    def _finalise(self):
        self._by_pos = {}
        for i, r in enumerate(self):
            for j, c in enumerate(r):
                c._bind(self._by_pos, i, j)

        # self.maxrow = len(self.rows)
        # self.maxcol = max([len(row) for row in self.rowrows])
        
    def _validate(self):
        # TODO make sure:
            # Nothing is repeated
            # There is no weird layout shite
            pass

    def make_sizer(self, parent):
        self._finalise()
        self._validate()

        gbs = wx.GridBagSizer()
        gbs.SetEmptyCellSize((0,0)) # lets us do lines and borders

        for rownum, row in enumerate(self):
            for colnum, cell in enumerate(row):
                if isinstance(cell, Cell):
                    cell.construct(gbs, parent)
            if row.expanding is not None:
                gbs.AddGrowableRow(rownum, row.expanding)

        for i, cd in enumerate(self.columndef.split('|')):
            try:
                expand = int(cd)
                gbs.AddGrowableCol(i, expand)
            except ValueError:
                # Not a number
                pass
        
        # Put it all in a box sizer, and put a margin around it
        margins = wx.BoxSizer(wx.VERTICAL)
        margins.Add(gbs, 1, wx.EXPAND | wx.ALL, policy.layout.border_margin)
        return margins


def make_sizer(parent, col, rows):
    return Layout(col, rows).make_sizer(parent)

class Panel(wx.Panel):
    '''Auto layout panel'''

    def __init__(self, parent, **kw):
        wx.Panel.__init__(parent, **kw)
        self.sizer = make_sizer(self.__wxlayout__)
        self.SetSizerAndFit(self.sizer)


class LabelledTextCell(Cell):
    # TESTING
    def __init__(self, name=None, labelalign='r', labelwidth=100, **kwargs):
        if isinstance(labelalign, str):
            labelalign = optionset.alignment.parse(labelalign)
        self.labelalign = labelalign
        self.labelwidth = labelwidth
        Cell.__init__(self, name, **kwargs)

    def create_content(self, parent):

        # TODO -- want to set labelwidth param
        content = wx.Panel(parent, wx.NewId())
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        content.label = wx.StaticText(content, wx.NewId(), self.name+':',
                                      style=self.labelalign | wx.ST_NO_AUTORESIZE)

        # lwidth = getattr(self, 'labelwidth', None)
        # if lwidth:
        w, h = content.label.GetSize()
        content.label.SetSize((self.labelwidth, h))
        content.label.SetMinSize((self.labelwidth, h))

        sizer.Add(content.label, 0, 
                  wx.ALIGN_CENTRE_VERTICAL | self.labelalign | wx.ALL | wx.RIGHT, 3)

        content.text = wx.TextCtrl(content, wx.NewId())  
        sizer.Add(content.text, 1, wx.EXPAND | wx.ALL, 0)
        content.SetAutoLayout(1)
        content.SetSizerAndFit(sizer)
        content.Layout()
        return content
        

def _testgui():

    class TestPanel(wx.Panel):

        textcell = Cell(cls=wx.TextCtrl)
        ltext = LabelledTextCell(labelwidth=50, labelalign='r')
        # cw = ColumnWidth(1)
        # rw = ExpandingRow(1)

        __wxlayout__ = Layout('1|2', (
            (textcell('text1')  | textcell('text2'), 1),
            (textcell('text3')  | span_right       , 2),
            (ltext('test')      | span_right       ),
            (ltext('blarg')     | span_right       ),
            # (span_down          | ltext('xx'))) TODO FIXME
            ))

            # TODO need to make these special class that are directives so that they can
            # be at the beginning of list


        def __init__(self, parent):
            wx.Panel.__init__(self, parent, id=wx.NewId())
            # self.text1 = idfree.TextCtrl(self)
            # self.text2 = idfree.TextCtrl(self)
            # self.text3 = idfree.TextCtrl(self)
            gbs = self.__wxlayout__.make_sizer(self)
            self.SetAutoLayout(1)
            self.SetSizerAndFit(gbs)
            self.Layout()

    
    class TestFrame(wx.Frame):

        def __init__(self, *args, **kwds):
            kwds["style"] = wx.DEFAULT_FRAME_STYLE
            wx.Frame.__init__(self, *args, **kwds)
            self.panel = TestPanel(self)


    class App(wx.App):
        # outputWindowClass = wxPyInformationalMessagesFrame
        def OnInit(self):
            frame = TestFrame(None, -1, "test")
            self.SetTopWindow(frame)
            frame.Show(True)
            return True


    # g = parse_format(format)
    # g._dump()
    # print g.cells

    app = App(0)
    app.MainLoop()


if __name__ == '__main__':
    _testgui()

