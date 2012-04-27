import wx
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as Canvas
# from matplotlib.backends.backend_wx import FigureCanvasWx as Canvas
from matplotlib.figure import Figure
# from matplotlib.collections import LineCollection
# from matplotlib.patches import Polygon, Circle, Rectangle
# from matplotlib.lines import Line2D

from pylab import figure, close, axes, subplot, show, savefig
import pylab
from numpy import arange, sin, pi
from matplotlib import _pylab_helpers

class MatplotPanel(Canvas):
    def __init__(self, parent):
        self.figure = Figure()
        Canvas.__init__(self, parent, wx.NewId(), self.figure)

        ax1 = self.figure.add_subplot(111)
        t = arange(0.0, 1.0, 0.01)
        ax1.plot(t, sin(2*pi*t))
        ax1.grid(True)
        ax1.set_title('A sine wave or two')
        # ax1.set_ylim( (-2,2) )
        # ax1.set_ylabel('1 Hz')

    def make_active(self):
        _pylab_helpers.Gcf.set_active(self.figure)



