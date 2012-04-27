import wx
import logging, autobind
import colorsys
from math import sin, cos, radians


log = logging.getLogger('newx.drawing')

class DrawPanel(wx.Panel):

    def __init__(self, parent):
        wx.Panel.__init__(self, parent, -1, style=wx.SUNKEN_BORDER)
        self.Bind(wx.EVT_PAINT, self.OnPaint)
        self.SetCursor(wx.StockCursor(wx.CURSOR_ARROW))
        self.font = wx.SystemSettings.GetFont(wx.SYS_DEFAULT_GUI_FONT)

    def Center(self, gc, x, y, text):
        w, h = gc.GetTextExtent(text)
        gc.DrawText(text, x-round(w/2.0), y-round(h/2.0))

    def OnPaint(self, evt):
        dc = wx.PaintDC(self)
        dc.Clear()
        gc = wx.GraphicsContext.Create(dc)
        gc.SetFont(self.font)
        self.width, self.height = self.GetSizeTuple()
        self.OnDraw(gc)

    def OnDraw(self, gc):
        gc.PushState()
        gc.SetPen(wx.Pen(wx.Colour(200,192,50,255)))
        gc.DrawRoundedRectangle(0, 0, -100, 100, 5)
        gc.PopState()

    def OnEraseBackground(self, evt):
        return True

if __name__ == '__main__':

    BASE  = 80.0    # sizes used in shapes drawn below
    BASE2 = BASE/2
    BASE4 = BASE/4
    
    class TestDrawPanel(DrawPanel):
        def OnDraw(self, gc):
            font = wx.SystemSettings.GetFont(wx.SYS_DEFAULT_GUI_FONT)
            font.SetWeight(wx.BOLD)
            gc.SetFont(font)

            # make a path that contains a circle and some lines, centered at 0,0
            path = gc.CreatePath()
            path.AddCircle(0, 0, BASE2)
            path.MoveToPoint(0, -BASE2)
            path.AddLineToPoint(0, BASE2)
            path.MoveToPoint(-BASE2, 0)
            path.AddLineToPoint(BASE2, 0)
            path.CloseSubpath()
            path.AddRectangle(-BASE4, -BASE4/2, BASE2, BASE4)


            # Now use that path to demonstrate various capbilites of the grpahics context
            gc.PushState()             # save current translation/scale/other state 
            gc.Translate(60, 75)       # reposition the context origin

            gc.SetPen(wx.Pen("navy", 1))
            gc.SetBrush(wx.Brush("pink"))

            # show the difference between stroking, filling and drawing
            for label, PathFunc in [("StrokePath", gc.StrokePath),
                                    ("FillPath",   gc.FillPath),
                                    ("DrawPath",   gc.DrawPath)]:
                if "wxGTK" in wx.PlatformInfo:
                    w, h = gc.GetTextExtent(label) # NYI in Cairo context
                else:
                    w, h = gc.GetTextExtent(label)

                gc.DrawText(label, -w/2, -BASE2-h)
                PathFunc(path)
                gc.Translate(2*BASE, 0)

                
            gc.PopState()              # restore saved state
            gc.PushState()             # save it again
            gc.Translate(60, 200)      # offset to the lower part of the window
            
            gc.DrawText("Scale", 0, -BASE2)
            gc.Translate(0, 20)

            # for testing clipping
            #gc.Clip(0, 0, 100, 100)
            #rgn = wx.RegionFromPoints([ (0,0), (75,0), (75,25,), (100, 25),
            #                            (100,100), (0,100), (0,0)  ])
            #gc.ClipRegion(rgn)
            #gc.ResetClip()
            
            gc.SetBrush(wx.Brush(wx.Colour(178,  34,  34, 128)))   # 128 == half transparent
            for cnt in range(8):
                gc.Scale(1.08, 1.08)    # increase scale by 8%
                gc.Translate(5,5)     
                gc.DrawPath(path)


            gc.PopState()              # restore saved state
            gc.PushState()             # save it again
            gc.Translate(400, 200)
            gc.DrawText("Rotate", 0, -BASE2)
            
            gc.Translate(0, 75)
            for angle in range(0, 360, 30):
                gc.PushState()         # save this new current state so we can pop back to 
                                       # it at the end of the loop
                r, g, b = [int(c * 255) for c in colorsys.hsv_to_rgb(float(angle)/360, 1, 1)]
                gc.SetBrush(wx.Brush(wx.Colour(r, g, b, 64)))

                # use translate to artfully reposition each drawn path
                gc.Translate(1.5 * BASE2 * cos(radians(angle)),
                             1.5 * BASE2 * sin(radians(angle)))

                # use Rotate to rotate the path
                gc.Rotate(radians(angle))

                # now draw it
                gc.DrawPath(path)
                gc.PopState()

            gc.PopState()

            gc.Translate(100, 100)
            gc.SetPen(wx.Pen(wx.Colour(200,192,50,255)))
            gc.DrawRoundedRectangle(-100, 100, -100, 100, 5)
            gc.PopState()


    class MainWindow(wx.Frame):
        def __init__(self, parent, id, title):
            wx.Frame.__init__(self, parent, -1, title, size = (500, 500),
                            style=wx.DEFAULT_FRAME_STYLE|wx.NO_FULL_REPAINT_ON_RESIZE)

            self.panel = TestDrawPanel(self)
            self.Show(True)

    
    class App(wx.App):
        def OnInit(self):
            frame = MainWindow(None, -1, "Testing")
            self.SetTopWindow(frame)
            return True

    app = App(0)
    app.MainLoop()


