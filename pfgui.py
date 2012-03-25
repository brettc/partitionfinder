import pickle
from GUI import Application, View, Document, Window, Cursor, rgb
from GUI import Window, Button, CheckBox, TextField, Label
from GUI.Files import FileType
from GUI.Geometry import pt_in_rect, offset_rect, rects_intersect
from GUI.StdColors import black, red

from partfinder import config, analysis_method, util, parser, reporter

class PFApp(Application):

    def __init__(self):
        Application.__init__(self)
        self.cfg_type = FileType(
            name = "PartitionFinder Document", 
            suffix = "cfg", 
            #mac_creator = "BLBE", mac_type = "BLOB", # These are optional
        )
        self.file_type = self.cfg_type
        # self. = Cursor("blob.tiff")
    
    def open_app(self):
        pass
        # self.new_cmd()

    def make_document(self, fileref):
        return PFDoc(file_type = self.cfg_type)

    def make_window(self, document):
        win = Window(size = (400, 400), document = document)
        win.shrink_wrap()

        view = PFView(model=document)
        win.place(view, left = 0, top = 0, right = 0, bottom = 0, sticky = 'nsew')
        win.show()

class PFView(View):

    def __init__(self, *args, **kwargs):
        View.__init__(self, *args, **kwargs)
        self.button = Button("Bing!", action = self.bing)
        self.label = Label(self.get_model().cfg.base_path)
        self.place(self.button, left=20, top=20)
        self.place(self.label, left=20, top=40)

    def bing(self):
        self.get_model().analyse()


    # def model_added(self, model):
        # self.label.text = model.cfg.base_path
    # def draw

class PFDoc(Document):

    def __init__(self, *args, **kwargs):
        Document.__init__(self, *args, **kwargs)
        self.cfg = config.Configuration()

    def new_contents(self):
        pass
        # self.blobs = []

    def read_contents(self, file):
        file.close()
        self.cfg.load_base_path(self.file.dir.get_path())

    # def read_contents(self, file):
        # print self.file
        # self.blobs = pickle.load(file)

    def write_contents(self, file):
        pass
        # pickle.dump(self.blobs, file)

    def analyse(self):
        # Now try processing everything....
        method = analysis_method.choose_method(self.cfg.search)
        rpt = reporter.TextReporter(self.cfg)
        anal = method(self.cfg, rpt, False, False)
        anal.analyse()

PFApp().run()
