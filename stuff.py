
class Stuff(object):
    def __init__(self, ident):
        self.ident = ident

class Factory(object):
    def __index__(self):
        self.thestuff = {}

    def get_thing_by_identifier(self, ident_or_path_name):
        if ident in self.thestuff:
            return self.thestuff[ident]

        newstuff = Stuff(ident)
        self.thestuff[ident] = Stuff(ident)

        return newstuff



knob

