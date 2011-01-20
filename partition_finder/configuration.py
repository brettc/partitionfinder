from pyparsing import Word, Dict, OneOrMore, alphas, nums, \
        Suppress, Optional, Group, stringEnd, delimitedList, \
        pythonStyleComment, ParseException, line, lineno, col


test1 = r"""
[partitions]
Gene1_pos1 = 1-500 501-555\3
Gene1_pos2 = 2-789\3
Gene1_pos3 = 3-789\3
Gene2_pos1 = 790-1449\3
Gene2_pos2 = 791-1449\3
Gene2_pos3 = 792-1449\3
Gene3_pos1 = 1450-2208\3
Gene3_pos2 = 1451-2208\3
Gene3_pos3 = 1452-2208\3

has_a_duplicate = 100-200 200-400
			start_stop_problem    =   790-300;


Gene3_pos1_2_3 = 1450-2208


[schemes]
allsame			= 		(Gene1_pos1, Gene1_pos2, Gene1_pos3, Gene2_pos1, Gene2_pos2, Gene2_pos3, Gene3_pos1, Gene3_pos2, Gene3_pos3)
by_gene 		= 		(Gene1_pos1, Gene1_pos2, Gene1_pos3) (Gene2_pos1, Gene2_pos2, Gene2_pos3) (Gene3_pos1, Gene3_pos2, Gene3_pos3)
1_2_3	 		= 		(Gene1_pos1, Gene2_pos1, Gene3_pos1) (Gene1_pos2, Gene2_pos2, Gene3_pos2) (Gene1_pos3, Gene2_pos2, Gene3_pos3)
1_2_3_by_gene	= 		(Gene1_pos1) (Gene1_pos2) (Gene1_pos3) (Gene2_pos1) (Gene2_pos2) (Gene2_pos3) (Gene3_pos1) (Gene3_pos2) (Gene3_pos3)
12_3 			= 		(Gene1_pos1, Gene1_pos2, Gene2_pos1, Gene2_pos2, Gene3_pos1, Gene3_pos2) (Gene1_pos3, Gene2_pos3, Gene3_pos3)
12_3_by_gene 	= 		(Gene1_pos1, Gene1_pos2) (Gene1_pos3) (Gene2_pos1, Gene2_pos2) (Gene2_pos3) (Gene3_pos1, Gene3_pos2) (Gene3_pos3)
"""

class ConfigurationError(Exception):
    def __init__(self, text, loc, msg):
        """Used for our own parsing problems"""
        self.line = line(loc, text)
        self.col = col(loc, text)
        self.lineno = lineno(loc, text)
        self.msg = msg

    def format_message(self):
        return "%s at line:%s, column:%s" % (self.msg, self.lineno, self.col)
        

class Configuration(object):
    """Loads the Configuration files and validates them"""
    def __init__(self):
        super(Configuration, self).__init__()

        # Some syntax that we need, but don't bother looking at
        SEMI = Suppress(";")
        EQUALS = Suppress("=")
        OPENB = Suppress("(")
        CLOSEB = Suppress(")")
        BACKSLASH = Suppress("\\")
        DASH = Suppress("-")

        # Partition Parsing
        column = Word(nums)
        partname = Word(alphas + '_-' + nums)
        partdef = (column("start") + DASH + column("end") + 
                   Optional(BACKSLASH +
                            column("step"))).setParseAction(self.part_action)
        partition = (partname("name") + EQUALS +
                     OneOrMore(partdef)).setParseAction(self.part_def_action) + Optional(SEMI)
        partlist = OneOrMore(partition)

        # Scheme Parsing
        schemename = Word(alphas + '_-' + nums)
        partnameref = partname.copy()
        partnameref.setParseAction(self.check_part_exists)
        scheme = Group(OneOrMore(Group(OPENB +
                                       delimitedList(partnameref("name")) +
                                       CLOSEB)))
        schemedef = Group(schemename + EQUALS + scheme) + Optional(SEMI)
        schemelist = OneOrMore(schemedef)

        self.config_parser = (Suppress("[partitions]") + partlist + 
                Suppress("[schemes]") + schemelist + stringEnd)

        # Setup internal lists
        self.parts = {}

    def part_action(self, part):
        """Turn the 2 or 3 tokens into a tuple, supplying a default if needed"""
        fromc = int(part.start)
        toc = int(part.end)
        if part.step:
            stepc = int(part.step)
        else:
            stepc = 1
        return fromc, toc, stepc

    def part_def_action(self, text, loc, part_def):
        """We have everything we need here to make a partition"""
        if part_def.name in self.parts:
            raise ConfigurationError(text, loc, "Repeated Partition Name '%s'" %
                                     part_def.name)
           # raise ParseException("", loc, "Repeated Def")
        self.parts[part_def.name] = part_def[1:]

    def check_part_exists(self, text, loc, partref):
        if partref.name not in self.parts:
            raise ConfigurationError(text, loc, "Partition %s not defined" %
                                     partref.name)

    def parse_file(self, fname):
        s = open(fname, 'r').read()
        return self.parse_configuration(s)

    def parse_configuration(self, s):
        try:
            # We ignore the return as we build everything in the actions
            self.config_parser.ignore(pythonStyleComment).parseString(s)
        except ConfigurationError, c:
            print c.format_message()

if __name__ == '__main__':
    c = Configuration()
    c.parse_configuration(test1)
    print c.parts.keys()







