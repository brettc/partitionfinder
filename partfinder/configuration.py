import logging
log = logging.getLogger("configuration")

import os

from pyparsing import (
    Word, Dict, OneOrMore, alphas, nums, Suppress, Optional, Group, stringEnd,
    delimitedList, pythonStyleComment, ParseException, line, lineno, col,
    Keyword)

class ConfigurationError(Exception):
    pass

class ParserError(Exception):
    def __init__(self, text, loc, msg):
        """Used for our own parsing problems"""
        self.line = line(loc, text)
        self.col = col(loc, text)
        self.lineno = lineno(loc, text)
        self.msg = msg

    def format_message(self):
        return "%s at line:%s, column:%s" % (self.msg, self.lineno, self.col)

class Parser(object):
    """Loads the Configuration files and validates them"""
    def __init__(self, config):
        # Create a configuration object. This is what we'll return 
        self.config = config

        # Some syntax that we need, but don't bother looking at
        SEMI = Suppress(";")
        EQUALS = Suppress("=")
        OPENB = Suppress("(")
        CLOSEB = Suppress(")")
        BACKSLASH = Suppress("\\")
        DASH = Suppress("-")

        ALIGN = Keyword("alignment")
        SCHEMES = Keyword("schemes")

        VARIABLE = Word(alphas + nums + "-_")
        VALUE = Word(alphas + nums + '_-/.\\')

        # General Section 
        general_def = (VARIABLE("name") + EQUALS +
                       VALUE("value")).setParseAction(self.def_variable_action) + Optional(SEMI)
        general = OneOrMore(general_def)

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

        self.config_parser = (
            general + 
            Suppress("[partitions]") + partlist + 
            Suppress("[schemes]") + schemelist + 
            stringEnd)

    def def_variable_action(self, text, loc, var_def):
        """Define a variable -- we'll check here if it is allowable"""
        if var_def.name == "alignment":
            self.config.alignment_file = var_def.value
        else:
            raise ParserError(text, loc, "'%s' is not an allowable setting" %
                                 var_def.name)

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
        if part_def.name in self.config.parts:
            raise ParserError(text, loc, "Repeated Partition Name '%s'" %
                                     part_def.name)

        # The actual partition definitions are defined as a list after the
        # first element (the first element is the name)
        log.debug("Found partition %s", part_def.name)
        self.config.parts[part_def.name] = part_def[1:]

    def check_part_exists(self, text, loc, partref):
        if partref.name not in self.config.parts:
            raise ConfigurationError(text, loc, "Partition %s not defined" %
                                     partref.name)

    def parse_file(self, fname):
        s = open(fname, 'r').read()
        try:
            self.parse_configuration(s)
        except ParserError, p:
            log.error(p.format_message())
            raise ConfigurationError()

    def parse_configuration(self, s):
        self.config_parser.ignore(pythonStyleComment).parseString(s)
        return self.config

class Configuration(object):
    """We use this to hold all of the configuration info"""
    def __init__(self, folder):
        self.folder = folder # Base directory
        self.config_file = os.path.join(self.folder, "partition_finder.cfg")

        if not os.path.exists(self.config_file) or not os.path.isfile(self.config_file):
            log.error("Configuration file '%s' does not exist",
                      self.config_file)
            raise ConfigurationError()
                
        self.parts = {}

    def load(self):
        p = Parser(self)
        log.debug("Loading configuration at '%s'", self.config_file)
        p.parse_file(self.config_file)


test1 = r"""

alignment = oeu

[partitions]
Gene1_pos1 = 1-789\3
Gene1_pos2 = 2-789\3
Gene1_pos3 = 3-789\3
Gene2_pos1 = 790-1449\3
Gene2_pos2 = 791-1449\3
Gene2_pos3 = 792-1449\3
Gene3_pos1 = 1450-2208\3
Gene3_pos2 = 1451-2208\3
Gene3_pos3 = 1452-2208\3

[schemes]
allsame			= 		(Gene1_pos1, Gene1_pos2, Gene1_pos3, Gene2_pos1, Gene2_pos2, Gene2_pos3, Gene3_pos1, Gene3_pos2, Gene3_pos3)
by_gene 		= 		(Gene1_pos1, Gene1_pos2, Gene1_pos3) (Gene2_pos1, Gene2_pos2, Gene2_pos3) (Gene3_pos1, Gene3_pos2, Gene3_pos3)
1_2_3	 		= 		(Gene1_pos1, Gene2_pos1, Gene3_pos1) (Gene1_pos2, Gene2_pos2, Gene3_pos2) (Gene1_pos3, Gene2_pos2, Gene3_pos3)
1_2_3_by_gene	= 		(Gene1_pos1) (Gene1_pos2) (Gene1_pos3) (Gene2_pos1) (Gene2_pos2) (Gene2_pos3) (Gene3_pos1) (Gene3_pos2) (Gene3_pos3)
12_3 			= 		(Gene1_pos1, Gene1_pos2, Gene2_pos1, Gene2_pos2, Gene3_pos1, Gene3_pos2) (Gene1_pos3, Gene2_pos3, Gene3_pos3)
12_3_by_gene 	= 		(Gene1_pos1, Gene1_pos2) (Gene1_pos3) (Gene2_pos1, Gene2_pos2) (Gene2_pos3) (Gene3_pos1, Gene3_pos2) (Gene3_pos3)
"""

if __name__ == '__main__':
    c = Parser()
    config = c.parse_configuration(test1)
    print config
    # print c.parts







