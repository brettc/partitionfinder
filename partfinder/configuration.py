import logging
log = logging.getLogger("config")

import os

from pyparsing import (
    Word, Dict, OneOrMore, alphas, nums, Suppress, Optional, Group, stringEnd,
    delimitedList, pythonStyleComment, ParseException, line, lineno, col,
    Keyword)

class ConfigurationError(Exception):
    pass

class ParserError(Exception):
    """Used for our own parsing problems"""
    def __init__(self, text, loc, msg):
        self.line = line(loc, text)
        self.col = col(loc, text)
        self.lineno = lineno(loc, text)
        self.msg = msg


    def format_message(self):
        return "%s at line:%s, column:%s" % (self.msg, self.lineno, self.col)

class Parser(object):
    """Parse configuration files

    The results are put into the configuration object
    """
    def __init__(self, config):
        self.init_config(config)
        self.init_grammar()

    def init_config(self, config):
        self.config = config
        config.parts = {}

    def init_grammar(self):
        """Set up the parsing classes
        Any changes to the grammar of the config file be done here.
        """
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
        # Just the assignment of variables
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

        # We've defined the grammar for each section. Here we just put it all
        # together
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
            raise ParserError(text, loc, "Partition %s not defined" %
                                     partref.name)

    def parse_file(self, fname):
        s = open(fname, 'r').read()
        try:
            self.parse_configuration(s)
        except ParserError, p:
            # Catch any parsing errors, print something out, then raise a
            # configuration error
            log.error(p.format_message())
            raise ConfigurationError

    def parse_configuration(self, s):
        self.config_parser.ignore(pythonStyleComment).parseString(s)

class Configuration(object):
    """We use this to hold all of the configuration info"""
    def __init__(self, base_path):
        self.base_path = os.path.abspath(base_path)
        if not os.path.exists(self.base_path) or \
           not os.path.isdir(self.base_path):
            log.error("No such folder: '%s'", self.base_path)
            raise ConfigurationError

        self.config_path = os.path.join(self.base_path, "partition_finder.cfg")
        self.output_path = os.path.join(self.base_path, "output")
        self.log_path = os.path.join(self.base_path, "partition_finder.log")

        # Add a log file in this folder
        self.init_log()

        # Ok, we know we have a folder...
        log.info("Using folder: '%s'", self.base_path)


    def load(self):
        if not os.path.exists(self.config_path) or \
           not os.path.isfile(self.config_path):
            log.error("Configuration file '%s' does not exist",
                      self.config_path)
            raise ConfigurationError

        p = Parser(self)
        log.debug("Loading configuration at '%s'", self.config_path)
        p.parse_file(self.config_path)

    def init_log(self):
        """Add a log file in the folder"""
        log_output = logging.FileHandler(self.log_path, mode='w')
        log_output.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
            '%(levelname)-8s | %(asctime)s | %(name)-10s | %(message)s',
            datefmt="%Y-%m-%d %H:%M:%S")
        log_output.setFormatter(formatter)
        logging.getLogger('').addHandler(log_output)

    def verify(self):
        """Check that the parts are consistent"""
        pass

    def has_changed(self):
        """check to see if the config has changed from when it was
        last run"""
        pass
        # Write a signature out somewhere...







