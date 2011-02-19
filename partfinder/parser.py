import logging
log = logging.getLogger("parser")

import os

from pyparsing import (
    Word, Dict, OneOrMore, alphas, nums, Suppress, Optional, Group, stringEnd,
    delimitedList, pythonStyleComment, ParseException, line, lineno, col,
    Keyword, ParserElement, ParseException)

# debugging
# ParserElement.verbose_stacktrace = True

from partition import Partition, PartitionError, all_partitions
from scheme import Scheme, SchemeError
from subset import Subset, SubsetError

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

    # These will get set in the configuration passed in
    required_variables = ['alignment_file']

    def __init__(self, config):
        # Config is filled out with objects that the parser creates
        self.config = config
        self.partitions = []
        self.schemes = []
        self.subsets = []
        self.init_grammar()

    def init_grammar(self):
        """Set up the parsing classes
        Any changes to the grammar of the config file be done here.
        """
        # Some syntax that we need, but don't bother looking at
        SEMIOPT = Optional(Suppress(";"))
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
        general_def = VARIABLE("name") + EQUALS + VALUE("value") + SEMIOPT
        general_def.setParseAction(self.define_variable)
        general = OneOrMore(Group(general_def))
        general.setParseAction(self.check_variables)

        # Partition Parsing
        column = Word(nums)
        partname = Word(alphas + '_-' + nums)
        partdef = column("start") + DASH + column("end") + Optional(BACKSLASH + column("step"))
        partdef.setParseAction(self.define_range)
        partdeflist = Group(OneOrMore(Group(partdef)))
        partition = partname("name") + EQUALS + partdeflist("parts") + SEMIOPT
        partition.setParseAction(self.define_partition)
        partlist = OneOrMore(Group(partition))


        # Scheme Parsing
        schemename = Word(alphas + '_-' + nums)
        partnameref = partname.copy() # Make a copy, cos we set a different action on it
        partnameref.setParseAction(self.check_part_exists)

        subset = Group(OPENB + delimitedList(partnameref("name")) + CLOSEB)
        subset.setParseAction(self.define_subset)

        scheme = Group(OneOrMore(subset))
        schemedef = schemename("name") + EQUALS + scheme("scheme") + SEMIOPT
        schemedef.setParseAction(self.define_schema)

        schemelist = OneOrMore(Group(schemedef))

        # We've defined the grammar for each section. Here we just put it all
        # together
        self.config_parser = (
            general
            + Suppress("[partitions]") + partlist
            + Suppress("[schemes]") + schemelist 
            + stringEnd
        )

    def define_variable(self, text, loc, var_def):
        """Define a variable -- we'll check here if it is allowable"""
        if var_def.name not in self.required_variables:
            raise ParserError(text, loc, "'%s' is not an allowable setting" %
                                 var_def.name)
        else:
            self.config.alignment_file = var_def.value
            log.debug("Setting '%s' to '%s'", var_def.name, var_def.value)

    def check_variables(self, text, loc, var_def):
        # Add the stuff to the configuration that was passed in
        # We should check that all the parameters are defined too...
        for var in self.required_variables:
            if not hasattr(self.config, var):
                raise ParserError(text, loc, "No '%s' defined in the configuration" % var)

    def define_range(self, part):
        """Turn the 2 or 3 tokens into integers, supplying a default if needed"""
        fromc = int(part.start)
        toc = int(part.end)
        if part.step:
            stepc = int(part.step)
        else:
            stepc = 1
        return [fromc, toc, stepc]

    def define_partition(self, text, loc, part_def):
        """We have everything we need here to make a partition"""
        try:
            # Creation adds it to set
            p = Partition(part_def.name, *tuple(part_def.parts))
            self.partitions.append(p)
        except PartitionError:
            raise ParserError(text, loc, "Error in '%s' can be found" % part_def.name)


    def check_part_exists(self, text, loc, partref):
        if partref.name not in all_partitions:
            raise ParserError(text, loc, "Partition %s not defined" %
                                     partref.name)

    def define_subset(self, text, loc, subset_def):
        try:
            # Get the partitions from the names
            parts = [all_partitions[nm] for nm in subset_def[0]]
            # create a subset
            self.subsets.append(Subset(*tuple(parts)))
        except SubsetError:
            raise ParserError(text, loc, "Error creating subset...")
    
    def define_schema(self, text, loc, scheme_def):
        try:
            # Clear out the subsets as we need to reuse it
            subs = tuple(self.subsets)
            self.subsets = []
            self.schemes.append(Scheme(scheme_def.name, *subs))
        except (SchemeError, SubsetError):
            raise ParserError(text, loc, "Error in '%s' can be found" %
                                     scheme_def.name)

    def parse_file(self, fname):
        s = open(fname, 'r').read()
        self.parse_configuration(s)

    def parse_configuration(self, s):
        self.result = self.config_parser.ignore(pythonStyleComment).parseString(s)


        self.config.partitions = self.partitions
        self.config.schemes = self.schemes

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    test_config = r"""
alignment_file = ./test.fas 

# schemes = section # Use the stuff defined below
# schemes = greedy

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
allsame         = (Gene1_pos1, Gene1_pos2, Gene1_pos3, Gene2_pos1, Gene2_pos2,
Gene2_pos3, Gene3_pos1, Gene3_pos2, Gene3_pos3)
by_gene         = (Gene1_pos1, Gene1_pos2, Gene1_pos3) (Gene2_pos1, Gene2_pos2, Gene2_pos3) (Gene3_pos1, Gene3_pos2, Gene3_pos3)
1_2_3           = (Gene1_pos1, Gene2_pos1, Gene3_pos1) (Gene1_pos2, Gene2_pos2, Gene3_pos2) (Gene1_pos3, Gene2_pos3, Gene3_pos3)
1_2_3_by_gene   = (Gene1_pos1) (Gene1_pos2) (Gene1_pos3) (Gene2_pos1) (Gene2_pos2) (Gene2_pos3) (Gene3_pos1) (Gene3_pos2) (Gene3_pos3)
12_3            = (Gene1_pos1, Gene1_pos2, Gene2_pos1, Gene2_pos2, Gene3_pos1, Gene3_pos2) (Gene1_pos3, Gene2_pos3, Gene3_pos3)
12_3_by_gene    = (Gene1_pos1, Gene1_pos2) (Gene1_pos3) (Gene2_pos1, Gene2_pos2) (Gene2_pos3) (Gene3_pos1, Gene3_pos2) (Gene3_pos3)
"""

    class Conf(object):
        pass
    c = Conf()
    p = Parser(c)
    try:
        p.parse_configuration(test_config)
    except ParserError as p:
        log.error(p.format_message())

    # for s in c.schemes:
        # print s.name
        # for ss in s.subsets:
            # print ss.subset_id

        

            # print ss.string_identifier
        # print name
    # else:
        # print p.schemes.subsets
    

