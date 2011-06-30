import logging
log = logging.getLogger("parser")

from pyparsing import (
    Word, OneOrMore, alphas, nums, Suppress, Optional, Group, stringEnd,
    delimitedList, pythonStyleComment, line, lineno, col, Keyword)

# debugging
# ParserElement.verbose_stacktrace = True

import partition, scheme, subset, phyml_models

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
    def __init__(self, settings):
        # For adding variables
        self.settings = settings

        # Use these to keep track of stuff that is going on in parser
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

        # Top Section
        FILENAME = Word(alphas + nums + '-_.')
        alignmentdef = Keyword('alignment') + EQUALS + FILENAME + SEMIOPT
        alignmentdef.setParseAction(self.set_alignment)

        branchdef = Keyword("branch_lengths") + EQUALS \
                + (Keyword("linked") | Keyword("unlinked")) + SEMIOPT
        branchdef.setParseAction(self.set_branchlengths)

        MODELNAME = Word(alphas + nums + '+')
        modellist = delimitedList(MODELNAME)
        modeldef = Keyword("models") + EQUALS + Group(
            (Keyword("all") | Keyword("mrbayes"))("predefined") | 
            Group(modellist)("userlist")) + SEMIOPT
        modeldef.setParseAction(self.set_models)
        topsection = alignmentdef + branchdef + modeldef

        # Partition Parsing
        column = Word(nums)
        partname = Word(alphas + '_-' + nums)
        partdef = column("start") + DASH + column("end") + Optional(BACKSLASH + column("step"))
        partdef.setParseAction(self.define_range)
        partdeflist = Group(OneOrMore(Group(partdef)))
        partition = Optional("charset") + partname("name") + EQUALS + partdeflist("parts") + SEMIOPT
        partition.setParseAction(self.define_partition)
        partlist = OneOrMore(Group(partition))
        partsection = Suppress("[partitions]") + partlist

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

        schemealgo = Keyword("search") + EQUALS + (
            Keyword("all") | Keyword("greedy") | Keyword("user"))
        schemealgo.setParseAction(self.set_scheme_algorithm)
        schemesection = \
                Suppress("[schemes]") + schemealgo + Optional(schemelist)

        # We've defined the grammar for each section. Here we just put it all
        # together
        self.config_parser = (topsection + partsection + schemesection + stringEnd)

    def set_alignment(self, text, loc, tokens):
        value = tokens[1]
        log.debug("Setting 'alignment' to %s", value)
        self.settings.alignment = value
        # TODO Make sure it is readable!
        # raise ParserError(text, loc, "No '%s' defined in the configuration" % var)
        
    def set_branchlengths(self, text, loc, tokens):
        value = tokens[1]
        log.debug("Setting 'branchlengths' to %s", value)
        self.settings.branchlengths = value

    def set_scheme_algorithm(self, text, loc, tokens):
        value = tokens[1]
        log.debug("Setting 'search_algorithm' to %s", value)
        self.settings.search_algorithm = value

    def set_models(self, text, loc, tokens):
        all_mods = set(phyml_models.get_all_models())
        mods = tokens[1]
        if mods.userlist:
            self.settings.models = []
            # It is a list of models
            for m in mods.userlist:
                if m not in all_mods:
                    raise ParserError(
                        text, loc, "'%s' is not a valid model" % m)
                self.settings.models.append(m)
            log.debug("Setting 'models' to given userlist containing %s", 
                      ", ".join(self.settings.models))
        else:
            modsgroup = mods.predefined
            if modsgroup == "all":
                self.settings.models = list(all_mods)
            else:
                pass
            log.debug("Setting 'models' to given predefined list called '%s'",
                      modsgroup)
                    

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
            p = partition.Partition(part_def.name, *tuple(part_def.parts))
            self.partitions.append(p)
        except partition.PartitionError:
            raise ParserError(text, loc, "Error in '%s' can be found" % part_def.name)


    def check_part_exists(self, text, loc, partref):
        if partref.name not in partition.all_partitions:
            raise ParserError(text, loc, "Partition %s not defined" %
                                     partref.name)

    def define_subset(self, text, loc, subset_def):
        try:
            # Get the partitions from the names
            parts = [partition.all_partitions[nm] for nm in subset_def[0]]
            # create a subset
            self.subsets.append(subset.Subset(*tuple(parts)))
        except subset.SubsetError:
            raise ParserError(text, loc, "Error creating subset...")
    
    def define_schema(self, text, loc, scheme_def):
        try:
            # Clear out the subsets as we need to reuse it
            subs = tuple(self.subsets)
            self.subsets = []
            self.schemes.append(scheme.Scheme(scheme_def.name, *subs))
        except (scheme.SchemeError, subset.SubsetError):
            raise ParserError(text, loc, "Error in '%s' can be found" %
                                     scheme_def.name)

    def parse_file(self, fname):
        s = open(fname, 'r').read()
        self.parse_configuration(s)

    def parse_configuration(self, s):
        self.result = self.config_parser.ignore(pythonStyleComment).parseString(s)

if __name__ == '__main__':
    logging.basicConfig()
    # logging.basicConfig(level=logging.DEBUG)
    # import config
    # config.initialise_temp()
    test_config = r"""
alignment = test.fas 



# models = JC 
#
models = all

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
search = all

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
    except ParserError, p:
        log.error(p.format_message())
    
    print c.__dict__
    

    # for s in c.schemes:
        # print s.name
        # for ss in s.subsets:
            # print ss.subset_id

        

            # print ss.string_identifier
        # print name
    # else:
        # print p.schemes.subsets
    

