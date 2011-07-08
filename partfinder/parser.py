import logging
log = logging.getLogger("parser")

from pyparsing import (
    Word, OneOrMore, alphas, nums, Suppress, Optional, Group, stringEnd,
    delimitedList, pythonStyleComment, line, lineno, col, Keyword, Or,
    NoMatch)

# debugging
# ParserElement.verbose_stacktrace = True

import partition, scheme, subset, phyml_models, config
from util import PartitionFinderError

# Only used internally
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
    def __init__(self, cfg):
        # For adding variables
        self.cfg = cfg

        # Use these to keep track of stuff that is going on in parser
        self.schemes = []
        self.subsets = []
        self.init_grammar()
        self.ignore_schemes = False
        # provide useful error messages when parsing settings with limited options 
        self.options = {
            'branchlengths': ['linked', 'unlinked'],
            'model_selection': ['AIC', 'AICc', 'BIC'],
            'search': ['all', 'user', 'greedy']
            }

    def init_grammar(self):
        """Set up the parsing classes
        Any changes to the grammar of the config file be done here.
        """
        # Some syntax that we need, but don't bother looking at
        SEMICOLON = (Suppress(";"))
        EQUALS = Suppress("=")
        OPENB = Suppress("(")
        CLOSEB = Suppress(")")
        BACKSLASH = Suppress("\\")
        DASH = Suppress("-")

        # Top Section
        FILENAME = Word(alphas + nums + '-_.')
        alignmentdef = Keyword('alignment') + EQUALS + FILENAME + SEMICOLON
        alignmentdef.setParseAction(self.set_alignment)

        def simple_option(name):
            opt = Keyword(name) + EQUALS + Word(alphas+nums) + SEMICOLON
            opt.setParseAction(self.set_simple_option)
            return opt

        branchdef = simple_option('branchlengths')

        MODELNAME = Word(alphas + nums + '+')
        modellist = delimitedList(MODELNAME)
        modeldef = Keyword("models") + EQUALS + Group(
            (Keyword("all") | Keyword("mrbayes") | Keyword("raxml"))("predefined") | 
            Group(modellist)("userlist")) + SEMICOLON
        modeldef.setParseAction(self.set_models)

    
        modseldef = simple_option("model_selection")
        topsection = alignmentdef + branchdef + modeldef + modseldef

        # Partition Parsing
        column = Word(nums)
        partname = Word(alphas + '_-' + nums)
        partdef = column("start") +\
                Optional(DASH + column("end")) +\
                Optional(BACKSLASH + column("step"))

        partdef.setParseAction(self.define_range)
        partdeflist = Group(OneOrMore(Group(partdef)))
        partition = Optional("charset") + partname("name") + \
                            EQUALS + partdeflist("parts") + SEMICOLON
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
        schemedef = schemename("name") + \
                            EQUALS + scheme("scheme") + SEMICOLON
        schemedef.setParseAction(self.define_schema)

        schemelist = OneOrMore(Group(schemedef))

        schemealgo = simple_option("search")
        schemesection = \
                Suppress("[schemes]") + schemealgo + Optional(schemelist)

        # We've defined the grammar for each section. Here we just put it all together
        self.config_parser = (topsection + partsection + schemesection + stringEnd)
        
    def set_alignment(self, text, loc, tokens):
        value = tokens[1]
        self.cfg.set_alignment_file(value)
        # TODO Make sure it is readable!
        # raise ParserError(text, loc, "No '%s' defined in the configuration" % var)
        #

    def set_simple_option(self, text, loc, tokens):
        try:
            self.cfg.set_option(tokens[0], tokens[1])
        except config.ConfigurationError:
            raise ParserError(text, loc, "Invalid option, see previous error")

        
    def set_models(self, text, loc, tokens):
        all_mods = set(phyml_models.get_all_models())
        mods = tokens[1]
        if mods.userlist:
            self.cfg.models = []
            # It is a list of models
            for m in mods.userlist:
                if m not in all_mods:
                    raise ParserError(
                        text, loc, "'%s' is not a valid model" % m)
                self.cfg.models.append(m)
            log.info("Setting 'models' to a userlist containing: %s", 
                      ", ".join(self.cfg.models))
        else:
            modsgroup = mods.predefined
            if modsgroup.lower() == "all":
                self.cfg.models = list(all_mods)
            elif modsgroup.lower() == "mrbayes":
                mrbayes_mods = set(phyml_models.get_mrbayes_models())
                self.cfg.models = list(mrbayes_mods)
            elif modsgroup.lower() == "raxml":
                self.cfg.models = phyml_models.get_raxml_models()
            else:
                pass
            log.info("Setting 'models' to '%s'", modsgroup)
                    

    def define_range(self, part):
        """Turn the 1, 2 or 3 tokens into integers, supplying a default if needed"""
        fromc = int(part.start)

        if part.end:
            toc = int(part.end)
        else:
            toc = fromc

        if part.step:
            stepc = int(part.step)
        else:
            stepc = 1
        return [fromc, toc, stepc]

    def define_partition(self, text, loc, part_def):
        """We have everything we need here to make a partition"""
        try:
            # Creation adds it to set
            p = partition.Partition(self.cfg, part_def.name, *tuple(part_def.parts))
        except partition.PartitionError:
            raise ParserError(text, loc, "Error in '%s' can be found" % part_def.name)


    def check_part_exists(self, text, loc, partref):
        if partref.name not in self.cfg.partitions:
            raise ParserError(text, loc, "Partition %s not defined" %
                                     partref.name)

    def define_subset(self, text, loc, subset_def):
        try:
            # Get the partitions from the names
            parts = [self.cfg.partitions[nm] for nm in subset_def[0]]

            # Keep a running list of these till we define the schema below
            self.subsets.append(subset.Subset(*tuple(parts)))
        except subset.SubsetError:
            raise ParserError(text, loc, "Error creating subset...")
    
    def define_schema(self, text, loc, scheme_def):
        try:
            # Clear out the subsets as we need to reuse it
            subs = tuple(self.subsets)
            self.subsets = []
            
            if self.ignore_schemes == False:
                scheme.Scheme(self.cfg, scheme_def.name, *subs)

        except (scheme.SchemeError, subset.SubsetError):
            raise ParserError(text, loc, "Error in '%s' can be found" %
                                     scheme_def.name)

    def parse_file(self, fname):
        #this just reads in the config file into 's'
        s = open(fname, 'r').read()
        self.parse_configuration(s)

    def parse_configuration(self, s):
        #parse the config cfg
        try:
            self.result = self.config_parser.ignore(pythonStyleComment).parseString(s)
        except ParserError, p:
            log.error(p.format_message())
            raise PartitionFinderError


if __name__ == '__main__':
    logging.basicConfig()
    # logging.basicConfig(level=logging.DEBUG)
    # import config
    # config.initialise_temp()
    test_config = r"""
# PartitionFinder configuration file
# Anything after a hash (#) is a comment, and is ignored. Feel free to add/remove these lines.
# email rob {dot} lanfear {at} gmail {dot} com for help

######## ALIGNMENT FILE ###########
#the name of your Phylip alignment, in the same file as this file
alignment = test.phy ;

######## BRANCHLENGTHS ###########
# 'linked' or 'unlinked'
branchlengths = linked 

######## MODELS OF EVOLUTION ###########
G# 'all' (all 56 usual models) or a list of models e.g. 'HKY, HKY+G, GTR, GTR+G'
# MRBAYES MODELS
# JC, F81, K80, HKY, SYM, GTR, JC+I, F81+I, K80+I, HKY+I, SYM+I, GTR+I, JC+G, F81+G, K80+G, HKY+G, SYM+G, GTR+G, JC+I+G, F81+I+G, K80+I+G, HKY+I+G, SYM+I+G, GTR+I+G
# RAXML MODELS
# GTR+G, GTR+I+G
models = all

########   PARTITIONS   ###########
# Define partitions as follows 'name = start-stop\gap_size'
# e.g. 'part_1 = 1-15\3' is the same as 'part_1 = 1 4 7 10 13' 
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

########     SCHEMES      #########
# 'all' (compares all possible schemes) or 'user' (searches schemes listed below) 
[schemes]
search = user

#user schemes listed below - only used if 'search = user'
allsame         = (Gene1_pos1, Gene1_pos2, Gene1_pos3, Gene2_pos1, Gene2_pos2, Gene2_pos3, Gene3_pos1, Gene3_pos2, Gene3_pos3)
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
    except InternalParserError, p:
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
    

