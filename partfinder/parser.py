# Copyright (C) 2012 Robert Lanfear and Brett Calcott
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details. You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# PartitionFinder also includes the PhyML program, the RAxML program, and the
# PyParsing library, all of which are protected by their own licenses and
# conditions, using PartitionFinder implies that you agree with those licences
# and conditions as well.

import logtools
log = logtools.get_logger(__file__)

from pyparsing import (
    Word, OneOrMore, alphas, nums, Suppress, Optional, Group, stringEnd,
    delimitedList, pythonStyleComment, line, lineno, col, Keyword,
    CaselessKeyword, ParseException )

# Use this for debugging
# ParserElement.verbose_stacktrace = True

import scheme
import subset
import subset_ops
import phyml_models
import raxml_models
import config
from util import PartitionFinderError


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
    def __init__(self, cfg):
        """Construct parser from the configuration
        """
        self.cfg = cfg

        # Use these to keep track of stuff that is going on in parser
        self.schemes = []
        self.current_subsets = []
        self.init_grammar()
        self.ignore_schemes = False

    def init_grammar(self):
        """Set up the parsing classes
        Any changes to the grammar of the config file be done here.
        """

        # Some syntax that we need, but don't care about
        SEMICOLON = (Suppress(";"))
        EQUALS = Suppress("=")

        # Top Section
        FILE_NAME = Word(alphas + nums + '-_.')
        alignment_def = Keyword('alignment') + EQUALS\
                        + FILE_NAME + SEMICOLON
        alignment_def.setParseAction(self.set_alignment)

        tree_def = Keyword('user_tree_topology') + EQUALS\
                   + FILE_NAME + SEMICOLON
        tree_def.setParseAction(self.set_user_tree)

        def simple_option(name):
            opt = Keyword(name) + EQUALS +\
                  Word(alphas + nums + '-_') + SEMICOLON
            opt.setParseAction(self.set_simple_option)
            return opt

        branch_def = simple_option('branchlengths')

        MODEL_NAME = Word(alphas + nums + '+')
        model_list = delimitedList(MODEL_NAME)
        model_def = Keyword("models") + EQUALS + Group((
            CaselessKeyword("all") |
            CaselessKeyword("mrbayes") |
            CaselessKeyword("raxml") |
            CaselessKeyword("beast") |
            CaselessKeyword("all_protein") |
            CaselessKeyword("all_protein_gamma") |
            CaselessKeyword("all_protein_gammaI"))("predefined") |
            Group(model_list)("userlist")) + SEMICOLON
        model_def.setParseAction(self.set_models)

        model_selection_def = simple_option("model_selection")
        top_section = alignment_def + Optional(tree_def) + branch_def + \
            model_def + model_selection_def

        # Data Block Parsing
        column = Word(nums)
        block_name = Word(alphas + '_-' + nums)
        block_def = column("start") +\
            Optional(Suppress("-") + column("end")) +\
            Optional(Suppress("\\") + column("step"))
        block_def.setParseAction(self.define_range)

        block_list_def = Group(OneOrMore(Group(block_def)))

        user_subset_def = Optional("charset") + block_name("name") + \
            EQUALS + block_list_def("parts") + SEMICOLON
        user_subset_def.setParseAction(self.define_user_subset)

        block_def_list = OneOrMore(Group(user_subset_def))
        block_section = Suppress("[data_blocks]") + block_def_list

        # Scheme Parsing
        scheme_name = Word(alphas + '_-' + nums)

        # Make a copy, cos we set a different action on it
        user_subset_ref = block_name.copy()
        user_subset_ref.setParseAction(self.check_block_exists)

        subset = Group(Suppress("(") +
                       delimitedList(user_subset_ref("name")) + Suppress(")"))
        subset.setParseAction(self.define_subset_grouping)

        scheme = Group(OneOrMore(subset))
        scheme_def = scheme_name("name") + \
            EQUALS + scheme("scheme") + SEMICOLON
        scheme_def.setParseAction(self.define_scheme)

        scheme_list = OneOrMore(Group(scheme_def))

        scheme_algo = simple_option("search")
        scheme_section = \
            Suppress("[schemes]") + scheme_algo + Optional(scheme_list)

        # We've defined the grammar for each section.
        # Here we just put it all together
        self.config_parser = (
            top_section + block_section + scheme_section + stringEnd)

    def set_alignment(self, text, loc, tokens):
        value = tokens[1]
        self.cfg.set_alignment_file(value)

    def set_user_tree(self, tokens):
        self.cfg.user_tree = tokens[1]

    def set_simple_option(self, text, loc, tokens):
        try:
            self.cfg.set_option(tokens[0], tokens[1])
        except config.ConfigurationError:
            raise ParserError(text, loc, "Invalid option in .cfg file")

    def define_range(self, text, loc, part):
        """Turn the 1, 2 or 3 tokens into integers: from:to/step
        """
        from_column = int(part.start)

        if part.end:
            to_column = int(part.end)
        else:
            to_column = from_column

        if part.step:
            step = int(part.step)
        else:
            step = 1

        if from_column > to_column:
            raise ParserError(
                text, loc, "Data block has must have beginning less than end "
                "(%d is great than %d)" % (from_column, to_column))
        return [from_column, to_column, step]

    def define_user_subset(self, text, loc, part_def):
        """Use a list of tuples with start,stop,step to produces columns"""

        # We need to convert to column definitions. Note that these are
        # zero based, which is not how they are specified in the config. So we
        # must do some fiddling to make sure they are right. In addition, we
        # use range(...) which excludes the final column, whereas the
        # definitions assume inclusive. Subtracting 1 from start deals with
        # both issues...
        columns = []
        description = []
        for start, stop, step in part_def.parts:
            columns.extend(range(start-1, stop, step))

            # Keep a description of this around
            description.append((start, stop, step))

        # Normalise it all
        column_set = set(columns)

        # If there was any overlap then these will differ...
        if len(columns) != len(column_set):
            raise ParserError(
                text, loc, "Block '%s' has internal overlap" % part_def.name)

        user_subset = subset.Subset(self.cfg, column_set)

        # TODO: Think about how we want to add descriptions.
        # Maybe this should be part of the __init__
        user_subset.add_description([part_def.name], [tuple(description)])
        self.cfg.user_subsets.append(user_subset)
        self.cfg.user_subsets_by_name[part_def.name] = user_subset

    def check_block_exists(self, text, loc, partref):
        if partref.name not in self.cfg.user_subsets_by_name:
            raise ParserError(text, loc, "Block %s not defined" %
                              partref.name)

    def define_subset_grouping(self, text, loc, subset_def):
        """These define initial groupings that users think are useful
        """
        try:
            # Get the partitions from the names
            subsets = [self.cfg.user_subsets_by_name[nm] for nm in subset_def[0]]

            # Keep a running list of these till we define the schema below
            self.current_subsets.append(subset_ops.merge_subsets(subsets))
        except subset.SubsetError:
            raise ParserError(text, loc, "Error creating subset...")

    def define_scheme(self, text, loc, scheme_def):
        try:
            subs = tuple(self.current_subsets)

            if not self.ignore_schemes:
                sch = scheme.Scheme(self.cfg, scheme_def.name, subs)
                self.cfg.user_schemes.add_scheme(sch)

        except (scheme.SchemeError, subset.SubsetError):
            raise ParserError(text, loc, "Error in '%s' can be found" %
                              scheme_def.name)
        finally:
            # Clear out the current_subsets as we need to reuse it
            self.current_subsets = []


    def parse_file(self, file_name):
        #this just reads in the config file into 's'
        s = open(file_name, 'rU').read()
        self.parse_configuration(s)

    def parse_configuration(self, s):
        """Parse a string as a configuration settings file"""
        try:
            self.result = self.config_parser.ignore(pythonStyleComment).\
                parseString(s)

        except ParserError, p:
            log.error(p.format_message())
            raise PartitionFinderError

        except ParseException, p:
            log.error("""There was a problem loading your .cfg file, please
                      check and try again""")
            log.error(str(p))

            # Let's see if there was something missing from the input file
            expectations = ["models", "search", "[schemes]", "[data_blocks]",
                            "model_selection", "branchlengths", "alignment"]
            missing = None
            for e in expectations:
                if p.msg.count(e):
                    missing = e

            if missing:
                log.info(
                    """It looks like the '%s' option might be missing or in the
                    wrong place. Or perhaps something is wrong in the lines just
                    before the '%s' option is missing. Please double check the
                    configuration file and try again""" % (missing, missing))
            else:
                log.info(
                    """The line causing the problem is this: '%s'. Please
                    check that line, and make sure it appears in the right
                    place in the config file. If it looks OK, try
                    double-checking the semi-colons on other lines""" % p.line)
            raise PartitionFinderError

    def set_models(self, text, loc, tokens):
        # TODO: Fix this ugly mess.
        if self.cfg.phylogeny_program == "phyml":
            self.phylo_models = phyml_models
        elif self.cfg.phylogeny_program == "raxml":
            self.phylo_models = raxml_models

        all_dna_mods = set(self.phylo_models.get_all_dna_models())
        all_protein_mods = set(self.phylo_models.get_all_protein_models())
        total_mods = all_dna_mods | all_protein_mods

        mods = tokens[1]
        DNA_mods = 0
        protein_mods = 0
        if mods.userlist:
            modlist = mods.userlist
            log.info("Setting 'models' to a user-specified list")
        else:
            modsgroup = mods.predefined
            if modsgroup.lower() == "all":
                modlist = list(all_dna_mods)
                DNA_mods += 1
            elif modsgroup.lower() == "mrbayes":
                modlist = set(phyml_models.get_mrbayes_models())
                DNA_mods += 1
            elif modsgroup.lower() == "beast":
                modlist = set(phyml_models.get_beast_models())
                DNA_mods += 1
            elif modsgroup.lower() == "raxml":
                modlist = set(phyml_models.get_raxml_models())
                DNA_mods += 1
            elif modsgroup.lower() == "all_protein":
                modlist = set(self.phylo_models.get_all_protein_models())
                protein_mods += 1
            elif modsgroup.lower() == "all_protein_gamma":
                if self.cfg.phylogeny_program == "raxml":
                    modlist = set(raxml_models.get_protein_models_gamma())
                    protein_mods += 1
                else:
                    log.error(
                        """The models option 'all_protein_gamma' is only
                        available with raxml  (the --raxml commandline option).
                        Please check and try again""")
                    raise ParserError
            elif modsgroup.lower() == "all_protein_gammaI":
                if self.cfg.phylogeny_program == "raxml":
                    modlist = set(raxml_models.get_protein_models_gammaI())
                    protein_mods += 1
                else:
                    log.error(
                        """The models option 'all_protein_gammaI' is only
                        available with raxml (the --raxml commandline option).
                        Please check and try again""")
                    raise ParserError
            else:
                pass

            # never include the LG4X model in predefined model lists
            # because it can't (yet) be used for partitionined analyses
            modlist = filter(lambda x: x.count("LG4X")==0, modlist)


            log.info("Setting 'models' to '%s'" % modsgroup)

        self.cfg.models = set()
        for m in modlist:
            if m not in total_mods:
                raise ParserError(
                    text, loc, "'%s' is not a valid model for phylogeny "
                               "program %s. Please check the lists of valid models in the"
                               " manual and try again" % (m, self.cfg.phylogeny_program))

            if m in all_dna_mods:
                DNA_mods += 1
            if m in all_protein_mods:
                protein_mods += 1

            self.cfg.models.add(m)

        log.info("The models included in this analysis are: %s",
                 ", ".join(self.cfg.models))
        self.cfg.model_count = len(self.cfg.models)

        # Check data type against the model list that we've got a sensible
        # model list
        if DNA_mods > 0 and protein_mods == 0 and self.cfg.datatype == "DNA":
            log.info("Setting datatype to 'DNA'")
        elif DNA_mods == 0 and protein_mods > 0 and self.cfg.datatype == "protein":
            log.info("Setting datatype to 'protein'")
        elif DNA_mods == 0 and protein_mods > 0 and self.cfg.datatype == "DNA":
            raise ParserError(
                text, loc, "The models list contains only models of amino acid change."
                           " PartitionFinder.py only works with nucleotide models (like the GTR model)."
                           " If you're analysing an amino acid dataset, please use PartitionFinderProtein,"
                           " which you can download here: www.robertlanfear.com/partitionfinder."
                           " The models line in the .cfg file is")
        elif DNA_mods > 0 and protein_mods == 0 and self.cfg.datatype == "protein":
            raise ParserError(
                text, loc, "The models list contains only models of nucleotide change."
                           " PartitionFinderProtein.py only works with amino acid models (like the WAG model)."
                           " If you're analysing a nucleotide dataset, please use PartitionFinder.py,"
                           " which you can download here: www.robertlanfear.com/partitionfinder"
                           " The models line in the .cfg file is")
        else:  # we've got a mixture of models.
            raise ParserError(
                text, loc, "The models list contains a mixture of protein and nucleotide models."
                           " If you're analysing a nucleotide dataset, please use PartitionFinder."
                           " If you're analysing an amino acid dataset, please use PartitionFinderProtein."
                           " You can download both of these programs from here: www.robertlanfear.com/partitionfinder"
                           " The models line in the .cfg file is")
