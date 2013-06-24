#Copyright (C) 2012 Robert Lanfear and Brett Calcott
#
#This program is free software: you can redistribute it and/or modify it
#under the terms of the GNU General Public License as published by the
#Free Software Foundation, either version 3 of the License, or (at your
#option) any later version.
#
#This program is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#General Public License for more details. You should have received a copy
#of the GNU General Public License along with this program.  If not, see
#<http://www.gnu.org/licenses/>. PartitionFinder also includes the PhyML
#program, the RAxML program, and the PyParsing library,
#all of which are protected by their own licenses and conditions, using
#PartitionFinder implies that you agree with those licences and conditions as well.

import logging
log = logging.getLogger("parser")

from pyparsing import (
    Word, OneOrMore, alphas, nums, Suppress, Optional, Group, stringEnd,
    delimitedList, pythonStyleComment, line, lineno, col, Keyword, Or,
    NoMatch, CaselessKeyword, ParseException, SkipTo)

# debugging
# ParserElement.verbose_stacktrace = True

import partition
import scheme
import subset
import phyml_models
import raxml_models
import config
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

        treedef = Keyword('user_tree_topology') + EQUALS + FILENAME + SEMICOLON
        treedef.setParseAction(self.set_user_tree)

        def simple_option(name):
            opt = Keyword(name) + EQUALS + Word(alphas + nums + '-_') + SEMICOLON
            opt.setParseAction(self.set_simple_option)
            return opt

        branchdef = simple_option('branchlengths')

        MODELNAME = Word(alphas + nums + '+')
        modellist = delimitedList(MODELNAME)
        modeldef = Keyword("models") + EQUALS + Group(
            (
                CaselessKeyword("all") | CaselessKeyword("mrbayes") | CaselessKeyword("raxml") |
                CaselessKeyword("beast") | CaselessKeyword("all_protein") |
                CaselessKeyword(
                    "all_protein_gamma") | CaselessKeyword("all_protein_gammaI")
            )("predefined") |
            Group(modellist)("userlist")) + SEMICOLON
        modeldef.setParseAction(self.set_models)

        modseldef = simple_option("model_selection")
        topsection = alignmentdef + Optional(treedef) + branchdef + \
            modeldef + modseldef

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
        partsection = Suppress("[data_blocks]") + partlist

        # Scheme Parsing
        schemename = Word(alphas + '_-' + nums)
        partnameref = partname.copy(
        )  # Make a copy, cos we set a different action on it
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
        self.config_parser = (
            topsection + partsection + schemesection + stringEnd)

    def set_alignment(self, text, loc, tokens):
        value = tokens[1]
        self.cfg.set_alignment_file(value)
        # TODO Make sure it is readable!
        # raise ParserError(text, loc, "No '%s' defined in the configuration" % var)
        #

    def set_user_tree(self, text, loc, tokens):
        self.cfg.user_tree = tokens[1]
        pass

    def set_simple_option(self, text, loc, tokens):
        try:
            self.cfg.set_option(tokens[0], tokens[1])
        except config.ConfigurationError:
            raise ParserError(text, loc, "Invalid option in .cfg file")

    def set_models(self, text, loc, tokens):
        if self.cfg.phylogeny_program == "phyml":
            self.phylo_models = phyml_models
        elif self.cfg.phylogeny_program == "raxml":
            self.phylo_models = raxml_models

        all_dna_mods = set(self.phylo_models.get_all_dna_models())
        all_protein_mods = set(self.phylo_models.get_all_protein_models())
        all_morphology_mods = set(self.phylo_models.get_all_morphology_models())
        total_mods = all_dna_mods.union(all_protein_mods)
        total_mods = total_mods.union(all_morphology_mods)

        mods = tokens[1]
        DNA_mods = 0
        prot_mods = 0
        if mods.userlist:
            modlist = mods.userlist
            log.info("Setting 'models' to a user-specified list")
        else:
            modsgroup = mods.predefined
            if modsgroup.lower() == "all":
                modlist = list(all_dna_mods)
                DNA_mods = DNA_mods + 1
            elif modsgroup.lower() == "mrbayes":
                modlist = set(phyml_models.get_mrbayes_models())
                DNA_mods = DNA_mods + 1
            elif modsgroup.lower() == "beast":
                modlist = set(phyml_models.get_beast_models())
                DNA_mods = DNA_mods + 1
            elif modsgroup.lower() == "raxml":
                modlist = set(phyml_models.get_raxml_models())
                DNA_mods = DNA_mods + 1
            elif modsgroup.lower() == "all_protein":
                modlist = set(self.phylo_models.get_all_protein_models())
                prot_mods = prot_mods + 1
            elif modsgroup.lower() == "all_protein_gamma":
                if self.cfg.phylogeny_program == "raxml":
                    modlist = set(raxml_models.get_protein_models_gamma())
                    prot_mods = prot_mods + 1
                else:
                    log.error("The models option 'all_protein_gamma' is only available with raxml"
                              ", (the --raxml commandline option). Please check and try again")
                    raise ParserError
            elif modsgroup.lower() == "all_protein_gammaI":
                if self.cfg.phylogeny_program == "raxml":
                    modlist = set(raxml_models.get_protein_models_gammaI())
                    prot_mods = prot_mods + 1
                else:
                    log.error("The models option 'all_protein_gammaI' is only available with raxml"
                              ", (the --raxml commandline option). Please check and try again")
                    raise ParserError
            else:
                pass
            log.info("Setting 'models' to '%s'", modsgroup)

        self.cfg.models = set()
        for m in modlist:
            if m not in total_mods:
                raise ParserError(
                    text, loc, "'%s' is not a valid model for phylogeny "
                               "program %s. Please check the lists of valid models in the"
                               " manual and try again" % (m, self.cfg.phylogeny_program))

            if m in all_dna_mods:
                DNA_mods = DNA_mods + 1
            if m in all_protein_mods:
                prot_mods = prot_mods + 1

            self.cfg.models.add(m)

        log.info("The models included in this analysis are: %s",
                 ", ".join(self.cfg.models))

        #check datatype against the model list that we've got a sensible model list
        if DNA_mods > 0 and prot_mods == 0 and self.cfg.datatype == "DNA":
            log.info("Setting datatype to 'DNA'")
        elif DNA_mods == 0 and prot_mods > 0 and self.cfg.datatype == "protein":
            log.info("Setting datatype to 'protein'")
        elif DNA_mods == 0 and prot_mods > 0 and self.cfg.datatype == "DNA":
            raise ParserError(
                text, loc, "The models list contains only models of amino acid change."
                " PartitionFinder.py only works with nucleotide models (like the GTR model)."
                " If you're analysing an amino acid dataset, please use PartitionFinderProtein,"
                " which you can download here: www.robertlanfear.com/partitionfinder."
                " The models line in the .cfg file is")
        elif DNA_mods > 0 and prot_mods == 0 and self.cfg.datatype == "protein":
            raise ParserError(
                text, loc, "The models list contains only models of nucelotide change."
                " PartitionFinderProtein.py only works with amino acid models (like the WAG model)."
                " If you're analysing a nucelotide dataset, please use PartitionFinder.py,"
                " which you can download here: www.robertlanfear.com/partitionfinder"
                " The models line in the .cfg file is")
        #else:  # we've got a mixture of models.
            #raise ParserError(
            #    text, loc, "The models list contains a mixture of protein and nucelotide models."
            #    " If you're analysing a nucelotide dataset, please use PartitionFinder."
            #    " If you're analysing an amino acid dataset, please use PartitionFinderProtein."
            #    " You can download both of these programs from here: www.robertlanfear.com/partitionfinder"
            #    " The models line in the .cfg file is")

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
            p = partition.Partition(
                self.cfg, part_def.name, *tuple(part_def.parts))
        except partition.PartitionError:
            raise ParserError(
                text, loc, "Error in '%s' can be found" % part_def.name)

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
                sch = scheme.Scheme(self.cfg, scheme_def.name, subs)
                self.cfg.user_schemes.add_scheme(sch)

        except (scheme.SchemeError, subset.SubsetError):
            raise ParserError(text, loc, "Error in '%s' can be found" %
                              scheme_def.name)

    def parse_file(self, fname):
        #this just reads in the config file into 's'
        s = open(fname, 'rU').read()
        self.parse_configuration(s)

    def parse_configuration(self, s):
        #parse the config cfg
        try:
            self.result = self.config_parser.ignore(
                pythonStyleComment).parseString(s)
        except ParserError, p:
            log.error(p.format_message())
            raise PartitionFinderError
        except ParseException, p:
            log.error("There was a problem loading your .cfg file, please check and try again")
            log.error(p)

            #let's see if there was something missing fro the input file
            expectations = ["models", "search", "[schemes]", "[data_blocks]",
                            "model_selection", "branchlengths", "alignment"]
            missing = None
            for e in expectations:
                if p.msg.count(e):
                    missing = e

            if missing:
                log.info("It looks like the '%s' option might be missing or in the wrong place" % (missing))
                log.info("Or perhaps something is wrong in the lines just before the '%s' option" % (missing))
                log.info("Please double check the .cfg file and try again")
            else:
                log.info(
                    "The line causing the problem is this: '%s'" % (p.line))
                log.info("Please check that line, and make sure it appears in the right place in the .cfg file.")
                log.info("If it looks OK, try double-checking the semi-colons on other lines in the .cfg file")
            raise PartitionFinderError
