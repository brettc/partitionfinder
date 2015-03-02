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
import pandas as pd
import os
import collections

log = logtools.get_logger()
from util import PartitionFinderError

_available_lists = ["dna",
                    "dna_total", 
                    "protein",
                    "protein_total", 
                    "beast", 
                    "mrbayes", 
                    "protein_gamma", 
                    "protein_gammaI"]

def load_models(the_config):
    HERE = os.path.abspath(os.path.dirname(__file__))
    all_models = pd.read_csv(os.path.join(HERE, 'models.csv'))

    # determine available models based on datatype and phylogeny program
    the_config.available_models = get_available_models(all_models, the_config)

    # check user models will run
    parse_user_models(the_config)

    log.info("This analysis will use the following %d models of molecular evolution"
             % len(the_config.models))
    log.info("%s" % ', '.join(the_config.models))


def get_available_models(all_models, the_config):
    # from the list of all models, which ones could we actually run
    if the_config.phylogeny_program == 'phyml':
        available_models = all_models[pd.notnull(all_models.phyml_commandline)]
    elif the_config.phylogeny_program == 'raxml':
        available_models = all_models[pd.notnull(all_models.raxml_commandline)]

    if the_config.datatype == 'DNA':
        available_models = available_models.query("datatype=='DNA'")
    elif the_config.datatype == 'protein':
        available_models = available_models.query("datatype=='protein'")
    elif the_config.datatype == 'morphology':
        available_models = available_models.query("datatype=='morphology'")

    if len(available_models) == 0:
        log.error("""Phylogeny program '%s' does not implement any models that deal 
                  with %s data. Please check and try again. For morphological data,
                  use RAxML (--raxml at the commandline)""" % 
                  (the_config.phylogeny_program, the_config.datatype))
        raise PartitionFinderError

    return available_models

def parse_user_models(the_config):

    # this will tell us if they entered any lists or models we can't use
    check_all_models_and_lists(the_config)

    check_for_duplicates(the_config.models)

    mod_list = check_model_lists(the_config.models)

    if mod_list:
        expand_model_list(the_config)

    # final check on models
    check_all_models(the_config)

def check_all_models(the_config):
    # everything has to be a model in the_config.available_models
    models = the_config.models
    allowed = set(the_config.available_models.name)

    problems = set(models).difference(allowed)

    if problems:
        log.error("""'%s' is/are not a valid model(s) for phylogeny program %s
                  and data type %s, please check and try again""" 
                  %(', '.join(problems), the_config.phylogeny_program, the_config.datatype))
        log.info("""If you are unsure which models are available, or why a model you think 
                 should work does not, please check the manual and the models.csv file 
                 (located in the /partfinder folder) for more information.""")
        raise PartitionFinderError

def check_all_models_and_lists(the_config):
    # everything has to be either a model in the_config.available_models
    # OR a valid option from the _available_models
    models = the_config.models
    allowed = set(_available_lists).union(set(the_config.available_models.name))

    problems = set(models).difference(allowed)

    if problems:
        log.error("""'%s' is/are not a valid model(s) or lists of models 
                  for phylogeny program %s and data type %s, 
                  please check and try again.""" 
                  %(', '.join(problems), the_config.phylogeny_program, the_config.datatype))

        log.info("""If you are unsure which models are available, or why a model you think 
                 should work does not, please check the manual and the models.csv file 
                 (located in the /partfinder folder) for more information.""")
        raise PartitionFinderError
        

def check_model_lists(models):

    mod_lists = set(models).intersection(set(_available_lists))

    if mod_lists and len(models)>1:
        log.error("""If you use a model list (you used '%s') you can only
                  specify a single list, and no other lists or models. Please 
                  check and try again.""" % ', '.join(mod_lists))
        raise PartitionFinderError

    return mod_lists


def expand_model_list(the_config):

    # by this point, we know that mod_list is a list of length 1
    mod_list = the_config.models[0]

    the_config.models = list(the_config.available_models.query("%s==1" % mod_list).name)

    if len(the_config.models)<1:
        log.error("""The model list '%s' is not a compatible with 
                  for phylogeny program %s and data type %s. 
                  There are no models in that list which work with 
                  that combination of program and data type. Please check and try again.""" 
                  %(mod_list, the_config.phylogeny_program, the_config.datatype))
        raise PartitionFinderError


def get_mrbayes_models():
    """
    Return a list of all models implemented in MrBayes. Thanks to Ainsley Seago
    for this.
    """
    mrbayes_base_models = ["JC", "F81", "K80", "HKY", "SYM", "GTR"]
    model_list = []
    for model in mrbayes_base_models:
        model_list.append(model)
        model_list.append("%s+I" % model)
        model_list.append("%s+G" % model)
        model_list.append("%s+I+G" % model)
    return model_list

def get_beast_models():
    """
    Return a list of all models implemented in BEAST v1.7.2.
    """
    beast_base_models = ["K80", "TrNef", "SYM", "HKY", "TrN", "GTR"]
    model_list = []
    for model in beast_base_models:
        model_list.append(model)
        model_list.append("%s+I" % model)
        model_list.append("%s+G" % model)
        model_list.append("%s+I+G" % model)
    return model_list

def get_raxml_models():
    """
    Return a list of all models implemented in RaxML. Thanks to Ainsley Seago
    for this.
    """
    model_list = ["GTR+G", "GTR+I+G"]
    return model_list


def check_for_duplicates(models):
    # model lists shouldn't contain duplicated models
    duplicates = [x for x, y in collections.Counter(models).items() if y > 1]
    if len(duplicates)>0:
        log.error("""There was a problem loading your list of models,
                  the following models seem to be duplicated: %s"""
                  % duplicates)
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
# cut down models to only those in the DNA/Protein and PhyML/RAxML models