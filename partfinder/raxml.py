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
log = logtools.get_logger()

import os
import sys
import fnmatch
import util
from database import DataLayout, DataRecord
from reporter import write_raxml_partitions
from config import the_config

from pyparsing import (
    Word, Literal, nums, Suppress, ParseException,
    SkipTo, OneOrMore, Regex, restOfLine, Optional
)

import raxml_models as models

_protein_letters = "ARNDCQEGHILKMFPSTWYV"
_dna_letters = "ATCG"
_morph_chars = "0123456789"

# This is set as the binary name because the previously compiled raxml had a
# bug when calculating site likelihoods, this needs to be changed back to
# "raxml" once a newer version without the bug is compiled.
_binary_name = 'raxml'
_binary_name_pthreads = 'raxml_pthreads'
if sys.platform == 'win32':
    _binary_name            += ".exe"
    _binary_name_pthreads   += ".exe"
if sys.platform == "linux" or sys.platform == "linux2":
    _binary_name            += ".linux"
    _binary_name_pthreads   += ".linux"

_raxml_binary = None
_raxml_pthreads_binary = None


def make_data_layout(cfg):
    if cfg.datatype == "protein":
        letters = _protein_letters
    elif cfg.datatype == "DNA":
        letters = _dna_letters
    elif cfg.datatype == "morphology":
        letters = _morph_chars
    return DataLayout(letters)


def run_raxml(command):
    global _raxml_binary
    if _raxml_binary is None:
        _raxml_binary = util.find_program(_binary_name)
    util.run_program(_raxml_binary, command)

def run_raxml_pthreads(command, cpus):
    global _raxml_pthreads_binary
    if _raxml_pthreads_binary is None:
        _raxml_pthreads_binary = util.find_program(_binary_name_pthreads)
    command = " ".join([command, "-T", str(cpus), " "])
    util.run_program(_raxml_pthreads_binary, command)


def write_partition_file(scheme, alignment_path):
    ''' Write a raxml partitions file to make an ML tree '''
    aln_dir, fname = os.path.split(alignment_path)
    partition_filepath = os.path.join(aln_dir, 'partitions.txt')
    partition_filehandle = open(partition_filepath, 'w')
    sorted_subsets = [sub for sub in scheme]
    sorted_subsets.sort(key=lambda sub: min(sub.columns), reverse=False)
    write_raxml_partitions(scheme, partition_filehandle, sorted_subsets, use_lg = True)
    return(partition_filepath)


def make_ml_topology(alignment_path, datatype, cmdline_extras, scheme, cpus):
    '''Make a ML tree to from a given partitioning scheme'''
    log.info("Estimating Maximum Likelihood tree with RAxML fast experimental tree search for %s", alignment_path)

    if(the_config.datatype != "morphology"):
        partition_file = write_partition_file(scheme, alignment_path)


    # First get the ML topology like this (-p is a hard-coded random number seed):
    # we do this to an accuracy of 10 log likelihood units with -e 10
    # we use the rapid ML option in RAxML -f E
    if datatype == "DNA":
        log.info("Using a separate GTR+G model for each data block")
        command = " -f E -s '%s' -m GTRGAMMA -O -n fastTREE -# 1 -p 123456789 -q '%s' -e 10 " % (
            alignment_path, partition_file)
    elif datatype == "protein":
        log.info("Using a separate LG+G model for each data block")
        command = " -f E -s '%s' -m PROTGAMMALG -O -n fastTREE -# 1 -p 123456789 -q '%s' -e 10 " % (
            alignment_path, partition_file)
    elif datatype == "morphology":
        model = models.get_model_commandline(the_config.models[0])
        log.info("Using the model specified in the .cfg file")
        command = "-f E -s %s %s -n fastTREE -p 123456789 %s" % (
            alignment_path, model, cmdline_extras)
    else:
        log.error("Unrecognised datatype: '%s'" % (datatype))
        raise util.PartitionFinderError


    # Force raxml to write to the dir with the alignment in it
    aln_dir, fname = os.path.split(alignment_path)
    command = ''.join([command, " -w '%s'" % os.path.abspath(aln_dir)])

    run_raxml_pthreads(command, cpus)
    alndir, aln = os.path.split(alignment_path)

    fast_tree_path = os.path.join(alndir, "RAxML_fastTree.fastTREE")

    # now we make the branch lengths with a partitioned model without rate multipliers
    if datatype == "DNA":
        log.info("Estimating GTR+G branch lengths on ML tree using all partitions")
        command = "-f e -s '%s' -t '%s' -m GTRGAMMA -O -n BLTREE -p 123456789 -q '%s' -w '%s' -e 1  " % (
            alignment_path, fast_tree_path, partition_file, os.path.abspath(alndir)) 
    elif datatype == "protein":
        log.info("Estimating LG+G branch lengths on ML tree using all partitions")
        command = "-f e -s '%s' -t '%s' -m PROTGAMMALG -O -n BLTREE -p 123456789 -q '%s' -w '%s' -e 1  " % (
            alignment_path, fast_tree_path, partition_file, os.path.abspath(alndir)) 
    elif datatype == "morphology":
        log.info("Estimating branch lengths on ML tree")
        command = "-f e -s '%s' -t '%s' %s -O -n BLTREE -p 123456789 -w '%s' -e 1  " % (
            alignment_path, fast_tree_path, model, os.path.abspath(alndir)) 
    else:
        log.error("Unrecognised datatype: '%s'" % (datatype))
        raise util.PartitionFinderError

    run_raxml_pthreads(command, cpus)
    tree_path = os.path.join(alndir, "RAxML_result.BLTREE")

    if not os.path.exists(tree_path):
        log.error("RAxML tree topology should be here but can't be be found: '%s'" % (tree_path))
        raise(util.PartitionFinderError)
    else:
        log.debug("RAxML tree with branch lengths ('%s') looks like this: ", tree_path)
        with open(tree_path, 'r') as fin:
            log.debug('%s', fin.read())

    log.info("ML topology estimation finished")

    return tree_path




def make_topology(alignment_path, datatype, cmdline_extras):
    '''Make a MP tree to start the analysis'''
    log.info("Making MP tree for %s", alignment_path)

    cmdline_extras = check_defaults(cmdline_extras)

    # First get the MP topology like this (-p is a hard-coded random number seed):
    if datatype == "DNA":
        command = "-y -s '%s' -m GTRGAMMA -n MPTREE -p 123456789 %s" % (
            alignment_path, cmdline_extras)
    elif datatype == "protein":
        command = "-y -s '%s' -m PROTGAMMALG -n MPTREE -p 123456789 %s" % (
            alignment_path, cmdline_extras)
    elif datatype == "morphology":
        # LOOK OUT: this relies on the assumption that we can only specify a single
        # model for morphology analyses...
        # choose a model for the data - necessary for RAxML to load the data properly
        model = models.get_model_commandline(the_config.models[0])
        command = "-y -s %s %s -K MK -n MPTREE -p 123456789 %s" % (
            alignment_path, model, cmdline_extras)
    else:
        log.error("Unrecognised datatype: '%s'" % (datatype))
        raise util.PartitionFinderError

    # Force raxml to write to the dir with the alignment in it
    aln_dir, fname = os.path.split(alignment_path)
    command = ''.join([command, " -w '%s'" % os.path.abspath(aln_dir)])

    run_raxml(command)
    dir, aln = os.path.split(alignment_path)
    tree_path = os.path.join(dir, "RAxML_parsimonyTree.MPTREE")

    if not os.path.exists(tree_path):
        log.error("RAxML tree topology should be here but can't be be found: '%s'" % (tree_path))
        raise(RaxmlError)
    else:
        log.debug("RAxML tree with branch lengths ('%s') looks like this: ", tree_path)
        with open(tree_path, 'r') as fin:
            log.debug('%s', fin.read())

    log.info("Topology estimation finished")

    return tree_path


def make_branch_lengths(alignment_path, topology_path, datatype, cmdline_extras):
    # Now we re-estimate branchlengths using a GTR+G model on the
    # (unpartitioned) dataset
    cmdline_extras = check_defaults(cmdline_extras)
    dir_path, fname = os.path.split(topology_path)
    tree_path = os.path.join(dir_path, 'topology_tree.phy')
    log.debug("Copying %s to %s", topology_path, tree_path)
    util.dupfile(topology_path, tree_path)
    os.remove(topology_path)  # saves headaches later...

    if datatype == "DNA":
        log.info("Estimating GTR+G branch lengths on tree using RAxML")
        command = "-f e -s '%s' -t '%s' -m GTRGAMMA -n BLTREE -w '%s' %s  " % (
            alignment_path, tree_path, os.path.abspath(dir_path), cmdline_extras)
    elif datatype == "protein":
        log.info("Estimating LG+G branch lengths on tree using RAxML")
        command = "-f e -s '%s' -t '%s' -m PROTGAMMALG -n BLTREE -w '%s' %s " % (
            alignment_path, tree_path, os.path.abspath(dir_path), cmdline_extras)
    elif datatype == "morphology":
        # LOOK OUT: this relies on the assumption that we can only specify a single
        # model for morphology analyses...
        # choose a model for the data - necessary for RAxML to load the data properly
        model = models.get_model_commandline(the_config.models[0])
        log.info("Estimating %s branch lengths on tree using RAxML", the_config.models[0])
        command = "-f e -s %s -t %s %s -K MK -n BLTREE -w %s %s" % (
                alignment_path, tree_path, model, os.path.abspath(dir_path), cmdline_extras)
    else:
        log.error("Unrecognised datatype: '%s'" % (datatype))
        raise util.PartitionFinderError

    run_raxml(command)


    dir, aln = os.path.split(alignment_path)
    tree_path = os.path.join(dir, "RAxML_result.BLTREE")

    if not os.path.exists(tree_path):
        log.error("RAxML tree topology should be here but can't be be found: '%s'" % (tree_path))
        raise util.PartitionFinderError
    else:
        log.debug("RAxML tree with branch lengths ('%s') looks like this: ", tree_path)
        with open(tree_path, 'r') as fin:
            log.debug('%s', fin.read())

    log.info("Branchlength estimation finished")

    # Now return the path of the final tree with branch lengths
    return tree_path


def check_defaults(cmdline_extras):
    # We use some sensible defaults, but allow users to override them with
    # extra cmdline options
    if cmdline_extras.count("-e") > 0:
        # then the user has specified a particular accuracy:
        accuracy = ""
    else:
        # we specify a default accuracy of 1 lnL unit
        accuracy = " -e 1.0 "

    # and we'll specify the -O option, so that the program doesn't exit if
    # there are undetermined seqs.  we'll put spaces at the start and end too,
    # just in case...
    cmdline_extras = ''.join([" ", cmdline_extras, accuracy, "-O "])
    return cmdline_extras


def analyse(model, alignment_path, tree_path, branchlengths, cmdline_extras):
    """Do the analysis -- this will overwrite stuff!"""

    # Move it to a new name to stop raxml stomping on different model analyses
    # dupfile(alignment_path, analysis_path)
    model_params = models.get_model_commandline(model)

    if branchlengths == 'linked':
        #constrain all branchlengths to be equal
        bl = ' -f B '
    elif branchlengths == 'unlinked':
        #let branchlenghts vary among subsets
        bl = ' -f e '
    else:
        # WTF?
        log.error("Unknown option for branchlengths: %s", branchlengths)
        raise util.PartitionFinderError

    cmdline_extras = check_defaults(cmdline_extras)

    # we can save memory on gappy alignments like this
    #if str(model).count('LG4')==0:
    #    cmdline_extras = ' '.join([cmdline_extras, '-U '])

    #raxml doesn't append alignment names automatically, like PhyML, let's do that here
    analysis_ID = raxml_analysis_ID(alignment_path, model)

    #force raxml to write to the dir with the alignment in it
    #-e 1.0 sets the precision to 1 lnL unit. This is all that's required here, and helps with speed.
    aln_dir, fname = os.path.split(alignment_path)
    command = " %s -s '%s' -t '%s' %s -n %s -w '%s' %s" % (
        bl, alignment_path, tree_path, model_params, analysis_ID, os.path.abspath(aln_dir), cmdline_extras)
    run_raxml(command)


def raxml_analysis_ID(alignment_path, model):
    dir, file = os.path.split(alignment_path)
    aln_name = os.path.splitext(file)[0]
    analysis_ID = '%s_%s.txt' % (aln_name, model)
    return analysis_ID


def make_tree_path(alignment_path):
    dir, aln = os.path.split(alignment_path)
    tree_path = os.path.join(dir, "RAxML_result.BLTREE")
    return tree_path


def make_output_path(alignment_path, model):
    analysis_ID = raxml_analysis_ID(alignment_path, model)
    dir, aln_file = os.path.split(alignment_path)
    stats_fname = "RAxML_info.%s" % (analysis_ID)
    stats_path = os.path.join(dir, stats_fname)
    tree_fname = "RAxML_result.%s" % (analysis_ID)
    tree_path = os.path.join(dir, tree_fname)
    return stats_path, tree_path


def remove_files(aln_path, model):
    '''remove all files from the alignment directory that are produced by raxml'''
    dir, file = os.path.split(aln_path)
    analysis_ID = raxml_analysis_ID(aln_path, model)
    dir = os.path.abspath(dir)
    fs = os.listdir(dir)
    fnames = fnmatch.filter(fs, '*%s*%s*' % ("RAxML", analysis_ID))

    pths = [os.path.join(dir, p) for p in fnames]
    util.delete_files(pths)


class RaxmlResult(DataRecord):
    pass


class Parser(object):
    def __init__(self, cfg):
        self.cfg = cfg

        if cfg.datatype == "protein":
            letters = _protein_letters
        elif cfg.datatype == "DNA":
            letters = _dna_letters
        elif cfg.datatype == "morphology":
            letters = "0123456789"
        else:
            log.error("Unknown datatype '%s', please check" % self.cfg.datatype)
            raise util.PartitionFinderError

        self.rate_indexes = self.cfg.data_layout.rate_indexes
        self.freq_indexes = self.cfg.data_layout.letter_indexes

        FLOAT = Word(nums + '.-').setParseAction(lambda x: float(x[0]))

        L = Word(letters, exact=1)
        COLON = Suppress(":")

        LNL_LABEL = Regex("Final GAMMA.+:") | Literal("Likelihood:")
        TIME_LABEL = Regex("Overall Time.+:") | Regex("Overall Time.+tion ")
        ALPHA_LABEL = Literal("alpha:")
        TREE_SIZE_LABEL = Literal("Tree-Length:")

        def labeled_float(label):
            return Suppress(SkipTo(label)) + Suppress(label) + FLOAT

        lnl = labeled_float(LNL_LABEL)
        lnl.setParseAction(self.set_lnl)

        seconds = labeled_float(TIME_LABEL)
        seconds.setParseAction(self.set_seconds)

        alpha = labeled_float(ALPHA_LABEL)
        alpha.setParseAction(self.set_alpha)

        tree_size = labeled_float(TREE_SIZE_LABEL)
        tree_size.setParseAction(self.set_tree_size)

        LG4X_LINE = "LG4X" + restOfLine
        lg4x = Optional(LG4X_LINE + LG4X_LINE)

        rate = Suppress("rate") + L + Suppress("<->") + L + COLON + FLOAT
        rate.setParseAction(self.set_rate)
        rates = OneOrMore(rate)

        freq = Suppress("freq pi(") + L + Suppress("):") + FLOAT
        freq.setParseAction(self.set_freq)
        freqs = OneOrMore(freq)

        LGM_LINE = "LGM" + restOfLine

        rate_block = Optional(LGM_LINE) + rates + freqs
        rate_block.setParseAction(self.rate_block)

        # Just look for these things
        self.root_parser = seconds + lnl + alpha + tree_size +\
            lg4x + OneOrMore(rate_block)

    def rate_block(self, tokens):
        self.current_block += 1

    def set_seconds(self, tokens):
        self.result.seconds = tokens[0]

    def set_lnl(self, tokens):
        self.result.lnl = tokens[0]

    def set_tree_size(self, tokens):
        self.result.site_rate = tokens[0]

    def set_alpha(self, tokens):
        self.result.alpha = tokens[0]

    def set_rate(self, tokens):
        # Ignore anything after the first rate block
        if self.current_block > 1:
            return
        basefrom, baseto, rate = tokens
        index = self.rate_indexes["%s_%s" % (basefrom, baseto)]
        self.result.rates[0, index] = rate

    def set_freq(self, tokens):
        # Ignore anything after the first rate block
        if self.current_block > 1:
            return
        base, rate = tokens
        index = self.freq_indexes[base]
        self.result.freqs[0, index] = rate

    def parse(self, text):
        log.debug("Parsing raxml output...")
        self.result = RaxmlResult(self.cfg)

        # The LGM produces multiple blocks of frequencies. We just want the
        # first one, so we must remember where we are.
        self.current_block = 1
        try:
            self.root_parser.parseString(text)
        except ParseException, p:
            log.error(str(p))
            raise util.ParseError

        log.debug("Result is %s", self.result)
        return self.result


def parse(text, cfg):
    the_parser = Parser(cfg)
    return the_parser.parse(text)


def fabricate(lnl):
    result = Parser('DNA')
    result.result = RaxmlResult()
    result.result.lnl = lnl
    result.result.tree_size = 0
    result.result.seconds = 0
    result.result.alpha = 0
    return result.result
