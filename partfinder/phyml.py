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
#program, the RAxML program, the PyParsing library, and the python-cluster library
#all of which are protected by their own licenses and conditions, using
#PartitionFinder implies that you agree with those licences and conditions as well.

"""Run phyml and parse the output"""

import logging
log = logging.getLogger("phyml")

import subprocess
import shlex
import os
import shutil
import sys
import util

from pyparsing import (
    Word, Literal, nums, Suppress, ParseException,
    SkipTo,
)

import phyml_models as models

_binary_name = 'phyml'
if sys.platform == 'win32':
    _binary_name += ".exe"

from util import PhylogenyProgramError


class PhymlError(PhylogenyProgramError):
    pass

def find_program():
    """Locate the binary ..."""
    pth = util.program_path
    pth = os.path.join(pth, _binary_name)

    log.debug("Checking for program %s in path %s", _binary_name, pth)
    if not os.path.exists(pth) or not os.path.isfile(pth):
        log.error("No such file: '%s'", pth)
        raise PhymlError
    log.debug("Found program %s at '%s'", _binary_name, pth)
    return pth

_phyml_binary = None

def run_phyml(command):
    global _phyml_binary
    if _phyml_binary is None:
        _phyml_binary = find_program()

    #turn off any memory checking in PhyML - thanks Jess Thomas for pointing out this problem
    command = "%s --no_memory_check" % (command)

    # Add in the command file
    log.debug("Running 'phyml %s'", command)
    command = "\"%s\" %s" % (_phyml_binary, command)

    # Note: We use shlex.split as it does a proper job of handling command
    # lines that are complex
    p = subprocess.Popen(
        shlex.split(command),
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)

    # Capture the output, we might put it into the errors
    stdout, stderr = p.communicate()
    # p.terminate()

    if p.returncode != 0:
        log.error("Phyml did not execute successfully")
        log.error("Phyml output follows, in case it's helpful for finding the problem")
        log.error("%s", stdout)
        log.error("%s", stderr)
        raise PhymlError


def dupfile(src, dst):
    # Make a copy or a symlink so that we don't overwrite different model runs
    # of the same alignment

    # TODO maybe this should throw...?
    try:
        if os.path.exists(dst):
            os.remove(dst)
        shutil.copyfile(src, dst)
    except OSError:
        log.error("Cannot link/copy file %s to %s", src, dst)
        raise PhymlError


def make_topology(alignment_path, datatype, cmdline_extras):
    '''Make a BioNJ tree to start the analysis'''
    log.info("Making BioNJ tree for %s", alignment_path)
    cmdline_extras = check_defaults(cmdline_extras)

    # First get the BioNJ topology like this:
    if datatype == "DNA":
        command = "-i '%s' -o n -b 0 %s" % (alignment_path, cmdline_extras)
    elif datatype == "protein":
        command = "-i '%s' -o n -b 0 -d aa %s" % (
            alignment_path, cmdline_extras)
    else:
        log.error("Unrecognised datatype: '%s'" % (datatype))
        raise(PhymlError)

    run_phyml(command)
    output_path = make_tree_path(alignment_path)
    return output_path


def make_branch_lengths(alignment_path, topology_path, datatype, cmdline_extras):
    # Now we re-estimate branchlengths using a GTR+I+G model on the (unpartitioned) dataset
    cmdline_extras = check_defaults(cmdline_extras)
    dir_path, fname = os.path.split(topology_path)
    tree_path = os.path.join(dir_path, 'topology_tree.phy')
    log.debug("Copying %s to %s", topology_path, tree_path)
    dupfile(topology_path, tree_path)

    if datatype == "DNA":
        log.info("Estimating GTR+I+G branch lengths on tree")
        command = "-i '%s' -u '%s' -m GTR -c 4 -a e -v e -o lr -b 0 %s" % (
            alignment_path, tree_path, cmdline_extras)
        run_phyml(command)
    if datatype == "protein":
        log.info("Estimating LG+F branch lengths on tree")
        command = "-i '%s' -u '%s' -m LG -c 1 -v 0 -f m -d aa -o lr -b 0 %s" % (
            alignment_path, tree_path, cmdline_extras)
        run_phyml(command)

    tree_path = make_tree_path(alignment_path)
    log.info("Branchlength estimation finished")

    # Now return the path of the final tree alignment
    return tree_path


def check_defaults(cmdline_extras):
    """We use some sensible defaults, but allow users to override them with extra cmdline options"""

    if cmdline_extras.count("--min_diff_lk_global") > 0:
        accuracy_global = ""
    else:
        accuracy_global = " --min_diff_lk_global 0.01 "
    if cmdline_extras.count("--min_diff_lk_local") > 0:
        accuracy_local = ""
    else:
        accuracy_local = " --min_diff_lk_local 0.01 "

    #we'll put spaces at the start and end too, just in case...
    cmdline_extras = ''.join(
        [" ", cmdline_extras, accuracy_local, accuracy_global, " "])

    return cmdline_extras


def analyse(model, alignment_path, tree_path, branchlengths, cmdline_extras):
    """Do the analysis -- this will overwrite stuff!"""

    # Move it to a new name to stop phyml stomping on different model analyses
    # dupfile(alignment_path, analysis_path)
    model_params = models.get_model_commandline(model)

    if branchlengths == 'linked':
        #constrain all branchlengths to be equal
        bl = ' --constrained_lens '
    elif branchlengths == 'unlinked':
        #let branchlenghts vary among subsets
        bl = ''
    else:
        # WTF?
        log.error("Unknown option for branchlengths: %s", branchlengths)
        raise PhymlError

    cmdline_extras = check_defaults(cmdline_extras)

    command = "--run_id %s -b 0 -i '%s' -u '%s' %s %s %s " % (
        model, alignment_path, tree_path, model_params, bl, cmdline_extras)
    run_phyml(command)

    # Now get rid of this -- we have the original elsewhere
    # os.remove(analysis_path)


def make_tree_path(alignment_path):
    pth, ext = os.path.splitext(alignment_path)
    return pth + ".phy_phyml_tree.txt"


def make_output_path(aln_path, model):
    # analyse_path = os.path.join(root_path, name + ".phy")
    pth, ext = os.path.splitext(aln_path)
    stats_path = "%s.phy_phyml_stats_%s.txt" % (pth, model)
    tree_path = "%s.phy_phyml_tree_%s.txt" % (pth, model)
    return stats_path, tree_path


def remove_files(aln_path, model):
    '''remove all files from the alignment directory that are produced by phyml'''
    fnames = make_output_path(aln_path, model)
    util.delete_files(fnames)


class PhymlResult(object):
    def __init__(self, lnl, tree_size, seconds):
        self.lnl = lnl
        self.seconds = seconds
        self.tree_size = tree_size
        self.rates = {} #placeholder, doesn't do anything but makes it compatible
        self.freqs = {} #placeholder, doesn't do anything but makes it compatible
        self.alpha = 0  #placeholder, doesn't do anything but makes it compatible

    def __str__(self):
        return "PhymlResult(lnl:%s, tree_size:%s, secs:%s)" % (self.lnl, self.tree_size, self.seconds)


class Parser(object):
    def __init__(self, datatype):
        FLOAT = Word(nums + '.-').setParseAction(lambda x: float(x[0]))
        INTEGER = Word(nums + '-').setParseAction(lambda x: int(x[0]))

        OB = Suppress("(")
        CB = Suppress(")")
        LNL_LABEL = Literal("Log-likelihood:")
        TREE_SIZE_LABEL = Literal("Tree size:")
        TIME_LABEL = Literal("Time used:")
        HMS = Word(nums + "hms")  # A bit rough...

        lnl = (LNL_LABEL + FLOAT("lnl"))
        tree_size = (TREE_SIZE_LABEL + FLOAT("tree_size"))
        time = (TIME_LABEL + HMS(
            "time") + OB + INTEGER("seconds") + Suppress("seconds") + CB)

        # Shorthand...
        def nextbit(label, val):
            return Suppress(SkipTo(label)) + val

        # Just look for these things
        self.root_parser = \
            nextbit(LNL_LABEL, lnl) +\
            nextbit(TREE_SIZE_LABEL, tree_size) +\
            nextbit(TIME_LABEL, time)

    def parse(self, text):
        log.debug("Parsing phyml output...")
        try:
            tokens = self.root_parser.parseString(text)
        except ParseException, p:
            log.error(str(p))
            raise PhymlError

        log.debug("Parsed LNL:      %s" % tokens.lnl)
        log.debug("Parsed TREESIZE: %s" % tokens.tree_size)
        log.debug("Parsed TIME:     %s" % tokens.time)

        return PhymlResult(lnl=tokens.lnl, tree_size=tokens.tree_size, seconds=tokens.seconds)


def parse(text, datatype):
    the_parser = Parser(datatype)
    return the_parser.parse(text)

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    import tempfile
    from alignment import TestAlignment
    import phyml_models
    alignment = TestAlignment("""
4 2208
spp1     CTTGAGGTTCAGAATGGTAATGAA------GTGCTGGTGCTGGAAGTTCAGCAGCAGCTCGGCGGCGGTATCGTACGTACCATCGCCATGGGTTCTTCCGACGGTCTGCGTCGCGGTCTGGATGTAAAAGACCTCGAGCACCCGATCGAAGTCCCAGTTGGTAAAGCAACACTGGGTCGTATCATGAACGTACTGGGTCAGCCAGTAGACATGAAGGGCGACATCGGTGAAGAAGAGCGTTGGGCT---------------ATCCACCGTGAAGCACCATCCTATGAAGAGCTGTCAAGCTCTCAGGAACTGCTGGAAACCGGCATCAAAGTTATCGACCTGATGTGTCCGTTTGCGAAGGGCGGTAAAGTTGGTCTGTTCGGTGGTGCGGGTGTAGGTAAAACCGTAAACATGATGGAGCTTATTCGTAACATCGCGATCGAGCACTCCGGTTATTCTGTGTTTGCGGGCGTAGGTGAACGTACTCGTGAGGGTAACGACTTCTACCACGAAATGACCGACTCCAACGTTATCGAT---------------------AAAGTTTCTCTGGTTTATGGCCAGATGAACGAGCCACCAGGTAACCGTCTGCGCGTTGCGCTGACCGGTCTGACCATGGCTGAGAAGTTCCGTGACGAAGGTCGCGACGTACTGCTGTTCGTCGATAACATCTATCGTTACACCCTGGCAGGTACTGAAGTTTCAGCACTGCTGGGTCGTATGCCTTCAGCGGTAGGTTACCAGCCGACTCTGGCGGAAGAAATGGGCGTTCGCATTCCAACGCTGGAAGAGTGTGATATCTGCCACGGCAGCGGCGCTAAAGCCGGTTCGAAGCCGCAGACCTGTCCTACCTGTCACGGTGCAGGCCAGGTACAGATGCGCCAGGGCTTCTTCGCTGTACAGCAGACCTGTCCACACTGCCAGGGCCGCGGTACGCTGATCAAAGATCCGTGCAACAAATGTCACGGTCATGGTCGCGTAGAGAAAACCAAAACCCTGTCCGTAAAAATTCCGGCAGGCGTTGATACCGGCGATCGTATTCGTCTGACTGGCGAAGGTGAAGCTGGTGAGCACGGCGCACCGGCAGGCGATCTGTACGTTCAGGTGCAGGTGAAGCAGCACGCTATTTTCGAGCGTGAAGGCAACAACCTGTACTGTGAAGTGCCGATCAACTTCTCAATGGCGGCTCTTGGCGGCGAGATTGAAGTGCCGACGCTTGATGGTCGCGTGAAGCTGAAAGTTCCGGGCGAAACGCAAACTGGCAAGCTGTTCCGTATGCGTGGCAAGGGCGTGAAGTCCGTGCGCGGCGGTGCACAGGGCGACCTTCTGTGCCGCGTGGTGGTCGAGACACCGGTAGGTCTTAACGAGAAGCAGAAACAGCTGCTCAAAGATCTGCAGGAAAGTTTTGGCGGCCCAACGGGTGAAAACAACGTTGTTAACGCCCTGTCGCAGAAACTGGAATTGCTGATCCGCCGCGAAGGCAAAGTACATCAGCAAACTTATGTCCATGGTGTGCCACAGGCTCCGCTGGCGGTAACCGGTGAAACGGAAGTGACCGGTACACAGGTGCGTTTCTGGCCAAGCCACGAAACCTTCACCAACGTAATCGAATTCGAATATGAGATTCTGGCAAAACGTCTGCGCGAGCTGTCATTCCTGAACTCCGGCGTTTCCATCCGTCTGCGCGATAAGCGTGAC---GGCAAAGAAGACCATTTCCACTATGAAGGTGGTATCAAGGCGTTTATTGAGTATCTCAATAAAAATAAAACGCCTATCCACCCGAATATCTTCTACTTCTCCACCGAA---AAAGACGGTATTGGCGTAGAAGTGGCGTTGCAGTGGAACGATGGTTTCCAGGAAAACATCTACTGCTTCACCAACAACATTCCACAGCGTGATGGCGGTACTCACCTTGCAGGCTTCCGTGCGGCGATGACCCGTACGCTGAACGCTTACATGGACAAAGAAGGCTACAGCAAAAAAGCCAAA------GTCAGCGCCACCGGTGATGATGCCCGTGAAGGCCTGATTGCCGTCGTTTCCGTGAAAGTACCGGATCCGAAATTCTCCTCTCAGACTAAAGACAAACTGGTCTCTTCTGAGGTGAAAACGGCGGTAGAACAGCAGATGAATGAACTGCTGAGCGAATACCTGCTGGAAAACCCGTCTGACGCCAAAATC
spp2     CTTGAGGTACAAAATGGTAATGAG------AGCCTGGTGCTGGAAGTTCAGCAGCAGCTCGGTGGTGGTATCGTACGTGCTATCGCCATGGGTTCTTCCGACGGTCTGCGTCGTGGTCTGGAAGTTAAAGACCTTGAGCACCCGATCGAAGTCCCGGTTGGTAAAGCAACGCTGGGTCGTATCATGAACGTGCTGGGTCAGCCGATCGATATGAAAGGCGACATCGGCGAAGAAGAACGTTGGGCG---------------ATTCACCGTGCAGCACCTTCCTATGAAGAGCTCTCCAGCTCTCAGGAACTGCTGGAAACCGGCATCAAAGTTATCGACCTGATGTGTCCGTTCGCGAAGGGCGGTAAAGTCGGTCTGTTCGGTGGTGCGGGTGTTGGTAAAACCGTAAACATGATGGAGCTGATCCGTAACATCGCGATCGAACACTCCGGTTACTCCGTGTTTGCTGGTGTTGGTGAGCGTACTCGTGAGGGTAACGACTTCTACCACGAAATGACCGACTCCAACGTTCTGGAT---------------------AAAGTATCCCTGGTTTACGGCCAGATGAACGAGCCGCCGGGAAACCGTCTGCGCGTTGCACTGACCGGCCTGACCATGGCTGAGAAATTCCGTGACGAAGGTCGTGACGTTCTGCTGTTCGTCGATAACATCTATCGTTATACCCTGGCCGGTACAGAAGTATCTGCACTGCTGGGTCGTATGCCTTCTGCGGTAGGTTATCAGCCGACGCTGGCGGAAGAGATGGGCGTTCGTATCCCGACGCTGGAAGAGTGCGACGTCTGCCACGGCAGCGGCGCGAAATCTGGCAGCAAACCGCAGACCTGTCCGACCTGTCATGGTCAGGGCCAGGTGCAGATGCGTCAGGGCTTCTTCGCCGTTCAGCAGACCTGTCCGCATTGTCAGGGGCGCGGTACGCTGATTAAAGATCCGTGCAACAAATGTCACGGTCACGGTCGCGTTGAGAAAACCAAAACCCTGTCGGTCAAAATCCCGGCGGGCGTGGATACCGGCGATCGTATTCGTCTGTCAGGAGAAGGCGAAGCGGGCGAACACGGTGCACCAGCAGGCGATCTGTACGTTCAGGTCCAGGTTAAGCAGCACGCCATCTTTGAGCGTGAAGGCAATAACCTGTACTGCGAAGTGCCTATTAACTTCACCATGGCAGCCCTCGGCGGCGAGATTGAAGTCCCGACGCTGGATGGCCGGGTGAATCTCAAAGTGCCTGGCGAAACGCAAACCGGCAAACTGTTCCGCATGCGCGGTAAAGGTGTGAAATCCGTGCGCGGTGGTGCTCAGGGCGACCTGCTGTGCCGCGTGGTGGTTGAAACACCAGTCGGGCTGAACGATAAGCAGAAACAGCTGCTGAAGGACCTGCAGGAAAGTTTTGGCGGACCAACGGGCGAGAAAAACGTGGTTAACGCCCTGTCGCAGAAGCTGGAGCTGGTTATTCAGCGCGACAATAAAGTTCACCGTCAGATCTATGCGCACGGTGTGCCGCAGGCTCCGCTGGCAGTGACCGGTGAGACCGAAAAAACCGGCACCATGGTACGTTTCTGGCCAAGCTATGAAACCTTCACCAACGTTGTCGAGTTCGAATACGAGATCCTGGCAAAACGTCTGCGTGAGCTGTCGTTCCTGAACTCCGGGGTTTCTATCCGTCTGCGTGACAAGCGTGAC---GGTAAAGAAGACCATTTCCACTACGAAGGCGGCATCAAGGCGTTCGTTGAGTATCTCAATAAGAACAAAACGCCGATCCACCCGAATATCTTCTACTTCTCCACCGAA---AAAGACGGTATTGGCGTCGAAGTAGCGCTGCAGTGGAACGACGGCTTCCAGGAAAACATCTACTGCTTCACCAACAACATCCCGCAGCGCGATGGCGGTACTCACCTTGCGGGCTTCCGCGCGGCGATGACCCGTACCCTGAACGCCTATATGGACAAAGAAGGCTACAGCAAAAAAGCCAAA------GTCAGCGCTACCGGCGACGATGCGCGTGAAGGCCTGATTGCCGTTGTCTCCGTGAAGGTTCCGGATCCGAAATTCTCCTCGCAGACCAAAGACAAACTGGTCTCCTCCGAGGTGAAAACCGCGGTTGAACAGCAGATGAATGAACTGCTGAACGAATACCTGCTGGAAAATCCGTCTGACGCGAAAATC
spp3     CTTGAGGTACAGAATAACAGCGAG------AAGCTGGTGCTGGAAGTTCAGCAGCAGCTCGGCGGCGGTATCGTACGTACCATCGCAATGGGTTCTTCCGACGGTCTGCGTCGTGGTCTGGAAGTGAAAGACCTCGAGCACCCGATCGAAGTCCCGGTAGGTAAAGCGACCCTGGGTCGTATCATGAACGTGCTGGGTCAGCCAATCGATATGAAAGGCGACATCGGCGAAGAAGATCGTTGGGCG---------------ATTCACCGCGCAGCACCTTCCTATGAAGAGCTGTCCAGCTCTCAGGAACTGCTGGAAACCGGCATCAAAGTTATCGACCTGATTTGTCCGTTCGCTAAGGGCGGTAAAGTTGGTCTGTTCGGTGGTGCGGGCGTAGGTAAAACCGTAAACATGATGGAGCTGATCCGTAACATCGCGATCGAGCACTCCGGTTACTCCGTGTTTGCAGGCGTGGGTGAGCGTACTCGTGAGGGTAACGACTTCTACCACGAGATGACCGACTCCAACGTTCTGGAC---------------------AAAGTTGCACTGGTTTACGGCCAGATGAACGAGCCGCCAGGTAACCGTCTGCGCGTAGCGCTGACCGGTCTGACCATCGCGGAGAAATTCCGTGACGAAGGCCGTGACGTTCTGCTGTTCGTCGATAACATCTATCGTTATACCCTGGCCGGTACAGAAGTTTCTGCACTGCTGGGTCGTATGCCATCTGCGGTAGGTTATCAGCCTACTCTGGCAGAAGAGATGGGTGTTCGTATCCCGACGCTGGAAGAGTGTGAAGTTTGCCACGGCAGCGGCGCGAAAAAAGGTTCTTCTCCGCAGACCTGTCCAACCTGTCATGGACAGGGCCAGGTGCAGATGCGTCAGGGCTTCTTCACCGTGCAGCAAAGCTGCCCGCACTGCCAGGGCCGCGGTACCATCATTAAAGATCCGTGCACCAACTGTCACGGCCATGGCCGCGTAGAGAAAACCAAAACGCTGTCGGTAAAAATTCCGGCAGGCGTGGATACCGGCGATCGTATCCGCCTTTCTGGTGAAGGCGAAGCGGGCGAGCACGGCGCACCTTCAGGCGATCTGTACGTTCAGGTTCAGGTGAAACAGCACCCAATCTTCGAGCGTGAAGGCAATAACCTGTACTGCGAAGTGCCGATCAACTTTGCGATGGCTGCGCTGGGCGGGGAAATTGAAGTGCCGACCCTTGACGGCCGCGTTAAGCTGAAGGTACCGAGCGAAACGCAAACCGGCAAGCTGTTCCGCATGCGCGGTAAAGGCGTGAAATCCGTACGCGGTGGCGCGCAGGGCGATCTGCTGTGCCGCGTCGTCGTTGAAACTCCGGTTAGCCTGAACGAAAAGCAGAAGAAACTGCTGCGTGATTTGGAAGAGAGCTTTGGCGGCCCAACGGGGGCGAACAATGTTGTGAACGCCCTGTCCCAGAAGCTGGAGCTGCTGATTCGCCGCGAAGGCAAAACCCATCAGCAAACCTACGTGCACGGTGTGCCGCAGGCTCCGCTGGCGGTCACCGGTGAAACCGAACTGACCGGTACCCAGGTGCGTTTCTGGCCGAGCCATGAAACCTTCACCAACGTCACCGAATTCGAATATGACATCCTGGCTAAGCGCCTGCGTGAGCTGTCGTTCCTGAACTCCGGCGTCTCTATTCGCCTGAACGATAAGCGCGAC---GGCAAGCAGGATCACTTCCACTACGAAGGCGGCATCAAGGCGTTTGTTGAGTACCTCAACAAGAACAAAACCCCGATTCACCCGAACGTCTTCTATTTCAGCACTGAA---AAAGACGGCATCGGCGTGGAAGTGGCGCTGCAGTGGAACGACGGCTTCCAGGAAAATATCTACTGCTTTACCAACAACATTCCTCAGCGCGACGGCGGTACTCACCTTGCGGGCTTCCGCGCGGCGATGACCCGTACCCTGAACGCCTATATGGACAAAGAAGGCTACAGCAAAAAAGCCAAA------GTGAGCGCCACCGGTGACGATGCGCGTGAAGGCCTGATTGCCGTAGTGTCCGTGAAGGTGCCGGATCCGAAGTTCTCTTCCCAGACCAAAGACAAACTGGTTTCTTCGGAAGTGAAATCCGCGGTTGAACAGCAGATGAACGAACTGCTGGCTGAATACCTGCTGGAAAATCCGGGCGACGCAAAAATT
spp4     CTCGAGGTGAAAAATGGTGATGCT------CGTCTGGTGCTGGAAGTTCAGCAGCAGCTGGGTGGTGGCGTGGTTCGTACCATCGCCATGGGTACTTCTGACGGCCTGAAGCGCGGTCTGGAAGTTACCGACCTGAAAAAACCTATCCAGGTTCCGGTTGGTAAAGCAACCCTCGGCCGTATCATGAACGTATTGGGTGAGCCAATCGACATGAAAGGCGACCTGCAGAATGACGACGGCACTGTAGTAGAGGTTTCCTCTATTCACCGTGCAGCACCTTCGTATGAAGATCAGTCTAACTCGCAGGAACTGCTGGAAACCGGCATCAAGGTTATCGACCTGATGTGTCCGTTCGCTAAGGGCGGTAAAGTCGGTCTGTTCGGTGGTGCGGGTGTAGGTAAAACCGTAAACATGATGGAGCTGATCCGTAACATCGCGGCTGAGCACTCAGGTTATTCGGTATTTGCTGGTGTGGGTGAGCGTACTCGTGAGGGTAACGACTTCTACCACGAAATGACTGACTCCAACGTTATCGAT---------------------AAAGTAGCGCTGGTGTATGGCCAGATGAACGAGCCGCCGGGTAACCGTCTGCGCGTAGCACTGACCGGTTTGACCATGGCGGAAAAATTCCGTGATGAAGGCCGTGACGTTCTGCTGTTCATCGACAACATCTATCGTTACACCCTGGCCGGTACTGAAGTATCAGCACTGCTGGGTCGTATGCCATCTGCGGTAGGCTATCAGCCAACGCTGGCAGAAGAGATGGGTGTGCGCATTCCAACACTGGAAGAGTGCGATGTCTGCCACGGTAGCGGCGCGAAAGCGGGGACCAAACCGCAGACCTGTCATACCTGTCATGGCGCAGGCCAGGTGCAGATGCGTCAGGGCTTCTTCACTGTGCAGCAGGCGTGTCCGACCTGTCACGGTCGCGGTTCAGTGATCAAAGATCCGTGCAATGCTTGTCATGGTCACGGTCGCGTTGAGCGCAGTAAAACCCTGTCGGTGAAAATTCCAGCAGGCGTGGATACCGGCGATCGCATTCGTCTGACCGGCGAAGGTGAAGCGGGCGAACAGGGCGCACCAGCGGGCGATCTGTACGTTCAGGTTTCGGTGAAAAAGCACCCGATCTTTGAGCGTGAAGATAACAACCTATATTGCGAAGTGCCGATTAACTTTGCGATGGCAGCATTGGGTGGCGAGATTGAAGTGCCGACGCTTGATGGGCGTGTGAACCTGAAAGTGCCTTCTGAAACGCAAACTGGCAAGCTGTTCCGCATGCGCGGTAAAGGCGTGAAATCGGTGCGTGGTGGTGCGGTAGGCGATTTGCTGTGTCGTGTGGTGGTGGAAACGCCAGTTAGCCTCAATGACAAACAGAAAGCGTTACTGCGTGAACTGGAAGAGAGTTTTGGCGGCCCGAGCGGTGAGAAAAACGTCGTAAACGCCCTGTCACAGAAGCTGGAGCTGACCATTCGCCGTGAAGGCAAAGTGCATCAGCAGGTTTATCAGCACGGCGTGCCGCAGGCACCGCTGGCGGTGTCCGGTGATACCGATGCAACCGGTACTCGCGTGCGTTTCTGGCCGAGCTACGAAACCTTCACCAATGTGATTGAGTTTGAGTACGAAATCCTGGCGAAACGCCTGCGTGAACTGTCGTTCCTGAACTCTGGCGTTTCGATTCGTCTGGAAGACAAACGCGAC---GGCAAGAACGATCACTTCCACTACGAAGGCGGCATCAAGGCGTTCGTTGAGTATCTCAACAAGAACAAAACCCCGATTCACCCAACGGTGTTCTACTTCTCGACGGAG---AAAGATGGCATTGGCGTGGAAGTGGCGCTGCAGTGGAACGATGGTTTCCAGGAAAACATCTACTGCTTCACCAACAACATTCCACAGCGCGACGGCGGTACGCACCTGGCGGGCTTCCGTGCGGCAATGACGCGTACGCTGAATGCCTACATGGATAAAGAAGGCTACAGCAAAAAAGCCAAA------GTCAGTGCGACCGGTGACGATGCGCGTGAAGGCCTGATTGCAGTGGTTTCCGTGAAAGTGCCGGATCCGAAATTCTCTTCTCAGACCAAAGATAAGCTGGTCTCTTCTGAAGTGAAATCGGCGGTTGAGCAGCAGATGAACGAACTGCTGGCGGAATACCTGCTGGAAAATCCGTCTGACGCGAAAATC
""")
    tmp = tempfile.mkdtemp()
    pth = os.path.join(tmp, 'test.phy')
    alignment.write(pth)
    tree_path = make_tree(pth)
    log.info("Tree is %s:", open(tree_path).read())

    for model in phyml_models.get_all_models():
        log.info("Analysing using model %s:" % model)
        out_pth = analyse(model, pth, tree_path)
        output = open(out_pth, 'rb').read()
        res = parse(output)
        log.info("Result is %s", res)

    # shutil.rmtree(tmp)
