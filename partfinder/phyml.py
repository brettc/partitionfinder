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
import csv

from pyparsing import (
    Word, Literal, nums, Suppress, ParseException,
    SkipTo,
)
from math import log as logarithm

import phyml_models as models

_binary_name = 'phyml'
if sys.platform == 'win32':
    _binary_name += ".exe"

from util import PhylogenyProgramError


class PhymlError(PhylogenyProgramError):
    def __init__(self, stderr, stdout):
        self.stderr = stderr
        self.stdout = stdout

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

def run_phyml(command, report_errors=True):
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
        if report_errors == True:
            log.error("Phyml did not execute successfully")
            log.error("Phyml output follows, in case it's helpful for finding the problem")
            log.error("%s", stdout)
            log.error("%s", stderr)
        raise PhymlError(stdout, stderr)


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

def likelihood_parser(phyml_lk_file):
    '''
    Takes a *_phyml_lk.txt file and returns a dictionary of sites and site
    likelihoods and a dictionary of sites and lists of likelihoods under
    different rate categories. If no rate categories are specified, it will
    return a dictionary with sites and likelihoods P(D|M) for each site.

    Here is an example of the first few lines of the file that it takes:

    Note : P(D|M) is the probability of site D given the model M (i.e., the
    site likelihood) P(D|M,rr[x]) is the probability of site D given the model
    M and the relative rate of evolution rr[x], where x is the class of rate to
    be considered.  We have P(D|M) = \sum_x P(x) x P(D|M,rr[x]).

    Site   P(D|M)          P(D|M,rr[1]=2.6534)   P(D|M,rr[2]=0.2289)   P(D|M,rr[3]=0.4957)   P(D|M,rr[4]=1.0697)   Posterior mean
    1      2.07027e-12     1.3895e-19            6.2676e-12            1.2534e-12            1.21786e-15           0.273422
    2      1.8652e-07      2.05811e-19           6.73481e-07           4.14575e-09           7.97623e-14           0.23049
    3      4.48873e-15     1.37274e-19           7.11221e-15           9.11826e-15           9.21848e-17           0.382265
    4      3.38958e-10     1.31413e-19           1.18939e-09           4.20659e-11           5.86537e-15           0.237972
    5      8.29969e-17     1.11587e-19           3.1672e-17            2.52183e-16           1.9722e-17            0.502077
    6      9.24579e-09     1.59891e-19           3.31101e-08           4.79946e-10           2.59524e-14           0.232669
    7      3.43996e-10     2.1917e-19            1.19544e-09           5.43128e-11           1.22969e-14           0.240455
    8      4.43262e-13     1.1447e-19            1.32148e-12           2.8874e-13            3.7386e-16            0.27685
    9      3.42513e-11     1.70149e-19           1.14227e-10           1.02103e-11           4.05239e-15           0.250765
    10     1.15506e-11     1.28107e-19           3.86378e-11           3.32642e-12           1.46151e-15           0.250024
    '''
    try:
        with open(str(phyml_lk_file)) as phyml_lk_file:
            # The phyml_lk files differ based on whether different rate
            # categories are estimated or not, this figures out which
            # file we are dealing with
            phyml_lk_file.next()
            line2 = phyml_lk_file.next()
            # Check to see if the file contains rate categories
            if line2[0] != "P":
                phyml_lk_file.next()

            # If it contains rate categories, we need to skip a few more lines
            else:
                for _ in xrange(4):
                    phyml_lk_file.next()
            # Read in the contents of the file and get rid of whitespace
            list_of_dicts = list(csv.DictReader(phyml_lk_file,
                delimiter = " ", skipinitialspace = True))
    except IOError:
        raise IOError("Could not find the likelihood file!")
    phyml_lk_file.close()

    # Right now, when the alignment is over 1,000,000 sites, PhyML
    # merges the site number with the site likelihood, catch that and
    # throw an error
    if len(list_of_dicts) > 999999:
        raise IOError("PhyML file cannot process more than 1 M sites")

    # The headers values change with each run so we need a list of them
    headers = []
    for k in list_of_dicts[0]:
        headers.append(k)
    # Sort the headers into alphabetical order
    headers.sort()

    # Check if the rate cateogories were estimated, if they weren't
    # just return the likelihood scores for each site, otherwise, return
    # site likelihoods and likelihoods under each rate category
    if len(headers) < 4:
        # Make a list of site log likelihoods
        likelihood_list = [[logarithm(float(site[headers[1]]))] for site in list_of_dicts]
        return likelihood_list

    else:
        # Make a list of site log likelihoods
        if list_of_dicts[0][headers[1]] == 'nan' or list_of_dicts[0][headers[1]] == 'inf':
            likelihood_list = None
            print "Whoopsies!"
            rate_list = None
            lk_rate_list = None
            lk_site_rate_list = None
            return likelihood_list, lk_rate_list, rate_list, lk_site_rate_list

        else:
            likelihood_list = [[logarithm(float(site[headers[1]]))] for site in list_of_dicts]

            # Make a rate list
            # print list_of_dicts[0][headers[len(headers) - 3]]
            # if list_of_dicts[0][headers[len(headers) - 3]] == 'nan' or list_of_dicts[0][headers[len(headers) - 3]] == 'inf':
            #     rate_list = None
            #     print "Whoopsies!"
            # else:
            rate_list = [[(logarithm(float(site[headers[len(headers) - 3]])))] for site in list_of_dicts]

            # Now make a list of lists of site likelihoods under different
            # rate categories
            lk_rate_list = []
            for i in list_of_dicts:
                ind_lk_list = []
                # Pull the likelihood from each rate category by calling the
                # appropriate key from "headers"
                for num in range(2, len(headers) - 3):
                    ind_lk_list.append(logarithm(float(i[headers[num]])))
                # Now add the list of likelihoods for the site to a master list
                lk_rate_list.append(ind_lk_list)

            # Now pull likelihoods and rates for a two dimensional list
            lk_site_rate_list = []
            for i in list_of_dicts:
                ind_lk_r_list = []
                ind_lk_r_list.append(logarithm(float(i[headers[1]])))
                ind_lk_r_list.append(logarithm(float(i[headers[len(headers) - 3]])))
                lk_site_rate_list.append(ind_lk_r_list)
            # Return both the list of site likelihoods and the list of lists of
            # likelihoods under different rate categories
            return likelihood_list, lk_rate_list, rate_list, lk_site_rate_list

program_name = "phyml"

def program():
    return program_name

def gen_per_site_stats(cfg, alignment_path, tree_path):
    if cfg.datatype == 'DNA':
        if cfg.branchlengths == 'linked':
            #command = "--run_id GTRGAMMA -b 0 -i '%s' -u '%s' -m JC69 --print_site_lnl --constrained_lens -f 0.25,0.25,0.25,0.25" % (alignment_path, tree_path)
            command = "--run_id GTRGAMMA -b 0 -i '%s' -u '%s' -m GTR --print_site_lnl --constrained_lens" % (alignment_path, tree_path)
        elif cfg.branchlengths == 'unlinked':
            command = "--run_id GTRGAMMA -b 0 -i '%s' -u '%s' -m GTR --print_site_lnl" % (
                alignment_path, tree_path)

    elif cfg.datatype == 'protein':
        if cfg.branchlengths == 'linked':
            command = "--run_id LGGAMMA -b 0 -d aa -i '%s' -u '%s' -m LG --print_site_lnl --constrained_lens" % (
                alignment_path, tree_path)
        elif cfg.branchlengths == 'unlinked':
            command = "--run_id LGGAMMA -b 0 -d aa -i '%s' -u '%s' -m LG --print_site_lnl" % (
                alignment_path, tree_path)
    # when we run phyml we don't report errors, because we don't want to clutter the output
    # the errors are still raised though. This is particular to running the kmeans output.
    run_phyml(command, report_errors=False)

def get_per_site_stats(phylip_file, cfg):
    # Retreive a list of the site likelihoods
    if cfg.datatype == 'DNA':
        phyml_lk_fname = ("%s_phyml_lk_GTRGAMMA.txt" % phylip_file)
    elif cfg.datatype == 'protein':
        phyml_lk_fname = ("%s_phyml_lk_LGGAMMA.txt" % phylip_file)
    # Open the phyml output and return the list of four lists (likelihoods,
    # rates, likelihoods_and_rates, and likelihoods under different rate
    # categories)
    return likelihood_parser(phyml_lk_fname)

def fabricate(lnl):
    result = PhymlResult(lnl, 0, 0)
    return result

def get_CIs(cfg):
    ci_list = []
    fname = os.path.join(cfg.base_path, 'rates.txt')
    the_cis = open(fname)
    for ci in the_cis.readlines():
        ci_list.append([logarithm(float(ci))])
    return ci_list

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
