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

import hashlib
import cPickle as pickle
import subset
from util import get_aic, get_aicc, get_bic
from scipy.stats import chi2 
from util import PartitionFinderError

class AnalysisError(PartitionFinderError):
    pass



def columnset_to_string(colset):
    s = list(colset)
    s.sort()
    # Add one, cos we converted to zero base...
    return ', '.join([str(x+1) for x in s])

def subset_unique_name(columns):
    """Return a unique string based on the subsets columns (which are unique)"""

    # Use pickle to give us a string from the object
    pickled_columns = pickle.dumps(columns, -1)

    # Now get an md5 hash from this. There is some vanishingly small chance that
    # we'll get the same thing. Google "MD5 Hash Collision"
    return hashlib.md5(pickled_columns).hexdigest()

def merge_fabricated_subsets(subset_list):
    '''Allows the merging of fabricated subsets and the preservation of their
    centroids and lnls'''
    columns = set()
    lnl = 0
    centroid = []

    # Figure out how many dimensions the centroid is
    centroid_dim = len(subset_list[0].centroid)
    for i in range(centroid_dim):
        centroid.append(0)

    for sub in subset_list:
        columns |= sub.column_set
        lnl += sub.best_lnl
        number = 0
        for observation in centroid:
            observation += sub.centroid[number]
            number += 1

    # Now just take the average of each centroid to be the centroid of the new
    # subset
    centroid = [x/len(subset_list) for x in centroid]

    new_sub = subset.Subset(sub.cfg, columns)

    # Add the centroid and sum of the lnls to the subset. TODO: create
    # functions to update these variables in the subset rather than messing
    # with them directly
    new_sub.centroid = centroid
    new_sub.lnl = lnl
    return new_sub


def merge_subsets(subset_list):
    """Take a set of subsets and merge them together"""
    columns = set()

    # We just need the columns
    names = []
    descriptions = []
    for sub in subset_list:
        columns |= sub.column_set
        descriptions.extend(sub.description)
        names.extend(sub.names)

    newsub = subset.Subset(sub.cfg, columns)
    # Only add the description if it isn't there (we might get back a cache
    # hit)
    if not newsub.names:
        newsub.add_description(names, descriptions)

    return newsub

def subsets_overlap(subset_list):
    columns = set()
    overlapping = []

    for sub in subset_list:
        # If the intersection is non-empty...
        ov = list(sub.column_set & columns)
        if ov:
            overlapping.append(ov)
        columns |= sub.column_set

    return ov

def check_against_alignment(full_subset, alignment, the_config):
    """Check the subset definition against the alignment"""

    alignment_set = set(range(0, alignment.sequence_length))
    leftout = alignment_set - full_subset.column_set
    if leftout:
        log.warning(
            "These columns are missing from the block definitions: %s",
            columnset_to_string(leftout))
        if the_config.no_ml_tree == False:
            log.error(
                "You cannot estimate a Maximum Likelihood (ML) starting tree"
                " (the default behaviour) when you have columns missing from"
                " your data block definitions, because the method we use "
                "to estimate the ML tree requires all sites in the alignment"
                " to be assigned to a data block. We recommend that you "
                " either remove the sites you don't want from your alignment"
                " or (if possible) include the missing sites in appropriate"
                " data blocks. Failing that, you can use the --no-ml-tree "
                " command line option. In this case, a NJ (PhyML) or MP"
                "(RaxML) starting tree will be estimated for your analysis. "
            )
            raise AnalysisError


def split_subset(a_subset, cluster_list):
    """Takes a subset and splits it according to a cluster list,
     then returns the subsets resulting from the split"""
    # Take each site from the first list and add it to a new
    subset_list = a_subset.columns
    subset_columns = []
    list_of_subsets = []
    for cluster in cluster_list:
        list_of_sites = []
        for site in cluster:
            list_of_sites.append(subset_list[site - 1])
        subset_columns.append(set(list_of_sites))

    tracker = 0
    for column_set in subset_columns:
        new_subset = subset.Subset(a_subset.cfg, column_set)
        list_of_subsets.append(new_subset)
        tracker += 1

    return list_of_subsets

def subset_list_score(list_of_subsets, the_config, alignment):
    """Takes a list of subsets and return the aic, aicc, or bic score"""

    lnL, sum_k, subs_len = subset_list_stats(list_of_subsets, the_config, alignment)

    if the_config.model_selection == 'aic':
        return get_aic(lnL, sum_k)
    elif the_config.model_selection == 'aicc':
        return get_aicc(lnL, sum_k, subs_len)
    elif the_config.model_selection == 'bic':
        return get_bic(lnL, sum_k, subs_len)


def subset_list_stats(list_of_subsets, the_config, alignment):
    """Takes a list of subsets and returns the lnL and the number of params"""
    sum_subset_k = 0
    lnL = 0
    subs_len = 0
    for sub in list_of_subsets:
        sum_subset_k += sub.best_params
        lnL += sub.best_lnl
        subs_len += len(sub.columns)
    # Grab the number of species so we know how many params there are
    num_taxa = len(alignment.species)
    # Linked brlens - only one extra parameter per subset
    if the_config.branchlengths == 'linked':
        sum_k = sum_subset_k + (len(list_of_subsets) - 1) + (
            (2 * num_taxa) - 3)
        log.debug("Total parameters from brlens: %d" %
                  ((2 * num_taxa) - 3))
        log.debug("Parameters from subset multipliers: %d" %
                  (len(list_of_subsets) -1))

    # Unlinked brlens - every subset has its own set of brlens
    elif the_config.branchlengths == 'unlinked':
        sum_k = sum_subset_k + (len(list_of_subsets) * (
            (2 * num_taxa) - 3))
        log.debug("Total parameters from brlens: %d" % ((
            2 * num_taxa) - 3) * (len(list_of_subsets)))

    log.debug("Grand_total_parameters: %d", sum_k)

    return lnL, sum_k, subs_len


def subset_list_score_diff(list1, list2, the_config, alignment):
    """Take two lists of subsets and return the score diff as list1 - list2"""
    list1_score = subset_list_score(list1, the_config, alignment)
    list2_score = subset_list_score(list2, the_config, alignment)

    score_diff = list1_score - list2_score

    return score_diff