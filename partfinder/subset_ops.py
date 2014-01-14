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

import hashlib
import cPickle as pickle
import subset
from util import get_aic, get_aicc, get_bic

import logging
log = logging.getLogger("subset_ops")

def subset_unique_name(columns):
    """Return a unique string based on the subsets columns (which are unique)"""

    # Use pickle to give us a string from the object
    pickled_columns = pickle.dumps(columns, -1)

    # Now get an md5 hash from this. There is some vanishingly small chance that
    # we'll get the same thing. Google "MD5 Hash Collision"
    return hashlib.md5(pickled_columns).hexdigest()


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

    for sub in subset_list:
        # If the intersection is non-empty...
        if sub.column_set & columns:
            return True
        columns |= sub.column_set

    return False


def has_missing(subset_list):
    return False


def split_subset(a_subset, cluster_list):
    """Takes a subset and splits it according to a cluster list,
     then returns the subsets resulting from the split"""
    # Take each site from the first list and add it to a new
    subset_list = a_subset.columns
    likelihood_list = a_subset.site_lnls_GTRG
    subset_columns = []
    site_likelihoods = []
    list_of_subsets = []
    for cluster in cluster_list:
        list_of_sites = []
        likelihood_for_site = []
        for site in cluster:
            list_of_sites.append(subset_list[site - 1])
            likelihood_for_site += (likelihood_list[site - 1])
        subset_columns.append(set(list_of_sites))
        site_likelihoods.append(likelihood_for_site)

    tracker = 0
    for column_set in subset_columns:
        new_subset = subset.Subset(a_subset.cfg, column_set)
        list_of_subsets.append(new_subset)
        new_subset.site_lnls_GTRG = site_likelihoods[tracker]
        tracker += 1

    return list_of_subsets

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

