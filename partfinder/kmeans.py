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

import time

import numpy as np
from sklearn.cluster import KMeans
from sklearn.preprocessing import scale
from collections import defaultdict
from alignment import SubsetAlignment
import util
from config import the_config
import entropy
import sys
from util import PartitionFinderError
import morph_tiger as mt

import subset_ops

# You can run kmeans in parallel, specify n_jobs as -1 and it will run
# on all cores available.
def kmeans(rate_array, number_of_ks, n_jobs):
    '''Take as input a list of sites, performs k-means clustering on
    sites and returns k centroids and a dictionary with k's as keys
    and lists of sites belonging to that k as values
    '''
    log.debug("Beginning k-means splitting")
    start = time.clock()

    # Create and scale an array for input into kmeans function
    array = scale(rate_array)

    # Call scikit_learn's k-means, use "k-means++" to find centroids
    # kmeans_out = KMeans(init='k-means++', n_init = 100)
    kmeans_out = KMeans(init='k-means++', n_clusters=number_of_ks,
            n_init=100, n_jobs=n_jobs, random_state = 2147483647)

    # Perform k-means clustering on the array of site likelihoods
    kmeans_out.fit(array)

    # Retrieve centroids
    centroids = kmeans_out.cluster_centers_
    # Add all centroids to a list to return
    centroid_list = [list(centroid) for centroid in centroids]

    # Retrieve a list with the cluster number for each site
    rate_categories = kmeans_out.labels_
    rate_categories = list(rate_categories)

    # Transpose the list of cluster numbers to a dictionary with
    # cluster numbers as keys and list of sites belonging to that
    # cluster as the value
    cluster_dict = defaultdict(list)
    for num in range(len(rate_categories)):
        cluster_dict[rate_categories[num]].append(num + 1)

    stop = time.clock()
    time_taken = "k-means splitting took %s seconds" % (stop - start)
    log.debug(time_taken)

    # Return centroids and dictionary with lists of sites for each k
    return centroid_list, dict(cluster_dict)

def rate_parser(rates_name):
    rates_list = []
    the_rates = open(rates_name)
    for rate in the_rates.readlines():
        rates_list.append([float(rate)])
    rate_array = np.array(rates_list)
    return rate_array


def get_per_site_stats(alignment, cfg, a_subset):
    if cfg.kmeans == 'entropy':
        sub_align = SubsetAlignment(alignment, a_subset)
        return entropy.sitewise_entropies(sub_align)
    elif cfg.kmeans == 'tiger' and cfg.datatype == 'morphology':
        sub_align = SubsetAlignment(alignment, a_subset)
        set_parts = mt.create_set_parts(sub_align)
        rates = mt.calculate_rates(set_parts)
        return rates
    else: #wtf
        log.error("Unkown option passed to 'kmeans'. Please check and try again")
        raise PartitionFinderError


def kmeans_split_subset(cfg, alignment, a_subset, tree_path,
                        n_jobs, number_of_ks=2):
    """Takes a subset and number of k's and returns
    subsets for however many k's are specified
    """
    # Get either entropies or TIGER rates
    per_site_stat_list = get_per_site_stats(alignment, cfg, a_subset)

    # Now store all of the per_site_stats with the subset
    a_subset.add_per_site_statistics(per_site_stat_list)
    log.debug("The per site statistics for the first 10 sites of subset %s are %s"
        % (a_subset.name, per_site_stat_list[0:10]))

    # Perform kmeans clustering on the per site stats
    kmeans_results = kmeans(per_site_stat_list, number_of_ks, n_jobs)
    centroids = kmeans_results[0]
    split_categories = kmeans_results[1]

    list_of_sites = []
    for k in range(len(split_categories)):
        list_of_sites.append(split_categories[k])

    log.debug("# split categories: %d" % len(split_categories))

    log.debug("Creating new subsets from k-means split")
    # Make the new subsets
    new_subsets = subset_ops.split_subset(a_subset, list_of_sites)

    # Now add the site_lnl centroid to each new subset
    marker = 0
    for s in new_subsets:
        s.centroid = centroids[marker]
        marker += 1

    return new_subsets
