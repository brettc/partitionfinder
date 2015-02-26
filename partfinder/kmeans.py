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
#
import logtools
log = logtools.get_logger()

import time

import numpy as np
from sklearn.cluster import KMeans
from collections import defaultdict
from util import PhylogenyProgramError
from alignment import SubsetAlignment
import util
import os
import subprocess
import shlex
from config import the_config

import subset_ops

# You can run kmeans in parallel, specify n_jobs as -1 and it will run
# on all cores available.
def kmeans(likelihood_list, number_of_ks, n_jobs):
    '''Take as input a list of sites, performs k-means clustering on
    sites and returns k centroids and a dictionary with k's as keys
    and lists of sites belonging to that k as values
    '''
    log.debug("Beginning k-means splitting")
    start = time.clock()

    # Create and scale an array for input into kmeans function
    array = np.array(likelihood_list)

    # Call scikit_learn's k-means, use "k-means++" to find centroids
    # kmeans_out = KMeans(init='k-means++', n_init = 100)
    kmeans_out = KMeans(init='k-means++', n_clusters=number_of_ks,
            n_init=100, n_jobs=n_jobs)
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

def entropy_calc(p):
    # Modify p to include only those elements that are not equal to 0
    p=p[p!=0]
    # The function returns the entropy result
    return np.dot(-p,np.log2(p))

def sitewise_entropies(alignment):
    if the_config.datatype == 'DNA':
        log.debug("Calculating DNA entropies")
        dna_states = "ACGT"
        dna_list = [np.sum(alignment.data == ord(nuc), axis = 0) for nuc in list(dna_states)]
        states = np.array(dna_list, dtype=float)

    elif the_config.datatype == 'protein':
        log.debug("Calculating protein entropies")
        aa_states = "ARNDCQEGHILKMFPSTWYV"
        amino_list = [np.sum(alignment.data == ord(aa), axis = 0) for aa in list(aa_states)]
        states = np.array(amino_list, dtype = float)

    states = states.T
    totals = np.sum(states, axis=1)
    totals.shape = len(states),1

    # for a column of all gaps, we'll have a zero total, so we just hack that here
    totals = np.where(totals==0, 1, totals)

    prob = states/totals

    column_entropy = [[entropy_calc(t)] for t in prob]

    return column_entropy

def rate_parser(rates_name):
    rates_list = []
    the_rates = open(rates_name)
    for rate in the_rates.readlines():
        rates_list.append([float(rate)])
    return rates_list

def run_rates(command, report_errors=True):
    program_name = "fast_TIGER"
    program_path = util.program_path
    program_path = os.path.join(program_path, program_name)
    # command = "\"%s\" %s" % (program_path, command)

    command = "\"%s\" %s" % (program_path, command)

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
            log.error("fast_TIGER did not execute successfully")
            log.error("fast_TIGER output follows, in case it's helpful for finding the problem")
            log.error("%s", stdout)
            log.error("%s", stderr)
        raise PhylogenyProgramError(stdout, stderr)

def sitewise_tiger_rates(cfg, phylip_file):
    if cfg.datatype == 'DNA':
        command = " dna " + phylip_file
    elif cfg.datatype == 'morphology':
        command = " morphology " + phylip_file
    run_rates(command, report_errors=False)
    rates_name = ("%s_r8s.txt" % phylip_file)
    return rate_parser(rates_name)

def get_per_site_stats(alignment, cfg, a_subset):
    if cfg.kmeans == 'entropy':
        sub_align = SubsetAlignment(alignment, a_subset)
        return sitewise_entropies(sub_align)
    elif cfg.kmeans == 'tiger':
        a_subset.make_alignment(cfg, alignment)
        phylip_file = a_subset.alignment_path
        return sitewise_tiger_rates(cfg, str(phylip_file))

def kmeans_split_subset(cfg, alignment, a_subset, tree_path,
                        n_jobs, number_of_ks=2):
    """Takes a subset and number of k's and returns
    subsets for however many k's are specified
    """
    # Get either entropies or TIGER rates
    per_site_stat_list = get_per_site_stats(alignment, cfg, a_subset)

    # Now store all of the per_site_stats with the subset
    a_subset.add_per_site_statistics(per_site_stat_list)

    log.debug("Site info list for subset %s is %s" % (a_subset.name, per_site_stat_list))

    # Perform kmeans clustering on the per site stats
    kmeans_results = kmeans(per_site_stat_list, number_of_ks, n_jobs)
    centroids = kmeans_results[0]
    split_categories = kmeans_results[1]

    list_of_sites = []
    for k in range(len(split_categories)):
        list_of_sites.append(split_categories[k])


    log.debug("Creating new subsets from k-means split")
    # Make the new subsets
    new_subsets = subset_ops.split_subset(a_subset, list_of_sites)

    # Now add the site_lnl centroid to each new subset
    marker = 0
    for s in new_subsets:
        s.centroid = centroids[marker]
        marker += 1

    return new_subsets


