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


def kmeans_split_subset(cfg, alignment, a_subset, tree_path,
                        n_jobs, number_of_ks=2):
    """Takes a subset and number of k's and returns
    subsets for however many k's are specified
    """
    a_subset.make_alignment(cfg, alignment)
    phylip_file = a_subset.alignment_path

    # Add option to output likelihoods, *raxml version takes more
    # modfying of the commands in the analyse function
    log.debug("Received subset, now gathering likelihoods")
    processor = cfg.processor

    # This is where we catch and deal with errors from PhyML/RAxML
    # we only catch very specific errors.
    try:
        processor.gen_per_site_stats(cfg, str(phylip_file),
            str(tree_path))
    except PhylogenyProgramError as e:
        error1 = "that consist entirely of undetermined values"
        if e.stdout.find(error1) != -1:
            log.warning("The program was unable to calculate site" +
                " likelihoods because of undetermined values for one or" +
                " more taxon, we will move to the next subset")
            a_subset.fabricated = True
            a_subset.analysis_error = "entirely undetermined values"
            return 1
        elif e.stdout.find("base frequency for state number") != 1:
            log.warning("The program was unable to calculate site" +
                " likelihoods because the frequency for one of the" +
                " nucleotides is equal to zero")
            a_subset.fabricated = True
            a_subset.analysis_error = "state frequency equal to zero"
            return 1
        elif e.stdout.find("Round_Optimize") != 1:
            log.warning("The program couldn't analyse this subset" +
                " because the optimisation failed")
            a_subset.fabricated = True
            a_subset.analysis_error = "Optimisation failed"
            return 1

        else:
            raise PhylogenyProgramError

    # Call processor to parse them likelihoods from the output file. NB these
    # can be site rates as well as likelihoods 
    per_site_statistics = processor.get_per_site_stats(phylip_file, cfg)

    # Now store all of the per_site_stats with the subset
    a_subset.add_per_site_statistics(per_site_statistics)

    # Now figure out which list the user wants and use that for the kmeans
    # splitting
    if cfg.kmeans_opt == 1:
        # Set the per_site_stat_list to site likelihoods only
        per_site_stat_list = per_site_statistics[0]
    if cfg.kmeans_opt == 2:
        # Set the per_site_stat_list to site rates only
        per_site_stat_list = per_site_statistics[2]
    if cfg.kmeans_opt == 3:
        # Set the per_site_stat_list to site rates and likelihoods together
        per_site_stat_list = per_site_statistics[3]
    if cfg.kmeans_opt == 4:
        # Set the per_site_stat_list to likelihoods under each gamma rate
        # category
        per_site_stat_list = per_site_statistics[1]

    log.debug("Site info list for subset %s is %s" % (a_subset.name, per_site_stat_list))

    # We use this variable in the case that, as a result of the split with
    # kmeans, the subset becomes unanalysable. In that instance, these will be
    # transferred to the new subset and the sum taken as a proxy for the
    # overal lnl
    a_subset.site_lnls_GTRG = per_site_statistics[0]

    # Perform kmeans clustering on the likelihoods
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


