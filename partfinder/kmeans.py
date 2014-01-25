'''functions for k-means splitting'''

import time
import os

import numpy as np
from sklearn.preprocessing import scale
from sklearn.cluster import KMeans
from collections import defaultdict
from util import PhylogenyProgramError

import logging
log = logging.getLogger("kmeans")
import subset_ops


# You can run kmeans in parallel, specify n_jobs as -1 and it will run
# on all cores available.
def kmeans(likelihood_list, number_of_ks=2, n_jobs=1):
    '''Take as input a list of sites, performs k-means clustering on
    sites and returns k centroids and a dictionary with k's as keys
    and lists of sites belonging to that k as values
    '''
    log.debug("Beginning k-means splitting")
    start = time.clock()
    all_rates_list = []
    for site in likelihood_list:
        lk_list = site
        all_rates_list.append(lk_list)

    # Create and scale an array for input into kmeans function
    array = np.array(all_rates_list)
    # array = scale(array)

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


def kmeans_split_subset(cfg, alignment, a_subset, tree_path, number_of_ks = 2):
    """Takes a subset and number of k's and returns
    subsets for however many k's are specified
    """
    a_subset.make_alignment(cfg, alignment)
    phylip_file = a_subset.alignment_path

    # Add option to output likelihoods, *raxml version takes more
    # modfying of the commands in the analyse function
    log.debug("Received subset, now gathering likelihoods")
    processor = cfg.processor

    # For some reason some instances can be analyzed using -f B but not -f g
    # in RAxML, this is to catch those instances and flag the subset as
    # fabricated to add to others later.
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
        else:
            raise PhylogenyProgramError

    # Call processor to parse them likelihoods from the output file. NB these
    # can be site rates as well as likelihoods ToDo: Change the name of this
    # to something other than "likelihood_list" since this isn't really the
    # best description of what it is...
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
    kmeans_results = kmeans(per_site_stat_list,
        number_of_ks)

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


def kmeans_wrapper(cfg, alignment, a_subset, tree_path, max_ks = 10):
    '''This function performs kmeans on
    a specified number of different k's on a dictionary parsed
    from a *_phyml_lk.txt file. It then calculates
    the within sum of squares (wss) of the k's and returns
    the optimal number of k's and a dictionary of
    the k's and sites belonging to the k's as output
    determined by the first time the wss a decrease that is 1/10
    of the initial decrease
    '''
    # Start variable counters
    # Currently we do kmeans with one k for the first cluster, this
    # can be modified to start with two k's. This will results in
    # a slightly larger amount of clusters, we should check the change
    # in AIC scores of both methods and choose one.
    a_subset.make_alignment(cfg, alignment)
    phylip_file = a_subset.alignment_path

    # Add option to output likelihoods, *raxml version takes more
    # modfying of the commands in the analyse function
    processor = cfg.processor

    try:
        processor.gen_per_site_stats(cfg, str(phylip_file),
            str(tree_path))
    except PhylogenyProgramError as e:
        log.error("Total bummer: %s" % e)
        return 1

    per_site_statistics = processor.get_per_site_stats(phylip_file, cfg)
    a_subset.add_per_site_statistics(per_site_statistics)
    likelihood_list = per_site_statistics[0]
    a_subset.site_lnls_GTRG = likelihood_list

    if cfg.kmeans_opt == 1:
        # Set likelihood_list to site likelihoods
        likelihood_list = per_site_statistics[0]
    elif cfg.kmeans_opt == 2:
        # Set likelihood_list to site rates
        likelihood_list = per_site_statistics[2]
    elif cfg.kmeans_opt == 3:
        raise ValueError("Search kmeans_wss cannot be used with kmeans" + 
            " option 3. Use option 1 or 2")
    elif cfg.kmeans_opt == 4:
        raise ValueError("Search kmeans_wss cannot be used with kmeans" + 
            " option 4. Use option 1 or 2")

    count = 1
    new_wss = 0
    list_of_wss = []
    # Calculate a bunch of kmeans on however many max_ks are determined
    for i in range(max_ks):
        # Run kmeans with each count
        site_categories = kmeans(likelihood_list, number_of_ks = count)[1]

        new_likelihood_lists = make_likelihood_list(likelihood_list,
            site_categories)
        # Set previous wss
        previous_wss = new_wss
        # Calculate new wss
        new_wss = within_sum_of_squares(new_likelihood_lists)
        list_of_wss.append(new_wss)
        # Keep the first wss value
        if count == 1:
            first_wss = new_wss
        # Find out the initial decrease
        if count == 2:
            initial_decrease = first_wss - new_wss
            log.info("The initial decrease is %s!" % initial_decrease)
        # Find out if the most recent decrease is less than 10% of the
        # initial decrease
        # print "New wss: " + str(new_wss)
        if count > 2:
            decrease = previous_wss - new_wss
            # print "decrease: " + str(decrease)
            if decrease < ((0.1) * initial_decrease):
                list_of_sites = []
                for k in site_categories:
                    list_of_sites.append(site_categories[k])
                new_subsets = subset_ops.split_subset(a_subset, list_of_sites)
                log.info(
                    "Found that %s categories minimizes within sum of squares"
                     % count)

                return new_subsets
        count += 1
    list_of_sites = []
    for k in site_categories:
        list_of_sites.append(site_categories[k])
    # Make new subsets
    new_subsets = subset_ops.split_subset(a_subset, list_of_sites)
    return new_subsets


def sum_of_squares(list_of_likelihoods):
    '''Input a list of likelihoods, and returns the sum of squares
    of the list
    '''
    sums_of_squares = 0
    mean_likelihood = sum(list_of_likelihoods)/len(list_of_likelihoods)
    for i in list_of_likelihoods:
        sums_of_squares += (i - mean_likelihood)**2
    return sums_of_squares


def within_sum_of_squares(likelihood_lists):
    '''Inputs a list that contains lists of likelihoods from different
    clusters, outputs the within sum of squares
    '''
    wss = 0
    # Call sum_of_squares on each list to get global wss
    for i in likelihood_lists:
        wss += sum_of_squares(i)
    return wss


def make_likelihood_list(likelihood_list, site_categories):
    '''Takes a likelihood_list and a dictionary with kmeans clusters
    as keys and a list of sites belonging to that cluster as values
    as input and returns a list of lists of likelihoods belonging to
    each cluster
    '''
    rate_list = []
    for cluster in site_categories:
        one_list = []
        for site in site_categories[cluster]:
            likelihood = (likelihood_list[site - 1])
            one_list.append(likelihood[0])
        rate_list.append(one_list)
    return rate_list

def kmeans_var_ks(cfg, a_subset, number_of_ks, likelihood_list): 
    '''Takes a desired number of k's and returns the newly
    split subsets into however many k's were desired
    '''
    site_categories = kmeans(likelihood_list, number_of_ks = number_of_ks)[1]
    new_likelihood_lists = make_likelihood_list(likelihood_list, site_categories)
    list_of_sites = []
    for k in site_categories:
        list_of_sites.append(site_categories[k])
    new_subsets = subset_ops.split_subset(a_subset, list_of_sites)
    return new_subsets
