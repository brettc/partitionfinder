
import time
import sys
import re
import os
# import util
import subprocess
import shlex
import csv

from math import log as logarithm
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
def kmeans(likelihood_list, number_of_ks = 2, n_jobs = 1):
    '''Take as input a dictionary made up of site numbers as keys
    and lists of rates as values, performs k-means clustering on
    sites and returns k centroids and a dictionary with k's as keys
    and lists of sites belonging to that k as values
    '''
    start = time.clock()
    all_rates_list = []
    for site in likelihood_list:
        lk_list = site
        log_rate_cat_list = [logarithm(lk) for lk in lk_list]
        all_rates_list.append(log_rate_cat_list)

    # Create and scale an array for input into kmeans function
    array = np.array(all_rates_list)
    array = scale(array)

    # Call scikit_learn's k-means, use "k-means++" to find centroids
    # kmeans_out = KMeans(init='k-means++', n_init = 100)
    kmeans_out = KMeans(init='k-means++', n_clusters = number_of_ks,
        n_init = 100)
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
    time_taken = "k-means took " + str(stop - start) + "seconds"
    log.info(time_taken)

    # Return centroids and dictionary with lists of sites for each k
    return centroid_list, dict(cluster_dict)

def kmeans_split_subset(cfg, alignment, a_subset, number_of_ks = 2):
    """Takes a subset and number of k's and returns
    subsets for however many k's are specified"""
    # pass
    a_subset.make_alignment(cfg, alignment)
    phylip_file = a_subset.alignment_path

    # Add option to output likelihoods, *raxml version takes more 
    # modfying of the commands in the analyse function
    processor = cfg.processor
    program_name = processor.program()

    try:
        # TO DO: still need to make this  call suitable to call RAxML as well,
        # use os.path.join for the items with a slash to that it works on windows
        processor.get_likelihoods("GTRGAMMA", str(phylip_file), 
            "./analysis/start_tree/topology_tree.phy")
    except PhylogenyProgramError as e:
        log.info("Total bummer: %s" % e)
        return 1

    likelihood_list = get_likelihood_list(cfg, phylip_file)

    split_categories = kmeans(likelihood_list, 
        number_of_ks)[1]
    
    list_of_sites = []
    for k in split_categories:
        list_of_sites.append(split_categories[k])

    # This is a quick fix for small clusters, probably can be more
    # sophisticated, it is for testing whether watching for small
    # clusters makes much of a difference during testing
    if number_of_ks == 2:
        for i in split_categories:
            print len(split_categories[i])
            if len(split_categories[i]) < 2:
                return 1

    # Make the new subsets
    new_subsets = subset_ops.split_subset(a_subset, list_of_sites)
    return new_subsets

def kmeans_wrapper(cfg, alignment, a_subset, max_ks = 10):
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
    print processor

    try:
        # TO DO: still need to make this  call suitable to call RAxML as well
        processor.get_likelihoods("GTR", str(phylip_file), 
            "./analysis/start_tree/filtered_source.phy_phyml_tree.txt")
    except Exception as e:
        log.info("Total bummer: %s" % e)
        return 1

    # os.path.join does nothing below. You should use it above. There
    # shouldn't be ANY forward slashes in the code (this will NOT work
    # on windows)
    phyml_lk_file = str(phylip_file) + "_phyml_lk_GTR.txt"

    # Open the phyml output and parse for input into the kmeans
    # function
    likelihood_list = phyml_likelihood_parser(
        phyml_lk_file)
    count = 1
    sum_wss = 0
    new_wss = 0
    list_of_wss = []
    # Calculate a bunch of kmeans on however many max_ks are determined
    for i in range(max_ks):
        # Run kmeans with each count
        site_categories = kmeans(likelihood_list, number_of_ks = count)[1]
        new_likelihood_lists = make_likelihood_list(likelihood_list, site_categories)
        # Set previous wss
        previous_wss = new_wss
        # Calculate new wss
        new_wss = wss(new_likelihood_lists)
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

def ss(list_of_likelihoods):
    '''Input a list of likelihoods, and returns the sum of squares
    of the list
    '''
    sums_of_squares = 0
    # Need to log transform the data as you do during kmeans
    log_list_of_likelihoods = []
    for i in list_of_likelihoods:
        log_list_of_likelihoods.append(logarithm(float(i)))
    mean_likelihood = sum(log_list_of_likelihoods)/len(log_list_of_likelihoods)
    for i in log_list_of_likelihoods:
        sums_of_squares += (i - mean_likelihood)**2
    return sums_of_squares

def wss(likelihood_lists):
    '''Inputs a list that contains lists of likelihoods from different
    clusters, outputs the within sum of squares
    '''
    within_sum_of_squares = 0
    # Call ss on each list to get global wss
    for i in likelihood_lists:
        within_sum_of_squares += ss(i)
    return within_sum_of_squares

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

def get_likelihood_list(cfg, phylip_file):
    '''Runs the appropriate processor to generate the site likelihood
    file, then parses the site likelihoods and returns them as a list
    '''
    phylip_file_split = os.path.split(phylip_file)
    processor = cfg.processor
    program_name = processor.program()
    # Figure out which program to use to calculate site likelihoods
    if program_name == 'phyml':
        phyml_lk_file = ("%s_phyml_lk_GTRGAMMA.txt" % phylip_file)
        # Open the phyml output and parse for input into the kmeans
        # function
        likelihood_list = processor.likelihood_parser(
            phyml_lk_file)[0]
    elif program_name == 'raxml':
        # Once the os.path.split is fixed, will need to change some of the
        # indexes into those lists
        subset_code = phylip_file_split[1].split(".")[0]
        raxml_lnl_file = os.path.join(phylip_file_split[0], 
            ("RAxML_perSiteLLs.%s_GTRGAMMA.txt" % subset_code))
        likelihood_list = processor.likelihood_parser(
            raxml_lnl_file)
    return likelihood_list



if __name__ == "__main__":
    # phylip_filename = sys.argv[1]
    # start = time.clock()
    # run_phyml("-i " + str(phylip_filename) + " -m GTR -o r --print_site_lnl -c 4")
    # stop = time.clock()
    # print "PhyML took " + str(stop-start) + " seconds!"
    # phyml_lk_file = str(phylip_filename) + "_phyml_lk.txt"
    # # rate_cat_likelhood_parser(phyml_lk_file)
    # phyml_likelihood_parser(phyml_lk_file)
    # likelihood_dictionary = phyml_likelihood_parser(phyml_lk_file)
    # print kmeans(likelihood_dictionary[1], number_of_ks = 4)

    # phylip_filename = sys.argv[1]
    # outfile_command = "testing"
    # run_raxml("-s " + str(phylip_filename) + " -m GTRGAMMA -n " + outfile_command + " -y -p 23456")
    # outfile_command2 = "testing2"
    # start = time.clock()
    # run_raxml("-s " + str(phylip_filename) + " -m GTRGAMMA -n " + outfile_command2 + " -f g -p 23456 -z RAxML_parsimonyTree." + outfile_command)
    # stop = time.clock()
    # print "RAxML took " + str(stop-start) + " seconds!"
    # likelihood_dict = raxml_likelihood_parser("RAxML_perSiteLLs." + outfile_command2)
    # print kmeans(likelihood_dict, number_of_ks = 5)

    phhml_lk_file = sys.argv[1]
    likelihood_list = phyml_likelihood_parser(phyml_lk_file)
    kmeans_wrapper(likelihood_list)



