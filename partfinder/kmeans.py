
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

def phyml_likelihood_parser(phyml_lk_file):
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
        # Make a list of site likelihoods
        likelihood_list = [[float(site[headers[1]])] for site in list_of_dicts]
        return likelihood_list

    else:
        # Make a list of site likelihoods
        likelihood_list = [[float(site[headers[1]])] for site in list_of_dicts]

        # Now make a list of lists of site likelihoods under different 
        # rate categories
        lk_rate_list = []
        for i in list_of_dicts:
            ind_lk_list = []
            # Pull the likelihood from each rate category by calling the 
            # appropriate key from "headers"
            for num in range(2, len(headers) - 3):
                ind_lk_list.append(float(i[headers[num]]))
            # Now add the list of likelihoods for the site to a master list
            lk_rate_list.append(ind_lk_list)

        # Return both the list of site likelihoods and the list of lists of
        # likelihoods under different rate categories
        return likelihood_list, lk_rate_list

def raxml_likelihood_parser(raxml_lnl_file):
    '''
    This function takes as input the RAxML_perSiteLLs* file from a RAxML -f g
    run, and returns a dictionary of sites and likelihoods to feed into kmeans.

    Note: the likelihoods are already logged, so either we should change the
    kmeans function and tell the PhyML parser to return log likelihoods or we
    should convert these log likelihoods back to regular likelihood scores
    '''
    # See if you can locate the file, then parse the second line
    # that contains the log likelihoods. If it isn't found
    # raise an error
    try:
        with open(str(raxml_lnl_file)) as raxml_lnl_file:
            line_num = 1
            for line in raxml_lnl_file.readlines():
                if line_num == 2:
                    site_lnl_list = line.split(" ")
                line_num += 1
    except IOError:
        raise IOError("Could not locate per site log likelihood file")

    # Get rid of the new line character and the first "tr1" from
    # the first element in the list
    site_lnl_list[0] = site_lnl_list[0].strip("tr1\t")
    site_lnl_list.pop(-1)

    # You have to take the antilog the output will be likelihoods
    # like the output from PhyML
    site_lk_list = [[10**float(site)] for site in site_lnl_list]

    raxml_lnl_file.close()
    return site_lk_list


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
    phylip_file_split = os.path.split(phylip_file)
    print phylip_file_split


    # Add option to output likelihoods, *raxml version takes more 
    # modfying of the commands in the analyse function
    processor = cfg.processor
    program_name = processor.program()
    print program_name

    try:
        # TO DO: still need to make this  call suitable to call RAxML as well,
        # use os.path.join for the items with a slash to that it works on windows
        processor.get_likelihoods("GTRGAMMA", str(phylip_file), 
            "./analysis/start_tree/topology_tree.phy")
    except PhylogenyProgramError as e:
        log.info("Total bummer: %s" % e)
        return 1

    if program_name == 'phyml':
        phyml_lk_file = (str(phylip_file) + 
        "_phyml_lk_GTR.txt")
        # Open the phyml output and parse for input into the kmeans
        # function
        likelihood_dictionary = phyml_likelihood_parser(
        phyml_lk_file)

    elif program_name == 'raxml':
        print phylip_file
        subset_code = phylip_file_split[1].split(".")[0]
        print subset_code
        raxml_tree_puzzler_file = os.path.join(phylip_file_split[0], 
            ("RAxML_perSiteLLs." + str(subset_code) + "_GTRGAMMA.txt"))
        print raxml_tree_puzzler_file
        likelihood_dictionary = raxml_likelihood_parser(raxml_tree_puzzler_file)

    split_categories = kmeans(likelihood_dictionary, 
        number_of_ks)[1]
    print split_categories
    
    list_of_sites = []
    for k in split_categories:
        list_of_sites.append(split_categories[k])

    # This is a quick fix for small clusters, probably can be more
    # sophisticated
    if number_of_ks == 2:
        for i in split_categories:
            print len(split_categories[i])
            if len(split_categories[i]) < 10:
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



