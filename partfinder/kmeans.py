
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

import logging
log = logging.getLogger("kmeans")



def rate_cat_likelhood_parser(phyml_lk_file):
    '''
    BRETT SAYS: Format Code so that it flows 80 columns. (your editor will do
       this for you)

    Note: A new parser has been written that uses csv module and also outputs a
    dictionary of site likelihoods.

    This function takes the path of the phyml_lk file and returns a dictionary
    with sites as keys and a list of likelihoods under different rate
    categories as as the values
    '''

    # BRETT SAYS: Try to use white space to separate logical pieces in a
    # function -- and associate comments

    # Open the file and assign it to a variable
    try:
        phyml_lk_file = open(phyml_lk_file, "r")
    except IOError:
        raise IOError("Could not find site likelihood file!")

    # BRETT SAYS: If you use good variable names, comments are less
    # necessary:
    # count -> line_number
    # site -> site_number
    #
    # Start the count of the line number
    line_number = 0
    # Begin site number count
    site_number = 1

    # BRETT SAYS: We know you're creating a dictionary! Tell us what will be in
    # it...
    # Create dictionary for output
    all_rates_dict = {}

    # Read the file line by line
    for line in phyml_lk_file.readlines():
        # BRETT SAYS: If you need a counter (like site), try this:
        # for line_number, line in enumerate(file.readlines()):

        line_number += 1
        if line_number > 7 and line_number < 1000007:
            # Split the data and objects into a list and remove
            # unnecessary whitespace
            rate_list = " ".join(line.split()).split(" ")

            # This could also be accomplished with regex
            # rate_list = re.sub(' +', ' ', line).split(" ")
            # Raise error if rate variation is not modeled
            if len(rate_list) < 4:
                raise IOError("Rate variation was not modeled in "
                    "PhyML analysis. Please choose > 1 rate category")

            # BRETT SAYS: Why are you doing this? You can just index into the
            # rate_list, why modify it?
            # rate_list[-1] == last element
            # rate_list[-2] == second to last element...
            #
            # Remove first two and last numbers from list
            rate_list.pop(0)
            rate_list.pop(0)
            rate_list.pop(-1)

            # Convert all likelihoods to floats
            # BRETT SAYS: Maybe do this: rate_list = [float(x) for x in rate_list]
            for i in range(len(rate_list)):
                rate_list[i] = float(rate_list[i])

            # Add site number and rate list to dictionary
            all_rates_dict[site] = rate_list
            site_number += 1


        # If there are more than 1,000,000 nts in a phyml_lk.txt
        # file, the numbers merge with the site likelihoods
        elif line_number >= 1000007:
            # Split the data
            rate_list = " ".join(line.split()).split(" ")
            # Remove last and first number from list
            # BRETT SAYS: if you need to modify lists, then you should learn
            # how to use *slicing*. Google it. You can do this x = x[1:-1]
            rate_list.pop(0)
            rate_list.pop(-1)
            # Convert all likelihoods to floats
            for i in range(len(rate_list)):
                rate_list[i] = float(rate_list[i])

            # Add site number and rate list to dictionary
            # BRETT SAYS: Why are you using a dict here? They are consecutive
            # numbers -- just append to a list instead. Get used to ZERO
            # indexed lists too.
            all_rates_dict[site_number] = rate_list
            site_number += 1
    # print all_rates_dict
    #
    phyml_lk_file.close()
    return all_rates_dict

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
        headers.pop(0)
        # Make a dictionary of sites and likelihoods
        like_dict = {}
        for i in list_of_dicts:
            like_dict[int(i[headers[-1]])] = [float(i[headers[0]])]
        phyml_lk_file.close()
        # Return a dictionary of sites and their likelihoods
        return like_dict
    else:
        # Get rid of the two values that csv imports that we don't need
        # along with the posterior means we won't use
        headers.pop(0)
        headers.pop(-1)
        headers.pop(-2)
        # Make a dictionary of sites and likelihoods
        like_dict = {}
        for i in list_of_dicts:
            like_dict[int(i[headers[-1]])] = [float(i[headers[0]])]
        # Now get rid of the site likelihoods from the header list
        headers.pop(0)
        # Figure out how many rate categories there are
        num_rate_cat = len(headers) - 1
        # Now make a dictionary with sites as the keys and a list
        # of likelihoods under the rate categories as values
        like_cat_dict = {}
        for i in list_of_dicts:
            # Make a list with the rate category likelihoods for each site
            like_list = []
            for p in range(num_rate_cat):
                like_list.append(float(i[headers[p]]))
            # Now add that list with it's corresponding site to a dictionary
            # BRETT SAYS: This looks complex.
            like_cat_dict[int(i[headers[-1]])] = like_list
        phyml_lk_file.close()
        # Return a tuple with the site likelihoods and the likelihoods under
        # different rate categories
        return like_dict, like_cat_dict

def raxml_likelihood_parser(raxml_lnl_file):
    '''
    This functiont takes as input the RAxML_perSiteLLs* file from a RAxML -f g
    run, and returns a dictionary of sites and likelihoods to feed into kmeans.

    Note: the likelihoods are already logged, so either we should change the
    kmeans function and tell the PhyML parser to return log likelihoods or we
    should convert these log likelihoods back to regular likelihood scores
    '''
    # See if you can locate the file, then parse the second line
    # that contains the log likelihoods. If it isn't found
    # raise an error
    try:
        # BRETT SAYS: SUPER ugly inconsistent variable name!
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
    num_sites = len(site_lnl_list)

    # Now make a dictionary with sites as values and likelihoods (in a list
    # for import into kmeans()) as the values
    site_num = 1
    site_lnl_dict = {}
    # Add the site and transform the log likelihoods back to likelihoods
    # so it is adequate to feed into kmeans()
    for i in site_lnl_list:
        site_lnl_dict[site_num] = [10**float(i)]
        site_num += 1
    perSiteLL_file.close()
    return site_lnl_dict


# You can run kmeans in parallel, specify n_jobs as -1 and it will run
# on all cores available.
def kmeans(dictionary, number_of_ks = 2, n_jobs = 1):
    '''Take as input a dictionary made up of site numbers as keys
    and lists of rates as values, performs k-means clustering on
    sites and returns k centroids and a dictionary with k's as keys
    and lists of sites belonging to that k as values
    '''
    # Take start time
    start = time.clock()
    all_rates_list = []
    for i in range(1, (len(dictionary) + 1)):
        rate_cat_list = dictionary[i]
        log_rate_cat_list = []
        # Natural log transform the data
        for rate in rate_cat_list:
            log_rate_cat_list.append(logarithm(rate))
        all_rates_list.append(log_rate_cat_list)
    # Create an array with the appropriate number of dimensions
    array = np.array(all_rates_list)
    # Use sklearns preprocessing to scale array
    array = scale(array)
    # Call scikit_learn's k-means, use "k-means++" to find centroids
    # kmeans_out = KMeans(init='k-means++', n_init = 100)
    kmeans_out = KMeans(init='k-means++', n_clusters = number_of_ks,
        n_init = 100)
    # Perform k-means clustering on the array of site likelihoods
    kmeans_out.fit(array)
    # Messing around with different outputs from KMeans object
    # params = kmeans_out.score(array)
    # print params
    # Retrieve centroids
    centroids = kmeans_out.cluster_centers_
    # Append all centroids to a list to return
    centroid_list = []
    for i in centroids:
        centroid_list.append(list(i))
    # Retrieve the cluster number for each site
    rate_categories = kmeans_out.labels_
    # Create dictionary of sites and the cluster they've been assigned to
    site_clus_dict = {}
    char_number = 1
    # Loop through list of rate categories and add them to dictionary
    # with associated character
    for i in rate_categories:
        site_clus_dict[char_number] = i
        char_number += 1
    # Create "defaultdict" to transpose dictionary
    cluster_dict = defaultdict(list)
    # Loop through site_clus_dict and create new key for each cluster
    # with values of lists of sites belonging to that cluster
    for k, v in site_clus_dict.iteritems():
        cluster_dict[v].append(k)
    # Convert defaultdict to regular dictionary
    cluster_dict = dict(cluster_dict)
    # Stop clock and output total time
    stop = time.clock()
    time_taken = "k-means took " + str(stop - start) + "seconds!"
    log.info(time_taken)
    # Return centroids and dictionary with lists of sites for each k
    return centroid_list, cluster_dict

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
    phylip_filename = sys.argv[1]
    outfile_command = "testing"
    run_raxml("-s " + str(phylip_filename) + " -m GTRGAMMA -n " + outfile_command + " -y -p 23456")
    outfile_command2 = "testing2"
    start = time.clock()
    run_raxml("-s " + str(phylip_filename) + " -m GTRGAMMA -n " + outfile_command2 + " -f g -p 23456 -z RAxML_parsimonyTree." + outfile_command)
    stop = time.clock()
    print "RAxML took " + str(stop-start) + " seconds!"
    likelihood_dict = raxml_likelihood_parser("RAxML_perSiteLLs." + outfile_command2)
    print kmeans(likelihood_dict, number_of_ks = 5)

