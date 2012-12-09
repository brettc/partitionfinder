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

from cluster import minkowski_distance, genmatrix, printmatrix
import subset
import scheme

import logging
log = logging.getLogger("cluster")


def sum_matrices(mat1, mat2, mat3):
    """sum three matrices of the same size"""
    result = []
    rowindex=0
    for row in mat1:
        result.append([sum(x) for x in zip(mat1[rowindex], mat2[rowindex], mat3[rowindex])] )
        rowindex += 1
    return result

def normalise_and_weight(matrix, normaliser, weight):
    """normalise and weight a 2D matrix (list of lists) using the passed value"""
    
    if normaliser==0:
        #then our matrix must have had a maximum distance of zero anyway, indicating that 
        #everything is already identical (e.g. same amino acid model)
        #in this case we don't need to do any work, because these distances won't count
        #when we sum the matrices anyway, so just return without worrying
        return matrix
        
    rowindex=0
    for row in matrix:
        colindex = 0 # keep track of where we are in the matrix
        for cell in row:
            matrix[rowindex][colindex] = cell*float(weight)/float(normaliser)
            colindex += 1
        rowindex += 1

    return matrix
    
    
def get_minmax(matrix):
    """return the min and max distances in a 2D list
    """

    maxdistance=None    
    mindistance=None
    rowindex=0
    for row in matrix:
        colindex = 0 # keep track of where we are in the matrix
        for cell in row:
            #ignore anything on or below the diagonal
            if (colindex > rowindex):
                if (cell<mindistance or mindistance is None): mindistance=cell
                if (cell>maxdistance or maxdistance is None): maxdistance=cell
            colindex += 1
        rowindex += 1

    return mindistance, maxdistance

def get_ranked_list(matrix, subsets):
    """return the closest subsets defined by a distance matrix 
    usually there will just be a pair that's closer than all other pairs
    BUT, it's feasible (if unlikely) that >2 subsets are equally close.
    This is possible if, e.g. all weights are zero. Then we just want to group all the 
    equally close subsets...
    
    So, we return a list of all the closest subsets
    """
    
    #let's make a dict keyed by the distance in the matrix, using setdefault to 
    #add things, in case there are subsets with identical distances
    distances ={}
    rowindex=0
    for row in matrix:
        colindex=0
        for cell in row:
            #ignore the diagonal and sub-diagonal (it's a symmetrical matrix)
            if (colindex > rowindex):
                #get any subs that we already know are that distance apart as a set
                #default to empty set if it's a new distance                
                subs = distances.setdefault(cell, set()) 
                #add subs that correspond to this cell
                subs.add(subsets[rowindex]) 
                subs.add(subsets[colindex])
            colindex += 1
        rowindex += 1
    
    ordered_subsets = []    

    unique_distances = distances.keys()
    unique_distances.sort()
    
    for d in unique_distances:
        ordered_subsets.append(list(distances[d]))

    return ordered_subsets



def get_closest(matrix, subsets):
    """return the closest subsets defined by a distance matrix 
    usually there will just be a pair that's closer than all other pairs
    BUT, it's feasible (if unlikely) that >2 subsets are equally close.
    This is possible if, e.g. all weights are zero. Then we just want to group all the 
    equally close subsets...
    
    So, we return a list of all the closest subsets
    """
    from random import randrange
    
    min_rows = set()
    min_cols = []
    mindistance=None
    rowindex=0
    for row in matrix:
        colindex=0
        for cell in row:
            #ignore the diagonal
            if (colindex != rowindex):
                if (cell<mindistance or mindistance is None): 
                    mindistance=cell
                    min_rows = set([rowindex])
                elif cell == mindistance:
                    min_rows.add(rowindex)
            colindex += 1
        rowindex += 1

    if len(min_rows)>2:
        log.warning("%d subsets are equally closely related. "
            "This is possible but unlikely if you have nonzero weights in the clustering"
            " algorithm. Double check to make sure..." %len(min_rows))

    subs = []
    subs_names = []
    for index in min_rows:
        subs.append(subsets[index])
        subs_names.append(subsets[index].full_name)

    log.info("Joining subsets: %s" %', '.join(subs_names))

    return subs

def get_distance_matrix(scheme, weights):
    #1. get the parameter lists for each subset
    subsets = [] #a list of subset names, so we know the order things appear in the list
    rates = [] #tree length
    freqs = [] #amino acid or base frequencies
    model = [] #model parameters e.g. A<->C
    
    for s in scheme.subsets:
        param_dict = s.get_param_values()
        subsets.append(s)
        rates.append([param_dict["rate"]])
        freqs.append(param_dict["freqs"])
        model.append(param_dict["model"])

    #2. for each parameter we get a matrix of euclidean distances
    rates_matrix = genmatrix(list=rates, combinfunc=minkowski_distance, symmetric=True, diagonal=0)
    freqs_matrix = genmatrix(list=freqs, combinfunc=minkowski_distance, symmetric=True, diagonal=0)
    model_matrix = genmatrix(list=model, combinfunc=minkowski_distance, symmetric=True, diagonal=0)


    #3. Normalise and weight those euclidean distances,
    min, max = get_minmax(rates_matrix)
    rates_matrix = normalise_and_weight(rates_matrix, max, weights["rate"])
    min, max = get_minmax(freqs_matrix)
    freqs_matrix = normalise_and_weight(freqs_matrix, max, weights["freqs"])
    min, max = get_minmax(model_matrix)
    model_matrix = normalise_and_weight(model_matrix, max, weights["model"])
    
    #4. sum the matrices
    distance_matrix = sum_matrices(rates_matrix, freqs_matrix, model_matrix)
    #printmatrix(distance_matrix)
    
    return distance_matrix

def get_closest_subsets(scheme, weights):
    """Find the closest subsets in a scheme
    """
    distance_matrix = get_distance_matrix(scheme, weights)

    subsets = [] #a list of subset names, so we know the order things appear in the list
    for s in scheme.subsets:
        subsets.append(s)
    
    #5. get the closest pair
    closest_subsets = get_closest(distance_matrix, subsets)
    #print closest_subsets
    
    return closest_subsets

def get_ranked_clustered_schemes(
    start_scheme, name_prefix, cfg):
    """The idea here is to take a scheme, and perform some analyses to find out how the 
    subsets in that scheme cluster.
    
    We then just return the list of schemes, ordered by closest to most distant in the 
    clustering space
    """
    subsets = [] #a list of subset names, so we know the order things appear in the list
    for s in start_scheme.subsets:
        subsets.append(s)
    
    distance_matrix = get_distance_matrix(start_scheme, cfg.cluster_weights)
    
    ranked_subset_groupings = get_ranked_list(distance_matrix, subsets)
        
    ranked_clustered_schemes = []
    counter = 1
    for g in ranked_subset_groupings:
        scheme_name = "%s_%d" %(name_prefix, counter)
        counter +=1
        scheme_g = make_clustered_scheme(start_scheme, scheme_name, g, cfg)
        ranked_clustered_schemes.append(scheme_g)
    
    return ranked_clustered_schemes
    
def make_clustered_scheme(start_scheme, scheme_name, subsets_to_cluster, cfg):
    
    #1. Create a new subset that merges the subsets_to_cluster
    newsub_parts = []
    for s in subsets_to_cluster:
        newsub_parts = newsub_parts + list(s.partitions)
    newsub = subset.Subset(*tuple(newsub_parts))
    
    #2. Then we define a new scheme with those merged subsets
    all_subs = [s for s in start_scheme.subsets]

    #pop out the subsets we're going to join together
    for s in subsets_to_cluster:
        all_subs.remove(s)
    
    #and now we add back in our new subset...
    all_subs.append(newsub)

    #and finally create the clustered scheme
    final_scheme = (scheme.Scheme(cfg, str(scheme_name), *tuple(all_subs)))

    return final_scheme    



def get_nearest_neighbour_scheme(
        start_scheme, scheme_name, cfg):
    """The idea here is to take a scheme, and perform some analyses to find a neighbouring
    scheme, where the neighbour has one less subset than the current scheme. 
    Really this is just progressive clustering, but specified to work well with PartitionFinder
    
    The weights argument allows us to assign different weights to different model parameters
    
    """
            
    #1. First we get the closest subsets, based on some weights. This will almost always
    #   be two subsets, but it's generalised so that it could be all of them...
    #   cluster weights is a dictionary of weights, keyed by: rate, freqs, model
    #   for the overall subset rate, the base/aminoacid frequencies, and the model parameters
    closest_subsets = get_closest_subsets(start_scheme, cfg.cluster_weights)

    scheme = make_clustered_scheme(start_scheme, scheme_name, closest_subsets, cfg)


    return scheme    
        
    