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

import subset_ops
import scheme
import numpy as np
import scipy.spatial.distance
import itertools
from util import PartitionFinderError

import logtools
log = logtools.get_logger()


def get_ranked_list(distance_matrix, subsets, N):
    """
    Return the N closest pairs of subsets in 'subsets' 
    """

    # TODO. If/when anaconda moves to numpy v1.8, we can use
    # a ~7x faster solution here, which uses np.argpartition()
    # see here: http://stackoverflow.com/questions/20540889/
    # extract-the-n-closest-pairs-from-a-numpy-distance-array
    N = int(N)
    closest = distance_matrix.argsort()[:N]
    n = len(subsets)
    ti = np.triu_indices(n, 1)
    r  = zip(ti[0][closest] + 1, ti[1][closest] + 1)

    # and we look up all the subsets that correspond to each distance
    # and add it to our ordered list of subset lists
    ordered_subsets = []
    for i, pair in enumerate(r):
        subset_group = [subsets[i-1] for i in pair]
        ordered_subsets.append(subset_group)

    return ordered_subsets

def get_manhattan_matrix(rates, freqs, model, alpha, weights):

    # get distances between all pairs for all nonzero weights
    # for each matrix we then normalise and weight it
    # by multiplying by the weight/max
    distance_arrays = []

    if weights["rate"] > 0:
        r_dists = scipy.spatial.distance.pdist(np.array(rates), 'cityblock')
        if np.amax(r_dists)>0:        
            norm = float(weights["rate"])/float(np.amax(r_dists))
            r_dists = np.multiply(r_dists, norm)
        distance_arrays.append(r_dists)
    if weights["freqs"] > 0:
        f_dists = scipy.spatial.distance.pdist(np.array(freqs), 'cityblock')
        if np.amax(f_dists)>0:
            norm = float(weights["freqs"])/float(np.amax(f_dists))
            f_dists = np.multiply(f_dists, norm)
        distance_arrays.append(f_dists)
    if weights["model"] > 0:
        m_dists = scipy.spatial.distance.pdist(np.array(model), 'cityblock')
        if np.amax(m_dists)>0:
            norm = float(weights["model"])/float(np.amax(m_dists))
            m_dists = np.multiply(m_dists, norm)
        distance_arrays.append(m_dists)
    if weights["alpha"] > 0:
        a_dists = scipy.spatial.distance.pdist(np.array(alpha), 'cityblock')
        if np.amax(a_dists)>0:
            norm = float(weights["alpha"])/float(np.amax(a_dists))
            a_dists = np.multiply(a_dists, norm)
        distance_arrays.append(a_dists)

    try: 
        final_dists = distance_arrays[0]
    except:
        # no distance arrays, something odd happened
        log.error("Problem calculating subset similarity. Please check that you included"
                " at least one non-zero weight in your clustering weights")
        raise PartitionFinderError

    for a in distance_arrays[1:]:
        try:
            final_dists = np.add(final_dists, a)
        except:
            log.error("Distance matrices from different parameters are not the same size")
            raise PartitionFinderError

    return final_dists

def get_distance_matrix(subsets, weights):

    #1. get the parameter lists for each subset
    rates = []  # tree length
    freqs = []  # amino acid or base frequencies
    model = []  # model parameters e.g. A<->C
    alpha = []  #alpha parameter of the gamma distribution of rates across sites

    for s in subsets:
        param_dict = s.get_param_values()
        rates.append([param_dict["rate"]])
        freqs.append(param_dict["freqs"])
        model.append(param_dict["model"])
        alpha.append([param_dict["alpha"]])

    # get pairwise manhattan distances between subsets (a square numpy array)
    final_dists = get_manhattan_matrix(rates, freqs, model, alpha, weights)

    return final_dists

def get_N_closest_subsets(subsets, cfg, N, distance_matrix = np.matrix([])):
    """Find the N most similar groups of subsets in a scheme
    """

    if not distance_matrix.any():
        distance_matrix = get_distance_matrix(subsets, cfg.cluster_weights)
    ranked_subset_groupings = get_ranked_list(distance_matrix, subsets, N)
    return ranked_subset_groupings



def get_closest_subset(sub, subsets, cfg, distance_matrix = np.matrix([])):
    """Find the most similar subsets to the focal subset 'sub'
    """
    if not distance_matrix.any():
        distance_matrix = get_distance_matrix(subsets, cfg.cluster_weights)
        distance_matrix = scipy.spatial.distance.squareform(distance_matrix)
        
    try:
        col = subsets.index(sub)
    except:
        log.error("Couldn't find the subset you were looking for")
        raise PartitionFinderError

    sub_col = distance_matrix[col,]

    sub_closest = np.min(sub_col[np.nonzero(sub_col)])

    closest_index = int(np.where(sub_col == sub_closest)[0])

    closest_subset = subsets[closest_index]

    return([sub, closest_subset])

def make_clustered_scheme(start_scheme, scheme_name, subsets_to_cluster, merged_sub, cfg):

    # 1. Then we define a new scheme with those merged subsets
    new_subsets = start_scheme.subsets - set(subsets_to_cluster)
    new_subsets.add(merged_sub)

    #3. Create the clustered scheme
    final_scheme = scheme.Scheme(cfg, str(scheme_name), new_subsets)

    return final_scheme


def make_split_scheme(start_scheme, scheme_name, subset_to_split, split_subsets, cfg):

    # 1. Then we define a new scheme with those merged subsets
    new_subsets = start_scheme.subsets - {subset_to_split}

    # 2. add all of the split subsets
    for s in split_subsets:
        new_subsets.add(s)

    #3. Create the clustered scheme
    final_scheme = scheme.Scheme(cfg, str(scheme_name), new_subsets)

    return final_scheme


def get_nearest_neighbour_scheme(start_scheme, scheme_name, cfg):
    """
    The idea here is to take a scheme, and perform some analyses to find a
    neighbouring scheme, where the neighbour has one less subset than the
    current scheme.  Really this is just progressive clustering, but specified
    to work well with PartitionFinder
    """

    # we use [0] becuase the function returns a ranked list of lists of length 1
    subsets = [s for s in start_scheme.subsets]
    closest_subsets = get_N_closest_subsets(subsets, cfg, 1)[0]

    merged_sub = subset_ops.merge_subsets(closest_subsets) 

    scheme = make_clustered_scheme(
        start_scheme, scheme_name, closest_subsets, merged_sub, cfg)

    return scheme




def update_c_matrix(c_matrix, sub_tuples, subsets, diffs):
    """
    Update a symmetric matrix of measurements between subsets.
    Each subset_tuple contains a new subset that is just merged
    from a collection of old subsets.
    """
    if len(c_matrix.shape) == 1:
        c_matrix = scipy.spatial.distance.squareform(c_matrix)

    for t, diff in itertools.izip(sub_tuples, diffs):
        old_subs = t[1]
        i = subsets.index(old_subs[0])
        j = subsets.index(old_subs[1])
        c_matrix[i,j] = c_matrix[j,i] = diff

    return c_matrix

def get_best_pair(c_matrix, best_change, subsets):

    log.debug("C matrix: %s", str(c_matrix))

    if len(c_matrix.shape) == 1:
        log.debug("C matrix shape was == 1")
        c_matrix = scipy.spatial.distance.squareform(c_matrix)
        log.debug("C matrix: %s", str(c_matrix))

    l = np.where(c_matrix==best_change)

    log.debug("C matrix dimensions: %s", str(c_matrix.shape))
    log.debug("Number of subsets: %d", len(subsets))
    log.debug("Location of best_change in c_matrix: %s", str(l))

    # this function is only called if the best_change is <0
    # so we can be sure that we are on an off-diagonal
    # still, we check and throw an error if there's an issue
    s1 = l[0][0]
    s2 = l[1][0] 

    if s1==s2:
        log.error("You can't merge a subset with itself, please check.")
        raise PartitionFinderError

    log.debug("Subsets to merge: %s and %s" %(s1, s2))
    sub1 = subsets[s1]
    sub2 = subsets[s2]
    return (sub1, sub2)

def reset_c_matrix(c_matrix, remove_list, add_list, subsets):

    if len(c_matrix.shape) == 1:
        c_matrix = scipy.spatial.distance.squareform(c_matrix)

    indices = []
    removals = []
    for r in remove_list:
        indices.append(subsets.index(r))
        removals = removals + r.names

    c_matrix = np.delete(c_matrix, indices, 1)
    c_matrix = np.delete(c_matrix, indices, 0)

    additions = []
    for a in add_list:
        additions = additions + a.names
        row = np.array((c_matrix.shape[0]) * [np.inf])
        col = c_matrix.shape[0] * [np.inf]
        col.append(0)
        col = np.array(col)[:, None]
        c_matrix = np.vstack((c_matrix, row))
        c_matrix = np.hstack((c_matrix, col))

    # we can only do this if we've removed the same stuff as we added, check
    if not additions.sort() == removals.sort():
        log.error("Removal and addition of subsets don't add up")
        log.error("Removing: %s", str(removals))
        log.error("Adding: %s", str(additions))
        raise PartitionFinderError

    return c_matrix

def reset_subsets(subsets, remove_list, add_list):

    log.debug("Updating subset list")
    log.debug("Original subset list: %s", str([s.name for s in subsets]))
    log.debug("Subsets to remove: %s", str([s.name for s in remove_list]))
    log.debug("Combined subsets to add: %s", str([s.name for s in add_list]))

    for r in remove_list:
        subsets.pop(subsets.index(r))
    for a in add_list:
        subsets.append(a)
    return subsets


def get_pairs_todo(closest_pairs, c_matrix, subsets):
    pairs_todo = []
    if len(c_matrix.shape) == 1:
        c_matrix = scipy.spatial.distance.squareform(c_matrix)

    for p in closest_pairs:
        i = subsets.index(p[0])
        j = subsets.index(p[1])
        if c_matrix[i,j] == np.inf:
            pairs_todo.append(p)

    return pairs_todo


