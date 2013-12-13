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

import subset
import subset_ops
import scheme
import numpy as np
import scipy.spatial.distance

import logging
log = logging.getLogger("cluster")


def get_ranked_list(distance_matrix, subsets, N):
    """
    Return the N closest pairs of subsets in 'subsets' 
    """

    # TODO. If/when anaconda moves to numpy v1.8, we can use
    # a ~7x faster solution here, which uses np.argpartition()
    # see here: http://stackoverflow.com/questions/20540889/
    # extract-the-n-closest-pairs-from-a-numpy-distance-array
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
        norm = float(weights["rate"])/float(np.amax(r_dists))
        r_dists = np.multiply(r_dists, norm)
        distance_arrays.append(r_dists)
    if weights["freqs"] > 0:
        f_dists = scipy.spatial.distance.pdist(np.array(freqs), 'cityblock')
        norm = float(weights["freqs"])/float(np.amax(f_dists))
        f_dists = np.multiply(f_dists, norm)
        distance_arrays.append(f_dists)
    if weights["model"] > 0:
        m_dists = scipy.spatial.distance.pdist(np.array(model), 'cityblock')
        norm = float(weights["model"])/float(np.amax(m_dists))
        m_dists = np.multiply(m_dists, norm)
        distance_arrays.append(m_dists)
    if weights["alpha"] > 0:
        a_dists = scipy.spatial.distance.pdist(np.array(alpha), 'cityblock')
        norm = float(weights["alpha"])/float(np.amax(a_dists))
        a_dists = np.multiply(a_dists, norm)
        distance_arrays.append(a_dists)

    try: 
        final_dists = distance_arrays[0]
    except:
        # no distance arrays, something odd happened
        log.error("Problem calculating subset similarity. Please check that you included"
                " at least one non-zero weight in your clustering weights")
        raise AnalysisError

    for a in distance_arrays[1:]:
        try:
            final_dists = np.add(final_dists, a)
        except:
            log.error("Distance matrices from different parameters are not the same size")
            raise AnalysisError

    return final_dists

def get_distance_matrix(start_scheme, weights):

    #1. get the parameter lists for each subset
    subsets = []  # a list of subsets
    rates = []  # tree length
    freqs = []  # amino acid or base frequencies
    model = []  # model parameters e.g. A<->C
    alpha = []  #alpha parameter of the gamma distribution of rates across sites

    for s in start_scheme.subsets:
        param_dict = s.get_param_values()
        subsets.append(s)
        rates.append([param_dict["rate"]])
        freqs.append(param_dict["freqs"])
        model.append(param_dict["model"])
        alpha.append([param_dict["alpha"]])

    # get pairwise manhattan distances between subsets (a square numpy array)
    final_dists = get_manhattan_matrix(rates, freqs, model, alpha, weights)

    return final_dists, subsets

def get_N_closest_subsets(start_scheme, cfg, N):
    """Find the N most similar groups of subsets in a scheme
    """
    distance_matrix, subsets = get_distance_matrix(start_scheme, cfg.cluster_weights)
    ranked_subset_groupings = get_ranked_list(distance_matrix, subsets, N)
    return ranked_subset_groupings


def make_clustered_scheme(start_scheme, scheme_name, subsets_to_cluster, merged_sub, cfg):

    # 1. Then we define a new scheme with those merged subsets
    new_subsets = start_scheme.subsets - set(subsets_to_cluster)
    new_subsets.add(merged_sub)

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
    closest_subsets = get_N_closest_subsets(start_scheme, cfg, 1)[0]

    merged_sub = subset_ops.merge_subsets(closest_subsets) 

    scheme = make_clustered_scheme(
        start_scheme, scheme_name, closest_subsets, merged_sub, cfg)

    return scheme
