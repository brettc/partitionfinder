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
import math
import scheme
import submodels
from analysis import Analysis, AnalysisError
from alignment import SubsetAlignment
import neighbour
import kmeans
from subset import Subset
import subset_ops
import entropy
from scipy import spatial
from scipy.misc import comb
import numpy as np
from config import the_config

log = logtools.get_logger()


class UserAnalysis(Analysis):
    def do_analysis(self):
        log.info("Performing User analysis")
        current_schemes = the_config.user_schemes
        scheme_count = len(current_schemes)
        subset_count = len(the_config.user_subsets)

        the_config.progress.begin(scheme_count, subset_count)
        if scheme_count > 0:
            for s in current_schemes:
                res = self.analyse_scheme(s)

                # Write out the scheme
                if not the_config.quick:
                    the_config.reporter.write_scheme_summary(s, res)
        else:
            log.error(
                "Search set to 'user', but no user schemes detected in .cfg file. Please check.")
            raise AnalysisError

        the_config.progress.end()

        the_config.reporter.write_best_scheme(self.results)


class StrictClusteringAnalysis(Analysis):
    """
    This analysis uses model parameters to guess at similar partitions, then
    just joins them together this is much less accurate than other methods, but
    a LOT quicker - it runs in order N time (where N is the number of initial
    datablocks), whereas the greedy algorithm is still N squared.
    """

    def do_analysis(self):
        log.info("Performing strict clustering analysis")

        partnum = len(the_config.user_subsets)
        subset_count = 2 * partnum - 1
        scheme_count = partnum
        the_config.progress.begin(scheme_count, subset_count)

        # Start with the most partitioned scheme
        start_description = range(partnum)
        start_scheme = scheme.create_scheme(
            the_config, "start_scheme", start_description)

        # Analyse our first scheme
        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        self.analyse_scheme(start_scheme)

        # Current scheme number
        cur_s = 2

        # Now we try out all clusterings of the first scheme, to see if we can
        # find a better one
        while True:
            log.info("***Strict clustering algorithm step %d of %d***" %
                     (cur_s - 1, partnum - 1))

            # Calculate the subsets which are most similar
            # e.g. combined rank ordering of euclidean distances
            # Could combine average site-rates, q matrices, and frequencies
            scheme_name = "step_%d" % (cur_s - 1)
            clustered_scheme = neighbour.get_nearest_neighbour_scheme(
                start_scheme, scheme_name, the_config)

            # Now analyse that new scheme
            cur_s += 1
            self.analyse_scheme(clustered_scheme)

            # Stop when we've analysed the scheme with all subsets combined...
            if len(set(clustered_scheme.subsets)) == 1:
                # ... then it's the scheme with everything together
                break
            else:
                # We keep going
                start_scheme = clustered_scheme

        the_config.progress.end()
        the_config.reporter.write_best_scheme(self.results)


class AllAnalysis(Analysis):
    def do_analysis(self):
        log.info("Performing complete analysis")
        partnum = len(the_config.user_subsets)

        scheme_count = submodels.count_all_schemes(partnum)
        subset_count = submodels.count_all_subsets(partnum)
        the_config.progress.begin(scheme_count, subset_count)

        # Iterate over submodels, which we can turn into schemes afterwards in the loop
        model_iterator = submodels.submodel_iterator([], 1, partnum)

        scheme_name = 1
        for m in model_iterator:
            s = scheme.model_to_scheme(m, scheme_name, the_config)
            scheme_name = scheme_name + 1
            res = self.analyse_scheme(s)

            # Write out the scheme
            if not the_config.quick:
                the_config.reporter.write_scheme_summary(s, res)

        the_config.reporter.write_best_scheme(self.results)


class GreedyAnalysis(Analysis):

    @logtools.log_info(log, "Performing Greedy Analysis")
    def do_analysis(self):
        '''A greedy algorithm for heuristic partitioning searches'''

        partnum = len(the_config.user_subsets)
        scheme_count = submodels.count_greedy_schemes(partnum)
        subset_count = submodels.count_greedy_subsets(partnum)

        the_config.progress.begin(scheme_count, subset_count)

        # Start with the most partitioned scheme, and record it.
        with logtools.indented(log, "*** Analysing starting scheme ***"):
            the_config.progress.begin(scheme_count, partnum)
            start_scheme = scheme.create_scheme(
                the_config, "start_scheme", range(partnum))
            start_result = self.analyse_scheme(start_scheme)
            start_score = start_result.score
            if not the_config.quick:
                the_config.reporter.write_scheme_summary(
                    self.results.best_scheme, self.results.best_result)

        subsets = [s for s in start_scheme.subsets]

        step = 1
        while len(set(start_scheme.subsets)) > 1:
            with logtools.indented(log, "***Greedy algorithm step %d***" % step):
                name_prefix = "step_%d" % (step)

                # get distances between subsets
                max_schemes = comb(len(start_scheme.subsets), 2)

                # this is a fake distance matrix, so that the greedy algorithm
                # can use all the tricks of the relaxed clustering algorithm
                dim = len(subsets)
                d_matrix = np.zeros((((dim*dim)-dim))/2)
                d_matrix[:] = np.inf

                if step == 1:
                    # Now initialise a change in info score matrix to inf
                    c_matrix = np.empty(d_matrix.shape)
                    c_matrix[:] = np.inf
                    c_matrix = spatial.distance.squareform(c_matrix)

                # 1. pick top N subset pairs from distance matrix
                cutoff = max_schemes # this defines the greedy algorithm: we look at all schemes

                closest_pairs = neighbour.get_N_closest_subsets(
                    subsets, the_config, cutoff, d_matrix)

                # 2. analyse subsets in top N that have not yet been analysed
                pairs_todo = neighbour.get_pairs_todo(closest_pairs, c_matrix, subsets)
                if len(pairs_todo)>0:
                    log.info("Analysing %d new subset pairs" % len(pairs_todo))
                    new_subs = []
                    sub_tuples = []
                    for pair in pairs_todo:
                        new_sub = subset_ops.merge_subsets(pair)
                        new_subs.append(new_sub)
                        sub_tuples.append((new_sub, pair))

                    the_config.progress.begin(scheme_count, len(new_subs))
                    self.analyse_list_of_subsets(new_subs)

                    # 3. for all K new subsets, update improvement matrix and find best pair
                    log.info("Finding the best partitioning scheme")
                    diffs = []
                    scheme_name = "step_%d" %(step)
                    for t in sub_tuples:
                        pair_merged = t[0]
                        pair = t[1]
                        new_scheme = neighbour.make_clustered_scheme(
                                start_scheme, scheme_name, pair, pair_merged, the_config)
                        r = self.analyse_scheme(new_scheme)
                        diff = r.score - start_score
                        diffs.append(diff)

                    c_matrix = neighbour.update_c_matrix(c_matrix, sub_tuples, subsets, diffs)


                # 4. Find the best pair of subsets, and build a scheme based on that
                # note that this matrix includes diagonals, which will all be zero
                # since this is equivalent to comparing a scheme to itself.
                # so we need to be careful to only proceed if we have a negative change
                # which indicates an improvement in the score
                best_change = np.amin(c_matrix)

                log.debug("Biggest improvement in info score: %s", str(best_change))

                if best_change>=0:
                    log.info("Found no schemes that improve the score, stopping")
                    break

                best_pair = neighbour.get_best_pair(c_matrix, best_change, subsets)

                best_merged = subset_ops.merge_subsets(best_pair)
                best_scheme = neighbour.make_clustered_scheme(
                    start_scheme, scheme_name, best_pair, best_merged, the_config)
                best_result = self.analyse_scheme(best_scheme)

                # the best change can get updated a fraction at this point
                # because calaculting the info score on the whole alignment
                # is a little different from doing it on the one subset
                best_change = self.results.best_score - start_score


                log.info("Best scheme combines subsets: '%s' and '%s'" %(best_pair[0].name, best_pair[1].name))


                log.info("The best scheme improves the %s score by %.2f to %.1f",
                    the_config.model_selection,
                    np.abs(best_change),
                    self.results.best_score)
                start_scheme = best_scheme
                start_score = best_result.score

                log.debug("Best pair: %s", str([s.name for s in best_pair]))
                log.debug("Merged into: %s", str([best_merged.name]))

                # 5. reset_c_matrix and the subset list
                c_matrix = neighbour.reset_c_matrix(c_matrix, list(best_pair), [best_merged], subsets)

                # we updated the subset list in a special way, which matches how we update the c matrix:
                subsets = neighbour.reset_subsets(subsets, list(best_pair), [best_merged])

                if not the_config.quick:
                    the_config.reporter.write_scheme_summary(
                        best_scheme, best_result)

                step += 1

        log.info("Greedy algorithm finished after %d steps" % step)
        log.info("Best scoring scheme is scheme %s, with %s score of %.3f"
                 % (self.results.best_scheme.name, the_config.model_selection,
                    self.results.best_score))

        the_config.reporter.write_best_scheme(self.results)

class RelaxedClusteringAnalysis(Analysis):
    '''
    A fast relaxed clustering algorithm for heuristic partitioning searches

    1. Rank subsets by their similarity (defined by clustering-weights)
    2. Analyse min(cluster-percent or cluster-max) most similar schemes
    3. Sequentially perform all groupings that imporve the AICc/BIC score, in order of improvement
    4. Analyse resulting scheme, iterate to 2.
    5. Quit if no improvements.
    '''

    def clean_scheme(self, start_scheme):
        # Here we look for and fix up subsets that are too small or don't have all states
        keep_going = 1
        merges = 0
        if keep_going == 1:
            with logtools.indented(log, "*** Checking subsets from scheme '%s' meet --min-subset-size and --all_states settings ***" %start_scheme.name):
                while keep_going > 0:

                    subsets = [s for s in start_scheme.subsets]

                    # sort the subsets, to keep results consistent over re-runs
                    subsets.sort(key = lambda x: 1.0/float(len(x.columns)))

                    # run through all subsets
                    for i, sub in enumerate(subsets):
                        found = 0
                        state_problems = self.alignment.check_state_probs(sub, the_config)

                        if  (
                                len(sub.columns) < the_config.min_subset_size or
                                state_problems == True
                            ):

                            # merge that subset with nearest neighbour
                            new_pair = neighbour.get_closest_subset(sub, subsets, the_config)

                            log.info("Subset '%s' will be merged with subset '%s'" %(new_pair[0].name, new_pair[1].name))
                            new_pair_merged = subset_ops.merge_subsets(new_pair)
                            start_scheme = neighbour.make_clustered_scheme(
                                    start_scheme, "cleaned_scheme", new_pair, new_pair_merged, the_config)
                            the_config.progress.begin(1, 1)
                            self.analyse_scheme(start_scheme)
                            subsets = [s for s in start_scheme.subsets]
                            merges = merges + 1
                            found = 1
                            break


                    # if we got to here, there were no subsets to merge
                    if found == 0:
                        keep_going = 0

                    if len(subsets) == 1:
                        log.error("The settings you have used for --all-states and/or --min-subset-size mean that all of your subsets have been merged into one prior to any analysis. Thus, no analysis is necessary. Please check and try again")
                        raise AnalysisError

                log.info("%d subsets merged because of --min-subset-size and/or --all-states settings" % merges)
        return(start_scheme)

    @logtools.log_info(log, "Performing relaxed clustering analysis")
    def do_analysis(self):

        # initialisation steps
        model_selection = the_config.model_selection
        partnum = len(the_config.user_subsets)

        if the_config.cluster_max == -987654321:
            the_config.cluster_max = max([1000, (10 * len(the_config.user_subsets))])
            log.info("Set rcluster-max to %d" %the_config.cluster_max) 

        scheme_count = submodels.count_relaxed_clustering_schemes(
            partnum, the_config.cluster_percent, the_config.cluster_max)
        subset_count = submodels.count_relaxed_clustering_subsets(
            partnum, the_config.cluster_percent, the_config.cluster_max)

        log.info("PartitionFinder will have to analyse %d subsets to"
                 " complete this analyses" % subset_count)
        the_config.progress.begin(scheme_count, subset_count)

        # Start with the most partitioned scheme, and record it.
        with logtools.indented(log, "*** Analysing starting scheme ***"):
            the_config.progress.begin(scheme_count, partnum)
            start_scheme = scheme.create_scheme(
                the_config, "start_scheme", range(partnum))
            start_result = self.analyse_scheme(start_scheme)
            start_score = start_result.score
            if not the_config.quick:
                the_config.reporter.write_scheme_summary(
                    self.results.best_scheme, self.results.best_result)



        subsets = [s for s in start_scheme.subsets]
        partnum = len(subsets)
        step = 1
        while True:
            with logtools.indented(log, "*** Relaxed clustering algorithm step %d of up to %d ***"
                % (step, partnum - 1)):

                # get distances between subsets
                max_schemes = comb(len(start_scheme.subsets), 2)
                log.info("Measuring the similarity of %d subset pairs" % max_schemes)
                d_matrix = neighbour.get_distance_matrix(subsets,
                    the_config.cluster_weights)

                if step == 1:
                    # Now initialise a change in info score matrix to inf
                    c_matrix = np.empty(d_matrix.shape)
                    c_matrix[:] = np.inf
                    c_matrix = spatial.distance.squareform(c_matrix)

                # 1. pick top N subset pairs from distance matrix
                cutoff = int(math.ceil(max_schemes * (the_config.cluster_percent * 0.01)))
                if cutoff <= 0: cutoff = 1
                if the_config.cluster_max != None and cutoff > the_config.cluster_max:
                    cutoff = the_config.cluster_max
                log.info("Choosing the %d most similar subset pairs" % cutoff)
                closest_pairs = neighbour.get_N_closest_subsets(
                    subsets, the_config, cutoff, d_matrix)

                # 2. analyse K subsets in top N that have not yet been analysed
                pairs_todo = neighbour.get_pairs_todo(closest_pairs, c_matrix, subsets)
                if len(pairs_todo)>0:
                    log.info("Analysing %d new subset pairs" % len(pairs_todo))
                    new_subs = []
                    sub_tuples = []
                    for pair in pairs_todo:
                        new_sub = subset_ops.merge_subsets(pair)
                        new_subs.append(new_sub)
                        sub_tuples.append((new_sub, pair))

                    the_config.progress.begin(scheme_count, len(new_subs))
                    self.analyse_list_of_subsets(new_subs)

                    # 3. for all K new subsets, update improvement matrix and find best pair
                    log.info("Finding the best partitioning scheme")
                    diffs = []
                    scheme_name = "step_%d" %(step)
                    for t in sub_tuples:
                        pair_merged = t[0]
                        pair = t[1]
                        new_scheme = neighbour.make_clustered_scheme(
                                start_scheme, scheme_name, pair, pair_merged, the_config)
                        r = self.analyse_scheme(new_scheme)
                        diff = r.score - start_score
                        diffs.append(diff)

                    c_matrix = neighbour.update_c_matrix(c_matrix, sub_tuples, subsets, diffs)

                # 4. Find the best pair of subsets, and build a scheme based on that
                # note that this matrix includes diagonals, which will all be zero
                # since this is equivalent to comparing a scheme to itself.
                # so we need to be careful to only proceed if we have a negative change
                # which indicates an improvement in the score
                best_change = np.amin(c_matrix)
                best_scheme = start_scheme

                if best_change>=0:
                    log.info("Found no schemes that improve the score, stopping")
                    break

                median_improvement = np.median(c_matrix[c_matrix<0])

                while best_change <= median_improvement:

                    best_pair = neighbour.get_best_pair(c_matrix, best_change, subsets)
                    best_merged = subset_ops.merge_subsets(best_pair)
                    best_scheme = neighbour.make_clustered_scheme(
                        start_scheme, scheme_name, best_pair, best_merged, the_config)
                    start_scheme = best_scheme

                    log.info("Combining subsets: '%s' and '%s'" %(best_pair[0].name, best_pair[1].name))
                    log.debug("This improves the %s score by: %s", the_config.model_selection, str(abs(best_change)))

                    # reset_c_matrix and the subset list
                    c_matrix = neighbour.reset_c_matrix(c_matrix, list(best_pair), [best_merged], subsets)

                    # we update the subset list in a way that means its structure tracks the c-matrix
                    subsets = neighbour.reset_subsets(subsets, list(best_pair), [best_merged])

                    best_change = np.amin(c_matrix)

                    if the_config.search == 'rcluster':
                        break
                        # otherwise we are using rclusterf, which continues in this loop
                        # i.e. with rcluster we just take the single best change


                # the best change can get updated a fraction at this point
                # because calaculting the info score on the whole alignment
                # is a little different from doing it on the one subset
                best_result = self.analyse_scheme(best_scheme)
                best_change = self.results.best_score - start_score


                log.info("The best scheme has %d subsets and improves the %s score by %.2f to %.1f",
                    len(best_scheme.subsets),
                    the_config.model_selection,
                    np.abs(best_change),
                    self.results.best_score)
                start_scheme = best_scheme
                start_score = best_result.score


                if not the_config.quick:
                    the_config.reporter.write_scheme_summary(
                        best_scheme, best_result)


                if len(set(start_scheme.subsets)) == 1:
                    break

                step += 1

        log.info("Relaxed clustering algorithm finished after %d steps" % step)
        log.info("Best scoring scheme is scheme %s, with %s score of %.3f"
                 % (self.results.best_scheme.name, model_selection,
                    self.results.best_score))

        if the_config.min_subset_size or the_config.all_states:
            best_scheme = self.clean_scheme(self.results.best_scheme)
            best_result = self.analyse_scheme(best_scheme)

            # scores after cleaning can be worse, so we reset these trackers...
            self.results.best_result = best_result
            self.results.best_score = best_result.score
            self.results.best_scheme = best_scheme
            log.info("Best scoring scheme after cleaning is scheme %s, with %s score of %.3f"
                     % (self.results.best_scheme.name, model_selection,
                        self.results.best_score))


        the_config.reporter.write_best_scheme(self.results)


class KmeansAnalysis(Analysis):


    def split_subsets(self, start_subsets, tree_path):
        split_subs = {}
        for i, sub in enumerate(start_subsets):

            # here we can test if the alignment has all states:
            state_probs = self.alignment.check_state_probs(sub, the_config)

            if  (
                    len(sub.columns) == 1 or
                    len(sub.columns) < the_config.min_subset_size or
                    state_probs == True or
                    sub.dont_split == True
                ):
                split_subs[sub] = [sub]
                log.info("Subset %d: %d sites not splittable" %(i+1, len(sub.columns)))
            else:
                split = kmeans.kmeans_split_subset(
                    the_config, self.alignment, sub, tree_path, n_jobs=self.threads)

                if split == 1:  # we couldn't analyse the big subset
                    sub.dont_split = True # never try to split this subset again
                    split_subs[sub] = [sub]  # so we keep it whole
                    log.info("Subset %d: %d sites not splittable" %(i+1, len(sub.columns)))

                elif len(split) == 1:
                    # in some cases (i.e. all site params are equal) kmeans
                    # cannot split subsets, so we get back the same as we put in
                    sub.dont_split = True # never try to split this subset again
                    split_subs[sub] = [sub]  # so we keep it whole
                    log.info("Subset %d: %d sites not splittable" %(i+1, len(sub.columns)))

                elif min([len(split[0].columns), len(split[1].columns)]) < the_config.min_subset_size:
                    # we don't split it if either of the two daughter subset was
                    # smaller than the minimum allowable size
                    sub.dont_split = True # never try to split this subset again
                    split_subs[sub] = [sub]  # so we keep it whole
                    log.info("Subset %d: %d sites not splittable" %(i+1, len(sub.columns)))

                elif (
                        self.alignment.check_state_probs(split[0], the_config) or
                        self.alignment.check_state_probs(split[1], the_config)
                     ):
                    # we don't split it if either of the daughter subsets has problems with
                    # the state frequencies
                    sub.dont_split = True # never try to split this subset again
                    split_subs[sub] = [sub]  # so we keep it whole
                    log.info("Subset %d: %d sites not splittable" %(i+1, len(sub.columns)))

                else:  # we could analyse the big subset
                    split_subs[sub] = split  # so we split it into >1
                    log.info("Subset %d: %d sites split into %d and %d"
                            %(i+1, len(sub.columns), len(split[0].columns), len(split[1].columns)))

        return split_subs

    def finalise_fabrication(self, start_subsets, step):

        fabricated_subsets = []
        for s in start_subsets:

            # here we put a sensible lower limit on the size of subsets
            if len(s.columns) < the_config.min_subset_size:
                s.fabricated = True
                log.debug("Subset %s with only %d sites found" %(s.subset_id, len(s.columns)))

            # here we can test if the alignment has all states:
            state_probs = self.alignment.check_state_probs(s, the_config)
            if state_probs:
                s.fabricated = True
                log.debug("Subset %s does not have all states in the alignment", s.subset_id)

            if s.fabricated:
                fabricated_subsets.append(s)
                log.debug("added %s to fabricated subset", s.name)

        if fabricated_subsets:
            with logtools.indented(log, "Finalising partitioning scheme"):
                log.debug("There are %d/%d fabricated subsets"
                          % (len(fabricated_subsets), len(start_subsets)))

                i = 1
                while fabricated_subsets:

                    all_subs = start_subsets

                    # occasionally subsets with all value == 0.0 are given a
                    # centroid of None by scikit-learn. The true entropy here
                    # is 0.0 for all sites, so the true centroid is 0.0
                    for s in all_subs:
                        if s.centroid == None:
                            s.centroid = [0.0]
                            log.debug("Fixed a subset with a centroid of None")
                            log.debug("The subset has %d columns" % len(s.columns))

                    s = fabricated_subsets.pop(0)

                    log.debug("Working on fabricated subset %s with %d sites" %(s.subset_id, len(s.columns)))
                    log.info("Finalising subset %d", i)
                    i = i+1

                    all_subs.remove(s)

                    centroid = s.centroid

                    best_match = None

                    # get closest subset to s
                    for sub in all_subs:

                        centroid_array = [sub.centroid, centroid]

                        euclid_dist = spatial.distance.pdist(centroid_array)

                        if euclid_dist < best_match or best_match is None:
                            best_match = euclid_dist
                            closest_sub = sub

                    # join s with closest_sub to make joined_sub
                    merged_sub = subset_ops.merge_subsets([s, closest_sub])

                    # remove closest sub
                    all_subs.remove(closest_sub)

                    # and if closest_sub was fabricated too, we remove it here
                    if fabricated_subsets.count(closest_sub):
                        fabricated_subsets.remove(closest_sub)

                    # analyse joined sub
                    self.analyse_list_of_subsets([merged_sub])

                    # here we put a sensible lower limit on the size of subsets
                    if len(merged_sub.columns)<the_config.min_subset_size:
                        merged_sub.fabricated = True

                    # if joined has to be fabricated, add to fabricated list
                    if merged_sub.fabricated:
                        fabricated_subsets.append(merged_sub)

                    all_subs.append(merged_sub)
        else:
            all_subs = start_subsets

        # now build a scheme from start_subs, and it should work
        final_scheme = scheme.Scheme(the_config, "final_scheme", all_subs)

        # return final scheme
        return final_scheme

    def build_new_subset_list(self, name_prefix, split_subs, start_subsets):
        new_scheme_subs = []
        for i, sub in enumerate(start_subsets):
            if len(sub.columns) == 1:
                new_scheme_subs.append(sub)
                log.debug("Split %d: parent subset has only one site, %s unchanged" %
                         (i+1, the_config.model_selection.upper()))

            elif sub.fabricated or sub.dont_split:
                new_scheme_subs.append(sub)
                log.debug("Split %d: parent subset couldn't be analysed, %s unchanged" %
                         (i+1, the_config.model_selection.upper()))

            else:  # compare split to un-split
                log.debug("Splitting new subset")

                # get list of split subsets from dictionary
                split_subsets = split_subs[sub]

                log.debug("# subs in this split: %d" % len(split_subsets))
                log.debug("dont_split: %s" % split_subsets[0].dont_split)

                # split subsets might be fabricated (i.e. unanalyseable)
                fabrications = 0
                for s in split_subsets:
                    fabrications = fabrications + s.fabricated

                if fabrications == 0:

                    score_diff = subset_ops.subset_list_score_diff(split_subsets, [sub], the_config, self.alignment)

                    log.info("Subset %d: split changed %s by: %.1f" %
                             (i+1, the_config.model_selection.upper(),
                              score_diff))

                    lnL, sum_k, subs_len = subset_ops.subset_list_stats([sub], the_config, self.alignment)

                    per_site_improvement = score_diff / subs_len

                    log.debug("Per site improvement: %.1f" % (per_site_improvement))

                    if score_diff < 0:
                        # We ONLY split the subset if the score improved and the LRT is significant
                        new_scheme_subs = new_scheme_subs + split_subsets
                    else:
                        sub.dont_split = True
                        new_scheme_subs.append(sub)
                elif fabrications == len(split_subsets):
                    # none of the split subsets worked, so don't analyse the parent again
                    sub.dont_split = True
                    new_scheme_subs.append(sub)
                    log.info("No splittable subsets")
                else:
                    new_scheme_subs = new_scheme_subs + split_subsets

        return new_scheme_subs

    def report(self, step):
        # Since the AIC will likely be better before we dealt with the
        # fabricated subsets, we need to set the best scheme and best result
        # to those from the last merged_scheme. TODO: add a variable to scheme
        # to take care of this problem so that the best AND analysable scheme
        # is the one that gets automatically flagged as the best scheme
        log.info("** Kmeans algorithm finished after %d steps **" % (step))
        log.info("Best scoring scheme has %d subsets and a %s score of %.3f"
                 % (
        len(self.results.best_scheme.subsets), the_config.model_selection,
        self.results.best_score))
        the_config.reporter.write_best_scheme(self.results)

        log.warning("Warning as of April 2016: We have noticed that the kmeans \
            algorithm does not perform well on some simulated datasets. \
            We are working on investigating and addressing this \
            but in the mean time we suggest being very cautious about using \
            this algorithm. At the very least, you should try other approaches \
            (e.g. partitioning by locus), and investigate your answers carefully \
            (both the trees and the partitioning schemes). If you have any \
            questions, please get in touch on the google group. Note that this \
            warning does not apply to cases where you are using models that have \
            an ascertainment bias for datasets that include only variable sites \
            as is often the case with morphological analyses."
            )



    def setup(self):

        log.warning("Warning as of April 2016: We have noticed that the kmeans \
            algorithm does not perform well on some simulated datasets. \
            We are working on investigating and addressing this \
            but in the mean time we suggest being very cautious about using \
            this algorithm. At the very least, you should try other approaches \
            (e.g. partitioning by locus), and investigate your answers carefully \
            (both the trees and the partitioning schemes). If you have any \
            questions, please get in touch on the google group. Note that this \
            warning does not apply to cases where you are using models that have \
            an ascertainment bias for datasets that include only variable sites \
            as is often the case with morphological analyses."
            )


        # set the default subset size to 100 for kmeans analyses
        if the_config.min_subset_size == False:
            the_config.min_subset_size = 100

        partnum = len(the_config.user_subsets)
        the_config.progress.begin(1, 1)

        # Start with the most partitioned scheme
        start_description = range(partnum)
        start_scheme = scheme.create_scheme(
            the_config, "start_scheme", start_description)

        site_max = sum([ len(s.columns) for s in start_scheme.subsets])

        if the_config.min_subset_size > site_max:
            log.error("The minimum subset size must be smaller than the \
                total number of sites you want to analyse. Your minimum \
                subset size is %d, and your alignment is %d sites. Please \
                check and try again." %(the_config.min_subset_size, site_max)
                )
            raise AnalysisError


        with logtools.indented(log, "**Analysing starting scheme (scheme %s)**" % start_scheme.name):
            start_result = self.analyse_scheme(start_scheme)

            if not the_config.quick:
                the_config.reporter.write_scheme_summary(start_scheme, start_result)

            tree_path = the_config.processor.make_tree_path(
                self.filtered_alignment_path)

        if the_config.kmeans == 'tiger' and the_config.datatype != 'morphology':
            log.error("You have selected kmeans and tiger \
                rates. This is an unsupported option for anything except \
                morphological data. The kmeans algorithm \
                now works with entropies, not TIGER rates.")
            raise AnalysisError

        return start_result, start_scheme, tree_path

    def one_kmeans_step(self, start_subsets, step, tree_path):
        name_prefix = "step_%d" % (step)

        # 1. Make split subsets

        with logtools.indented(log, "Splitting subsets using k-means"):
            split_subs = self.split_subsets(start_subsets, tree_path)

            # 2. Analyse split subsets (this to take advantage of parallelisation)
            subs = []

            # make a list from the dictionary
            for vals in split_subs.values():
                subs.extend(vals)


        log.debug("%d subsets successfully split" %(len(subs) - len(start_subsets)))

        with logtools.indented(log, "Calculating scores of all new subsets that can be analysed"):
            self.analyse_list_of_subsets(subs)

            # 3. Build new list of subsets
            new_scheme_subs = self.build_new_subset_list(name_prefix, split_subs, start_subsets)


        # 4. Are we done yet?
        if len(new_scheme_subs) == len(list(start_subsets)):
            log.info("""Could not improve %s score. Kmeans algorithm finished."""
                     % (the_config.model_selection))
            done = True
        else:
            n_splits = len(new_scheme_subs) - len(start_subsets)

            if n_splits > 1:
                t = 'subsets'
            else:
                t = 'subset'
            log.info("""The %s score of %d %s
                     improved when split"""
                     % (the_config.model_selection, n_splits, t))

            start_subsets = new_scheme_subs

            done = False

        return done, start_subsets

    def reassign_invariant_sites(self, subsets):

        #TODO add a skip:
        #if(len(subsets)==1):
        #   return(subsets)

        # get entropies for whole alignment for this subset
        onesub = subset_ops.merge_subsets(subsets)
        entropies = entropy.sitewise_entropies(SubsetAlignment(self.alignment, onesub))

        # find nearest site for each invariant site
        # replacements is a dict of: key: invariant col; value: replacement col,
        # e.g.
        # {512: 513, 514: 513, 515: 513, 516: 517}
        replacements = entropy.get_replacement_sites(entropies, onesub.columns)

        # now make a dict of the CURRENT subsets: key: site; value: subset
        sch_dict = {}
        for i, sub in enumerate(subsets):
            for site in sub.columns:
                sch_dict[site] = i

        # then reassign the sites as necessary based on replacements
        for r in replacements:
            sch_dict[r] = sch_dict[replacements[r]]

        # now build subsets according to the new sites
        sub_dict = {} # this gives us the subsets to build
        for k, v in sch_dict.iteritems():
            sub_dict.setdefault(v, []).append(k)

        new_subsets = []
        for s in sub_dict:
            n = Subset(the_config, set(sub_dict[s]))
            new_subsets.append(n)

        return(new_subsets)


    @logtools.log_info(log, "Performing k-means Analysis")
    def do_analysis(self):
        '''A kmeans algorithm for heuristic partitioning searches'''

        start_result, start_scheme, tree_path = self.setup()

        step = 0

        start_subsets = list(start_scheme.subsets) # we only work on lists of subsets

        self.analyse_list_of_subsets(start_subsets)

        # now we suppress ExternalProgramError for the rest of the algorithm
        the_config.suppress_errors = True

        for s in start_subsets:
            if s.fabricated:
                log.error("""One or more of your starting datablocks could not
                          be analysed. Please check your data and try again.
                          One way to fix this is to join your small datablocks
                          together into larger datablocks""")
                raise AnalysisError

        while True:
            step += 1
            with logtools.indented(log, "***k-means algorithm step %d***"
                    % step):
                done, start_subsets = self.one_kmeans_step(
                    start_subsets, step, tree_path)

            if done:
                break

        # Ok, we're done. we just need deal with fabricated subsets
        final_scheme = self.finalise_fabrication(start_subsets, step)

        # Finally, for krmeans, we put the invariant sites back with their
        # nearest variable neighbours
        if the_config.search == 'krmeans':
            log.info("Reassigning invariant sites for krmeans algorithm")
            # the definition of krmeans is that we reassign the zero entropies
            final_subsets = self.reassign_invariant_sites(final_scheme.subsets)
            final_scheme = scheme.Scheme(the_config, "final_scheme_reassigned", final_subsets)

        log.info("Analysing final scheme")

        final_result = self.analyse_scheme(final_scheme)

        self.report(step)

        if not the_config.quick:
            the_config.reporter.write_scheme_summary(final_scheme, final_result)


        return(final_scheme)


def choose_method(search):
    if search == 'all':
        method = AllAnalysis
    elif search == 'user':
        method = UserAnalysis
    elif search == 'greedy':
        method = GreedyAnalysis
    elif search == 'hcluster':
        method = StrictClusteringAnalysis
    elif search == 'rcluster':
        method = RelaxedClusteringAnalysis
    elif search == 'rclusterf':
        method = RelaxedClusteringAnalysis
    elif search == 'kmeans' or search == 'krmeans':
        method = KmeansAnalysis
    else:
        log.error("Search algorithm '%s' is not recognised", search)
        raise AnalysisError
    return method
