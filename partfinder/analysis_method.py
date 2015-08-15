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
log = logtools.get_logger()

import math
import scheme
import submodels
from analysis import Analysis, AnalysisError
import neighbour
import kmeans
import itertools
import subset_ops
from scipy import spatial, exp
from scipy.misc import comb
import numpy as np
from config import the_config

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

        # Start with the most partitioned scheme
        start_description = range(partnum)
        start_scheme = scheme.create_scheme(
            the_config, "start_scheme", start_description)

        with logtools.indented(log, "Analysing starting scheme (scheme %s)" %
                          start_scheme.name):
            start_result = self.analyse_scheme(start_scheme)

            if not the_config.quick:
                the_config.reporter.write_scheme_summary(start_scheme, start_result)

        step = 1

        # Now we try out all lumpings of the current scheme, to see if we can
        # find a better one and if we do, we just keep going
        while True:
            with logtools.indented(log, "***Greedy algorithm step %d***" % step):
                name_prefix = "step_%d" % (step)
                old_best_score = self.results.best_score

                # Make a list of all the new subsets, and get them analysed
                # We do them in blocks of 10K, to avoid memory overload
                lumped_subset_iterator = itertools.combinations(start_scheme.subsets, 2)
                new_subs = []
                log.info("Building subsets")
                for subset_grouping in lumped_subset_iterator:
                    new_sub = subset_ops.merge_subsets(subset_grouping)
                    if not new_sub.is_done:
                        new_subs.append(new_sub)
                    if len(new_subs)>9999:
                        log.info("Analysing 10,000 subsets")
                        self.analyse_list_of_subsets(new_subs)
                        new_subs = []
                        log.info("Building more subsets")

                # analyse what's left, and clean out list
                log.info("Analysing %d subsets"%(len(new_subs)))
                self.analyse_list_of_subsets(new_subs)
                new_subs = []

                log.info("Analysing schemes")
                # we repeat the iterator, for memory efficiency
                lumped_subset_iterator = itertools.combinations(start_scheme.subsets, 2)
                sch_num = 1
                for subset_grouping in lumped_subset_iterator:
                    # could do this without another merge, but this seems most robust
                    new_sub = subset_ops.merge_subsets(subset_grouping)
                    scheme_name = "%s_%d" % (name_prefix, sch_num)
                    lumped_scheme = neighbour.make_clustered_scheme(
                        start_scheme, scheme_name, subset_grouping, new_sub, the_config)
                    sch_num = sch_num + 1

                    new_result = self.analyse_scheme(lumped_scheme)
                    log.debug("Difference in %s: %.1f" %
                            (the_config.model_selection,
                            (new_result.score - old_best_score)))

                if self.results.best_score != old_best_score:
                    log.info("""Analysed all schemes for this step. The best
                        scheme changed the %s score by %.1f units.""" % (
                        the_config.model_selection,
                        self.results.best_score - old_best_score
                    ))

                    self.results.best_scheme.name = "step_%d" % step

                    if not the_config.quick:
                        the_config.reporter.write_scheme_summary(
                            self.results.best_scheme, self.results.best_result)

                    # Now we find out which is the best lumping we know of for
                    # this step
                    start_scheme = self.results.best_scheme
                else:
                    log.info("""Analysed all schemes for this step and found no
                        schemes that improve the score, stopping""")
                    break

            # We're done if it's the scheme with everything together
            if len(set(lumped_scheme.subsets)) == 1:
                break

            step += 1

        log.info("Greedy algorithm finished after %d steps" % step)
        log.info("Best scoring scheme is scheme %s, with %s score of %.3f"
                 % (self.results.best_scheme.name, the_config.model_selection,
                    self.results.best_score))

        the_config.reporter.write_best_scheme(self.results)


class RelaxedClusteringAnalysis(Analysis):
    '''
    A relaxed clustering algorithm for heuristic partitioning searches

    1. Rank subsets by their similarity (defined by clustering-weights)
    2. Analyse min(cluster-percent or cluster-max) most similar schemes
    3. Take the scheme that improves the AIC/BIC score the most
    4. Quit if no improvements.
    '''

    def do_analysis(self):
        log.info("Performing relaxed clustering analysis")

        # initialisation steps
        model_selection = the_config.model_selection
        partnum = len(the_config.user_subsets)

        scheme_count = submodels.count_relaxed_clustering_schemes(
            partnum, the_config.cluster_percent, the_config.cluster_max)
        subset_count = submodels.count_relaxed_clustering_subsets(
            partnum, the_config.cluster_percent, the_config.cluster_max)

        log.info("PartitionFinder will have to analyse %d subsets to"
                 " complete this analyses" % subset_count)
        the_config.progress.begin(scheme_count, subset_count)

        # Start with the most partitioned scheme, and record it.
        log.info("*** Analysing starting scheme ***")
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
        while True:
            log.info("*** Relaxed clustering algorithm step %d of up to %d ***"
                % (step, partnum - 1))

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
            if the_config.cluster_max != None and cutoff>the_config.cluster_max:
                cutoff = the_config.cluster_max
            log.debug("Choosing the %d most similar subset pairs" % cutoff)
            closest_pairs = neighbour.get_N_closest_subsets(
                subsets, the_config, cutoff, d_matrix)

            # 2. analyse K subsets in top N that have not yet been analysed
            pairs_todo = neighbour.get_pairs_todo(closest_pairs, c_matrix, subsets)
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
            best_change = np.amin(c_matrix)

            best_pair = neighbour.get_best_pair(c_matrix, best_change, subsets)

            best_merged = subset_ops.merge_subsets(best_pair)
            best_scheme = neighbour.make_clustered_scheme(
                start_scheme, scheme_name, best_pair, best_merged, the_config)
            best_result = self.analyse_scheme(best_scheme)

            # the best change can get updated a fraction at this point
            best_change = self.results.best_score - start_score

            if best_change>=0:
                log.info("Found no schemes that improve the score, stopping")
                break

            log.info("The best scheme improves the %s score by %.2f to %.1f",
                the_config.model_selection,
                np.abs(best_change),
                self.results.best_score)
            start_scheme = best_scheme
            start_score = best_result.score


            # 5. reset_c_matrix and the subset list
            c_matrix = neighbour.reset_c_matrix(c_matrix, list(best_pair), [best_merged], subsets)
            subsets = neighbour.reset_subsets(subsets, list(best_pair), [best_merged])

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



        the_config.reporter.write_best_scheme(self.results)

class HybridAnalysis(Analysis):
    def split_subsets(self, start_subsets, tree_path):
        split_subs = {}
        for sub in start_subsets:
            if len(sub.columns) == 1:
                split_subs[sub] = [sub]
            else:
                split = kmeans.kmeans_split_subset(
                    the_config, self.alignment, sub, tree_path, n_jobs=self.threads)

                if split == 1:  # we couldn't analyse the big subset
                    sub.dont_split = True # never try to split this subset again
                    split_subs[sub] = [sub]  # so we keep it whole
                else:  # we could analyse the big subset
                    split_subs[sub] = split  # so we split it into >1

        return split_subs

    def finalise_fabrication(self, start_subsets, step):

        fabricated_subsets = []
        for s in start_subsets:
            if s.fabricated:
                fabricated_subsets.append(s)

        if fabricated_subsets:
            log.info("""Finalising partitioning scheme, by incorporating
                     the subsets that couldn't be analysed with their
                     nearest neighbours""")
            log.debug("There are %d/%d fabricated subsets"
                      % (len(fabricated_subsets), len(start_subsets)))

            while fabricated_subsets:

                all_subs = start_subsets
                s = fabricated_subsets.pop(0)
                all_subs.remove(s)

                centroid = s.centroid
                best_match = None

                # get closest subset to s
                for sub in all_subs:
                    centroid_array = [sub.centroid, centroid]

                    print centroid_array

                    euclid_dist = spatial.distance.pdist(centroid_array)

                    if euclid_dist < best_match or best_match is None:
                        best_match = euclid_dist
                        closest_sub = sub

                # join s with closest_sub to make joined_sub
                merged_sub = subset_ops.merge_subsets([s, closest_sub])

                # remove closest sub
                all_subs.remove(closest_sub)

                # analyse joined sub
                self.analyse_list_of_subsets([merged_sub])

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
        for sub in start_subsets:
            if len(sub.columns) == 1 or sub.fabricated or sub.dont_split:
                new_scheme_subs.append(sub)

            else:  # compare split to un-split
                # get list of split subsets from dictionary
                split_subsets = split_subs[sub]


                # split subsets might be fabricated (i.e. unanalyseable)
                fabrications = 0
                for s in split_subsets:
                    fabrications = fabrications + s.fabricated

                if fabrications == 0:
                    split_score = subset_ops.subset_list_score(split_subsets, the_config, self.alignment)
                    unsplit_score = subset_ops.subset_list_score([sub], the_config, self.alignment)

                    score_diff = split_score - unsplit_score
                    log.debug("Difference in %s: %.1f" %
                             (the_config.model_selection.upper(),
                              score_diff))
                    if score_diff < 0:
                        new_scheme_subs = new_scheme_subs + split_subsets
                    else:
                        sub.dont_split = True
                        new_scheme_subs.append(sub)
                elif fabrications == len(split_subsets):
                    # none of the split subsets worked, so don't analyse the parent again
                    sub.dont_split = True
                    new_scheme_subs.append(sub)
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

    def setup(self):
        partnum = len(the_config.user_subsets)
        the_config.progress.begin(1, 1)

        # Start with the most partitioned scheme
        start_description = range(partnum)
        start_scheme = scheme.create_scheme(
            the_config, "start_scheme", start_description)

        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        start_result = self.analyse_scheme(start_scheme)

        the_config.reporter.write_scheme_summary(start_scheme, start_result)

        tree_path = the_config.processor.make_tree_path(
            self.filtered_alignment_path)

        return start_result, start_scheme, tree_path

    def one_kmeans_step(self, start_subsets, step, tree_path):
        name_prefix = "step_%d" % (step)

        # 1. Make split subsets
        log.info("Splitting %d subsets using K-means" % len(start_subsets))
        split_subs = self.split_subsets(start_subsets, tree_path)

        # 2. Analyse split subsets (this to take advantage of parallelisation)
        subs = []

        # make a list from the dictionary
        for vals in split_subs.values():
            subs.extend(vals)

        self.analyse_list_of_subsets(subs)

        # 3. Build new list of subsets
        new_scheme_subs = self.build_new_subset_list(name_prefix, split_subs, start_subsets)


        # 4. Are we done yet?
        if len(new_scheme_subs) == len(list(start_subsets)):
            log.info("""The %s score of 0 subsets
                     improved when split. Algorithm finished."""
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

    def do_analysis(self):
        '''A kmeans algorithm for heuristic partitioning searches'''

        log.info("Performing hybrid analysis")
        log.info("This will run a K-means analysis "
                 "followed by a relaxed clustering analysis")

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

        log.info("Analysing final scheme")
        final_result = self.analyse_scheme(final_scheme)

        the_config.user_subsets = list(final_scheme.subsets)

        log.info("Performing relaxed clustering analysis")

        # initialisation steps
        model_selection = the_config.model_selection
        partnum = len(the_config.user_subsets)

        scheme_count = submodels.count_relaxed_clustering_schemes(
            partnum, the_config.cluster_percent, the_config.cluster_max)
        subset_count = submodels.count_relaxed_clustering_subsets(
            partnum, the_config.cluster_percent, the_config.cluster_max)

        log.info("PartitionFinder will have to analyse %d subsets to"
                 " complete this analyses" % subset_count)
        the_config.progress.begin(scheme_count, subset_count)

        # Start with the most partitioned scheme, and record it.
        log.info("*** Analysing starting scheme ***")
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
        while True:
            log.info("*** Relaxed clustering algorithm step %d of up to %d ***"
                % (step, partnum - 1))

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
            if the_config.cluster_max != None and cutoff>the_config.cluster_max:
                cutoff = the_config.cluster_max
            log.debug("Choosing the %d most similar subset pairs" % cutoff)
            closest_pairs = neighbour.get_N_closest_subsets(
                subsets, the_config, cutoff, d_matrix)

            # 2. analyse K subsets in top N that have not yet been analysed
            pairs_todo = neighbour.get_pairs_todo(closest_pairs, c_matrix, subsets)
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
            best_change = np.amin(c_matrix)

            best_pair = neighbour.get_best_pair(c_matrix, best_change, subsets)

            best_merged = subset_ops.merge_subsets(best_pair)
            best_scheme = neighbour.make_clustered_scheme(
                start_scheme, scheme_name, best_pair, best_merged, the_config)
            best_result = self.analyse_scheme(best_scheme)

            # the best change can get updated a fraction at this point
            best_change = self.results.best_score - start_score

            if best_change>=0:
                log.info("Found no schemes that improve the score, stopping")
                break

            log.info("The best scheme improves the %s score by %.2f to %.1f",
                the_config.model_selection,
                np.abs(best_change),
                self.results.best_score)
            start_scheme = best_scheme
            start_score = best_result.score

            if not the_config.quick:
                the_config.reporter.write_scheme_summary(
                    best_scheme, best_result)

            # 5. reset_c_matrix and the subset list
            c_matrix = neighbour.reset_c_matrix(c_matrix, list(best_pair), [best_merged], subsets)
            subsets = neighbour.reset_subsets(subsets, list(best_pair), [best_merged])

            if len(set(start_scheme.subsets)) == 1:
                break

            step += 1

        log.info("Relaxed clustering algorithm finished after %d steps" % step)
        log.info("Best scoring scheme is scheme %s, with %s score of %.3f"
                 % (self.results.best_scheme.name, model_selection,
                    self.results.best_score))



        the_config.reporter.write_best_scheme(self.results)



class KmeansAnalysis(Analysis):


    def split_subsets(self, start_subsets, tree_path):
        split_subs = {}
        for sub in start_subsets:
            if len(sub.columns) == 1:
                split_subs[sub] = [sub]
            else:
                split = kmeans.kmeans_split_subset(
                    the_config, self.alignment, sub, tree_path, n_jobs=self.threads)


                if split == 1:  # we couldn't analyse the big subset
                    sub.dont_split = True # never try to split this subset again
                    split_subs[sub] = [sub]  # so we keep it whole
                elif len(split) == 1:
                    # in some cases (i.e. all site params are equal) kmeans
                    # cannot split subsets, so we get back the same as we put in
                    sub.dont_split = True # never try to split this subset again
                    split_subs[sub] = [sub]  # so we keep it whole
                else:  # we could analyse the big subset
                    split_subs[sub] = split  # so we split it into >1
                    log.info("..subset of %d sites split into %d and %d sites"
                            %(len(sub.columns), len(split[0].columns), len(split[1].columns)))

        return split_subs

    def finalise_fabrication(self, start_subsets, step):

        fabricated_subsets = []
        for s in start_subsets:
            if s.fabricated:
                fabricated_subsets.append(s)

        if fabricated_subsets:
            log.info("""Finalising partitioning scheme, by incorporating
                     the subsets that couldn't be analysed with their
                     nearest neighbours""")
            log.debug("There are %d/%d fabricated subsets"
                      % (len(fabricated_subsets), len(start_subsets)))

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

                # analyse joined sub
                self.analyse_list_of_subsets([merged_sub])

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
        for sub in start_subsets:
            if len(sub.columns) == 1 or sub.fabricated or sub.dont_split:
                new_scheme_subs.append(sub)

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

                    log.info("..split changed %s by: %.1f" %
                             (the_config.model_selection.upper(),
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

    def setup(self):
        partnum = len(the_config.user_subsets)
        the_config.progress.begin(1, 1)

        # Start with the most partitioned scheme
        start_description = range(partnum)
        start_scheme = scheme.create_scheme(
            the_config, "start_scheme", start_description)

        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        start_result = self.analyse_scheme(start_scheme)

        the_config.reporter.write_scheme_summary(start_scheme, start_result)

        tree_path = the_config.processor.make_tree_path(
            self.filtered_alignment_path)

        if the_config.kmeans == 'tiger':
            try:
                from _tiger import TigerDNA
                the_config.TigerDNA = TigerDNA
            except:
                log.error("Couldn't find compiled tiger code.")
                log.error("You have selected kmeans and tiger \
                    rates. This is an unsupported option, if you still wish to use \
                    this option, you must compile the tiger code.")
                log.error("Once you compile the tiger code, this option will work. \
                    But please note that this is an \
                    unsupported option. For empirical work we recommend using \
                    entropy calculations for site rates, which is the default \
                    behaviour for the kmeans algorithm in PF2.")
                raise AnalysisError
        else:
            the_config.TigerDNA = None


        return start_result, start_scheme, tree_path

    def one_kmeans_step(self, start_subsets, step, tree_path):
        name_prefix = "step_%d" % (step)

        # 1. Make split subsets
        log.info("Splitting subsets using k-means")
        split_subs = self.split_subsets(start_subsets, tree_path)

        # 2. Analyse split subsets (this to take advantage of parallelisation)
        subs = []

        # make a list from the dictionary
        for vals in split_subs.values():
            subs.extend(vals)

        log.info("%d subsets successfully split" %(len(subs) - len(start_subsets)))

        log.info("Comparing scores of split to unsplit subsets")
        self.analyse_list_of_subsets(subs)

        # 3. Build new list of subsets
        new_scheme_subs = self.build_new_subset_list(name_prefix, split_subs, start_subsets)


        # 4. Are we done yet?
        if len(new_scheme_subs) == len(list(start_subsets)):
            log.info("""The %s score of 0 subsets
                     improved when split. Algorithm finished."""
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

        log.info("Analysing final scheme")
        final_result = self.analyse_scheme(final_scheme)

        self.report(step)

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
    elif search == 'kmeans':
        method = KmeansAnalysis
    elif search == 'hybrid':
        method = HybridAnalysis
    else:
        log.error("Search algorithm '%s' is not yet implemented", search)
        raise AnalysisError
    return method
