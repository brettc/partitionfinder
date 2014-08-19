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

log = logtools.get_logger(__file__)

import math
import scheme
import submodels
from analysis import Analysis, AnalysisError
import neighbour
import kmeans
import itertools
import subset_ops
from scipy import spatial
from scipy.misc import comb
import warnings
import numpy as np

from util import PhylogenyProgramError


class UserAnalysis(Analysis):
    def do_analysis(self):
        log.info("Performing User analysis")
        current_schemes = self.cfg.user_schemes
        scheme_count = len(current_schemes)
        subset_count = len(self.cfg.user_subsets)

        self.cfg.progress.begin(scheme_count, subset_count)
        if scheme_count > 0:
            for s in current_schemes:
                res = self.analyse_scheme(s)

                # Write out the scheme
                self.cfg.reporter.write_scheme_summary(s, res)
        else:
            log.error(
                "Search set to 'user', but no user schemes detected in .cfg file. Please check.")
            raise AnalysisError

        self.cfg.progress.end()

        self.cfg.reporter.write_best_scheme(self.results)


class StrictClusteringAnalysis(Analysis):
    """
    This analysis uses model parameters to guess at similar partitions, then
    just joins them together this is much less accurate than other methods, but
    a LOT quicker - it runs in order N time (where N is the number of initial
    datablocks), whereas the greedy algorithm is still N squared.
    """

    def do_analysis(self):
        log.info("Performing strict clustering analysis")

        partnum = len(self.cfg.user_subsets)
        subset_count = 2 * partnum - 1
        scheme_count = partnum
        self.cfg.progress.begin(scheme_count, subset_count)

        # Start with the most partitioned scheme
        start_description = range(partnum)
        start_scheme = scheme.create_scheme(
            self.cfg, "start_scheme", start_description)

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
                start_scheme, scheme_name, self.cfg)

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

        self.cfg.progress.end()
        self.cfg.reporter.write_best_scheme(self.results)


class AllAnalysis(Analysis):
    def do_analysis(self):
        log.info("Performing complete analysis")
        partnum = len(self.cfg.user_subsets)

        scheme_count = submodels.count_all_schemes(partnum)
        subset_count = submodels.count_all_subsets(partnum)
        self.cfg.progress.begin(scheme_count, subset_count)

        # Iterate over submodels, which we can turn into schemes afterwards in the loop
        model_iterator = submodels.submodel_iterator([], 1, partnum)

        scheme_name = 1
        for m in model_iterator:
            s = scheme.model_to_scheme(m, scheme_name, self.cfg)
            scheme_name = scheme_name + 1
            res = self.analyse_scheme(s)

            # Write out the scheme
            self.cfg.reporter.write_scheme_summary(s, res)

        self.cfg.reporter.write_best_scheme(self.results)


class GreedyAnalysis(Analysis):

    @logtools.log_info(log, "Performing Greedy Analysis")
    def do_analysis(self):
        '''A greedy algorithm for heuristic partitioning searches'''

        partnum = len(self.cfg.user_subsets)
        scheme_count = submodels.count_greedy_schemes(partnum)
        subset_count = submodels.count_greedy_subsets(partnum)

        self.cfg.progress.begin(scheme_count, subset_count)

        # Start with the most partitioned scheme
        start_description = range(partnum)
        start_scheme = scheme.create_scheme(
            self.cfg, "start_scheme", start_description)

        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        log.push()
        start_result = self.analyse_scheme(start_scheme)
        self.cfg.reporter.write_scheme_summary(start_scheme, start_result)
        log.pop()

        step = 1

        # Now we try out all lumpings of the current scheme, to see if we can
        # find a better one and if we do, we just keep going
        while True:
            log.info("***Greedy algorithm step %d***" % step)
            log.push()
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
                    start_scheme, scheme_name, subset_grouping, new_sub, self.cfg)
                sch_num = sch_num + 1

                new_result = self.analyse_scheme(lumped_scheme)
                log.debug("Difference in %s: %.1f" %
                          (self.cfg.model_selection,
                           (new_result.score - old_best_score)))

            if self.results.best_score != old_best_score:
                log.info("""Analysed all schemes for this step. The best
                    scheme changed the %s score by %.1f units.""" % (
                    self.cfg.model_selection,
                    self.results.best_score - old_best_score
                ))

                self.results.best_scheme.name = "step_%d" % step
                self.cfg.reporter.write_scheme_summary(
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

            log.pop()
            step += 1

        #log.pop()
        log.info("Greedy algorithm finished after %d steps" % step)
        log.info("Best scoring scheme is scheme %s, with %s score of %.3f"
                 % (self.results.best_scheme.name, self.cfg.model_selection,
                    self.results.best_score))

        self.cfg.reporter.write_best_scheme(self.results)
        #log.pop()


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
        model_selection = self.cfg.model_selection
        partnum = len(self.cfg.user_subsets)
        nseq = len(self.alignment.species)

        scheme_count = submodels.count_relaxed_clustering_schemes(
            partnum, self.cfg.cluster_percent, self.cfg.cluster_max)
        subset_count = submodels.count_relaxed_clustering_subsets(
            partnum, self.cfg.cluster_percent, self.cfg.cluster_max)

        log.info("PartitionFinder will have to analyse %d subsets to"
                 " complete this analyses" % subset_count)
        self.cfg.progress.begin(scheme_count, subset_count)

        # Start with the most partitioned scheme, and record it.
        log.info("*** Analysing starting scheme ***")
        self.cfg.progress.begin(scheme_count, partnum)
        start_scheme = scheme.create_scheme(
            self.cfg, "start_scheme", range(partnum))
        start_result = self.analyse_scheme(start_scheme)
        start_score = start_result.score
        if not self.cfg.quick:
            self.cfg.reporter.write_scheme_summary(
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
                self.cfg.cluster_weights)

            if step == 1:
                # Now initialise a change in info score matrix to inf
                c_matrix = np.empty(d_matrix.shape)
                c_matrix[:] = np.inf
                c_matrix = spatial.distance.squareform(c_matrix)

            # 1. pick top N subset pairs from distance matrix
            cutoff = int(math.ceil(max_schemes * (self.cfg.cluster_percent * 0.01)))
            if cutoff <= 0: cutoff = 1
            if self.cfg.cluster_max != None and cutoff>self.cfg.cluster_max:
                cutoff = self.cfg.cluster_max
            log.debug("Choosing the %d most similar subset pairs" % cutoff)
            closest_pairs = neighbour.get_N_closest_subsets(
                subsets, self.cfg, cutoff, d_matrix)

            # 2. analyse K subsets in top N that have not yet been analysed
            pairs_todo = neighbour.get_pairs_todo(closest_pairs, c_matrix, subsets)
            log.info("Analysing %d new subset pairs" % len(pairs_todo))
            new_subs = []
            sub_tuples = []
            for pair in pairs_todo:
                new_sub = subset_ops.merge_subsets(pair)
                new_subs.append(new_sub)
                sub_tuples.append((new_sub, pair))

            self.cfg.progress.begin(scheme_count, len(new_subs))
            self.analyse_list_of_subsets(new_subs)

            # 3. for all K new subsets, update improvement matrix and find best pair
            log.info("Finding the best partitioning scheme")
            diffs = []
            min_diff = np.inf
            scheme_name = "step_%d" %(step)
            for t in sub_tuples:
                pair_merged = t[0]
                pair = t[1]
                new_scheme = neighbour.make_clustered_scheme(
                        start_scheme, scheme_name, pair, pair_merged, self.cfg)
                r = self.analyse_scheme(new_scheme)
                diff = r.score - start_score
                diffs.append(diff)

            c_matrix = neighbour.update_c_matrix(c_matrix, sub_tuples, subsets, diffs)

            # 4. Find the best pair of subsets, and build a scheme based on that
            best_change = np.amin(c_matrix)

            best_pair = neighbour.get_best_pair(c_matrix, best_change, subsets)

            best_merged = subset_ops.merge_subsets(best_pair)
            best_scheme = neighbour.make_clustered_scheme(
                start_scheme, scheme_name, best_pair, best_merged, self.cfg)                
            best_result = self.analyse_scheme(best_scheme)

            # the best change can get updated a fraction at this point
            best_change = self.results.best_score - start_score

            if best_change>=0:
                log.info("Found no schemes that improve the score, stopping")
                break

            log.info("The best scheme improves the %s score by %.2f to %.1f",
                self.cfg.model_selection, 
                np.abs(best_change),
                self.results.best_score)
            start_scheme = best_scheme
            start_score = best_result.score

            if not self.cfg.quick:
                self.cfg.reporter.write_scheme_summary(
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



        self.cfg.reporter.write_best_scheme(self.results)


class KmeansAnalysis(Analysis):

    @logtools.log_info(log, "Performing k-means Analysis")
    def do_analysis(self):
        '''A greedy algorithm for heuristic partitioning searches'''

        partnum = len(self.cfg.user_subsets)
        self.cfg.progress.begin(1, 1)

        # Start with the most partitioned scheme
        start_description = range(partnum)
        start_scheme = scheme.create_scheme(
            self.cfg, "start_scheme", start_description)

        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        log.push()
        start_result = self.analyse_scheme(start_scheme)
        self.cfg.reporter.write_scheme_summary(start_scheme, start_result)
        log.pop()

        processor = self.cfg.processor
        alignment_path = self.filtered_alignment_path
        tree_path = processor.make_tree_path(alignment_path)

        step = 0
        while True:
            step += 1
            log.info("***k-means algorithm step %d***" % step)
            log.push()
            name_prefix = "step_%d" % (step)

            # 1. Make split subsets
            split_subs = {}
            for sub in start_scheme.subsets:
                if len(sub.columns) == 1:
                    split_subs[sub] = [sub]
                else:
                    split = kmeans.kmeans_split_subset(self.cfg,
                                                       self.alignment,
                                                       sub,
                                                       tree_path)
                    
                    if split == 1: # we couldn't analyse the big subset
                        split_subs[sub] = [sub] # so we keep it whole
                    else: # we could analyse the big subset
                        split_subs[sub] = split # so we split it into >1

            # 2. Analyse split subsets (this to take advantage of parallelisation)
            subs = []
            for key in split_subs:
                subs = subs + split_subs[key]
            self.analyse_list_of_subsets(subs)

            # 3. Build new scheme
            new_scheme_subs = [] 
            for sub in start_scheme.subsets:

                if len(sub.columns) == 1 or sub.fabricated:
                    new_scheme_subs.append(sub)
                else: # compare split to un-split
                    split_subsets = split_subs[sub]
                    split_scheme = neighbour.make_split_scheme(
                        start_scheme, name_prefix, sub, split_subsets, self.cfg)

                    new_result = self.analyse_scheme(split_scheme)
                    score_diff = new_result.score - start_result.score
                    log.info("Difference in %s: %.1f" %
                                (self.cfg.model_selection.upper(),
                                score_diff))
                    if score_diff < 0:
                        new_scheme_subs = new_scheme_subs + split_subsets
                    else:
                        new_scheme_subs.append(sub)


            if len(new_scheme_subs) == len(list(start_scheme.subsets)):
                log.info("""Analysed all subsets, but couldn't Find
                            a split that improved the score. Quitting.""")
                break
            else:
                n_splits = len(new_scheme_subs) - len(start_scheme.subsets)
                old_best_score = start_result.score
                start_scheme = scheme.Scheme(self.cfg, name_prefix, new_scheme_subs)
                start_result = self.analyse_scheme(start_scheme)


                log.info("""Analysed all subsets. Found %d subsets which can be split.
                    New scheme changes the %s score by %.1f units.""" % (
                    n_splits,
                    self.cfg.model_selection,
                    self.results.best_score - old_best_score
                ))

                self.cfg.reporter.write_scheme_summary(
                    start_scheme, start_result)


        fabricated_subsets = start_scheme.get_fabricated_subsets()

        if fabricated_subsets:
            log.info("Finalising partitioning scheme")
            log.debug("There are %d/%d fabricated subsets" 
                    %(len(fabricated_subsets), len(start_scheme.subsets)))

            while fabricated_subsets:

                s = fabricated_subsets[0]
                centroid = s.centroid
                best_match = None

                # get closest subset to s
                all_subs = list(start_scheme)
                all_subs.remove(s)
                for sub in all_subs:
                    centroid_array = [sub.centroid, centroid]
                    warnings.simplefilter('ignore', DeprecationWarning)
                    euclid_dist = spatial.distance.pdist(centroid_array)
                    if euclid_dist < best_match or best_match is None:
                        best_match = euclid_dist
                        closest_sub = sub

                # create, analyse, report the new scheme
                scheme_name = "step_%d" %(step)
                merged_sub = subset_ops.merge_fabricated_subsets([s, closest_sub])

                start_scheme = neighbour.make_clustered_scheme(
                    start_scheme, scheme_name, [s, closest_sub], merged_sub, self.cfg)
                start_result = self.analyse_scheme(start_scheme)
                self.cfg.reporter.write_scheme_summary(start_scheme, start_result)

                fabricated_subsets = start_scheme.get_fabricated_subsets()


        # Since the AIC will likely be better before we dealt with the
        # fabricated subsets, we need to set the best scheme and best result
        # to those from the last merged_scheme. TODO: add a variable to scheme
        # to take care of this problem so that the best AND analysable scheme
        # is the one that gets automatically flagged as the best scheme
        self.results.best_scheme = start_scheme
        self.results.best_result = start_result


        log.info("** Kmeans algorithm finished after %d steps **" % (step))
        log.info("Best scoring scheme has %d subsets and a %s score of %.3f"
                 % (len(self.results.best_scheme.subsets), self.cfg.model_selection,
                    self.results.best_score))

        self.cfg.reporter.write_best_scheme(self.results)




class oldKmeansAnalysis(Analysis):
    def do_analysis(self):
        # Copied and pasted from greedy analysis
        partnum = len(self.cfg.user_subsets)
        self.cfg.progress.begin(1, 1)

        # Start with the most partitioned scheme
        start_description = range(partnum)
        start_scheme = scheme.create_scheme(
            self.cfg, "start_scheme", start_description)

        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        old_score = self.analyse_scheme(start_scheme)

        # Get first scheme
        best_scheme = start_scheme
        subset_index = 0
        all_subsets = list(best_scheme.subsets)

        processor = self.cfg.processor
        alignment_path = self.filtered_alignment_path
        tree_path = processor.make_tree_path(alignment_path)
        best_result = self.analyse_scheme(best_scheme)
        fabricated_subsets = []
        step = 1

        while subset_index < len(all_subsets):
            log.info("Best scheme has %s score of %.2f and %d subset(s)"
                     % (self.cfg.model_selection.upper(), best_result.score,
                        len(best_scheme.subsets)))

            log.info("***Kmeans algorithm step %d***" % step)
            step += 1

            current_subset = all_subsets[subset_index]

            log.info("Analysing subset of %d sites" %
                     len(current_subset.columns))

            # First check if the subset is large enough to split, if it isn't,
            # move to the next subset
            if len(current_subset.columns) == 1:
                log.info("This subset cannot be split further")
                subset_index += 1
                continue

            if current_subset.fabricated:
                log.info(
                    "This subset cannot be split further because %s cannot analyse it",
                    self.cfg.phylogeny_program)
                subset_index += 1
                fabricated_subsets.append(current_subset)
                continue

            split_subsets = kmeans.kmeans_split_subset(self.cfg,
                                                       self.alignment,
                                                       current_subset,
                                                       tree_path)


            # kmeans_split_subset will return a 1 and flag the subset as
            # fabricated if for some reason it raises a PhylogenyProgramError,
            # this it to catch those fabricated subsets
            if split_subsets == 1:
                subset_index += 1
                fabricated_subsets.append(current_subset)
                continue


            # Take a copy
            updated_subsets = all_subsets[:]

            # Replace the current one with the split one
            # Google "slice assignments"
            # This list is the key to avoiding recursion. It expands to contain
            # all of the split subsets by replacing them with the split ones
            updated_subsets[subset_index:subset_index + 1] = split_subsets

            test_scheme = scheme.Scheme(self.cfg, str(step - 1),
                                        updated_subsets)

            new_result = self.analyse_scheme(test_scheme)

            if new_result.score < best_result.score:
                best_scheme = test_scheme
                best_result = new_result

                # Change this to the one with split subsets in it. Note that
                # the subset_index now points a NEW subset, one that was split
                all_subsets = updated_subsets

                # record each scheme that's an improvement
                self.cfg.reporter.write_scheme_summary(
                    self.results.best_scheme, self.results.best_result)

                if len(split_subsets) == 2:
                    log.info(
                        "Splitting subset into %d:%d sites improved the %s score"
                        % (len(split_subsets[0].columns),
                           len(split_subsets[1].columns),
                           self.cfg.model_selection))

                    for s in split_subsets:
                        m = [x % 3 for x in s.columns]
                        l = float(len(s.columns))
                        props = [(float(m.count(1)) / l),
                                 (float(m.count(2)) / l),
                                 (float(m.count(0)) / l)]
                        log.info("%d subset has 1st, 2nd, 3rd props: %s" % (
                        len(s.columns), str(props)))



            else:
                log.info("Splitting this subset did not improve the %s score",
                         self.cfg.model_selection.upper())
                # Move to the next subset in the all_subsets list
                subset_index += 1

        log.info("Best scheme has %s score of %.2f and %d subset(s)"
                 % (self.cfg.model_selection.upper(), best_result.score,
                    len(best_scheme.subsets)))

        if fabricated_subsets:
            log.info("Finalising partitioning scheme")
            log.info("This involves cleaning up small subsets which %s "
                     "can't analyse", self.cfg.phylogeny_program)

        # Now join the fabricated subsets back up with other subsets
        while fabricated_subsets:
            log.info("***Kmeans algorithm step %d***" % step)
            step += 1

            # Take the first subset in the list (to be "popped" off later)
            s = fabricated_subsets[0]
            centroid = s.centroid

            best_match = None

            # Take a list copy of the best scheme
            scheme_list = list(best_scheme)
            scheme_list.remove(s)
            # Loop through the subsets in the best scheme and find the one
            # with the nearest centroid
            for sub in scheme_list:
                centroid_array = [sub.centroid, centroid]
                # euclid_dist = abs(sub.centroid[0] - centroid[0])
                warnings.simplefilter('ignore', DeprecationWarning)
                euclid_dist = spatial.distance.pdist(centroid_array)
                if euclid_dist < best_match or best_match is None:
                    best_match = euclid_dist
                    closest_sub = sub

            # Now merge those subsets
            merged_sub = subset_ops.merge_fabricated_subsets([s, closest_sub])

            # Remove the offending subset from the fabricated subset list
            fabricated_subsets.pop(0)
            # If the closest subset happens to be "fabricated" as well, take
            # it out of the fabricated_subsets list
            if closest_sub in fabricated_subsets:
                fabricated_subsets.remove(closest_sub)

            # Get rid of the two subsets that were merged from the best_scheme
            scheme_list.remove(closest_sub)

            # Now add the new subset to the scheme and see if the new subset
            # can be analyzed
            scheme_list.append(merged_sub)
            merged_scheme = scheme.Scheme(self.cfg, str(step - 1), scheme_list)

            merged_result = self.analyse_scheme(merged_scheme)
            # If it can be analyzed, move the algorithm forward, if it can't
            # be analyzed add it to the list of fabricated_subsets
            for new_subs in merged_scheme:
                if new_subs.fabricated and new_subs not in fabricated_subsets:
                    fabricated_subsets.append(new_subs)
            best_scheme = merged_scheme
            best_result = merged_result

        # Since the AIC will likely be better before we dealt with the
        # fabricated subsets, we need to set the best scheme and best result
        # to those from the last merged_scheme. TODO: add a variable to scheme
        # to take care of this problem so that the best AND analysable scheme
        # is the one that gets automatically flagged as the best scheme
        self.results.best_scheme = best_scheme
        self.results.best_result = best_result

        self.cfg.reporter.write_scheme_summary(
            self.results.best_scheme, self.results.best_result)

        log.info("** Kmeans algorithm finished after %d steps **" % (step - 1))
        log.info("Best scoring scheme is scheme %s, with %s score of %.3f"
                 % (self.results.best_scheme.name, self.cfg.model_selection,
                    self.results.best_score))

        self.cfg.reporter.write_best_scheme(self.results)


class KmeansAnalysisWrapper(Analysis):
    def do_analysis(self):
        # Copied and pasted from greedy analysis
        partnum = len(self.cfg.user_subsets)
        scheme_count = submodels.count_greedy_schemes(partnum)
        subset_count = submodels.count_greedy_subsets(partnum)

        self.cfg.progress.begin(scheme_count, subset_count)

        # Start with the most partitioned scheme
        start_description = range(partnum)
        start_scheme = scheme.create_scheme(
            self.cfg, "start_scheme", start_description)

        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        old_score = self.analyse_scheme(start_scheme)

        # Get first scheme
        best_scheme = start_scheme
        subset_index = 0

        processor = self.cfg.processor
        alignment_path = self.filtered_alignment_path
        tree_path = processor.make_tree_path(alignment_path)

        split_subsets = []
        for a_subset in start_scheme:
            how_many = kmeans.kmeans_wrapper(self.cfg, self.alignment,
                                             a_subset, tree_path)
            split_subsets += how_many
        split_scheme = scheme.Scheme(self.cfg, "split_scheme", split_subsets)
        best_result = self.analyse_scheme(best_scheme)
        split_score = self.analyse_scheme(split_scheme)
        if split_score.score < best_result.score:
            best_scheme = split_scheme
            log.info("Initial splits generated superior scheme")
        all_subsets = list(best_scheme.subsets)

        fabricated_subsets = []
        step = 1

        while subset_index < len(all_subsets):
            log.info("Best scheme has %s score of %.2f and %d subset(s)"
                     % (self.cfg.model_selection.upper(), best_result.score,
                        len(best_scheme.subsets)))

            log.info("***Kmeans algorithm step %d***" % step)
            step += 1

            current_subset = all_subsets[subset_index]

            log.info("Analysing subset of %d sites",
                     len(current_subset.columns))

            # First check if the subset is large enough to split, if it isn't,
            # move to the next subset
            if len(current_subset.columns) == 1:
                log.info("This subset cannot be split further")
                subset_index += 1
                continue

            if current_subset.fabricated:
                log.info(
                    "This subset cannot be split further because %s cannot "
                    "analyse it",
                    self.cfg.phylogeny_program)
                subset_index += 1
                fabricated_subsets.append(current_subset)
                continue

            split_subsets = kmeans.kmeans_split_subset(self.cfg,
                                                       self.alignment,
                                                       current_subset,
                                                       tree_path)


            # kmeans_split_subset will return a 1 and flag the subset as
            # fabricated if for some reason it raises a PhylogenyProgramError,
            # this it to catch those fabricated subsets
            if split_subsets == 1:
                subset_index += 1
                fabricated_subsets.append(current_subset)
                continue

            for each_subset in split_subsets:
                log.info("Subset resulting from split is %d sites long",
                         len(each_subset.columns))

            # Take a copy
            updated_subsets = all_subsets[:]

            # Replace the current one with the split one
            # Google "slice assignments"
            # This list is the key to avoiding recursion. It expands to contain
            # all of the split subsets by replacing them with the split ones
            updated_subsets[subset_index:subset_index + 1] = split_subsets

            test_scheme = scheme.Scheme(self.cfg, str(step - 1),
                                        updated_subsets)

            new_result = self.analyse_scheme(test_scheme)

            if new_result.score < best_result.score:
                best_scheme = test_scheme
                best_result = new_result

                # Change this to the one with split subsets in it. Note that
                # the subset_index now points a NEW subset, one that was split
                all_subsets = updated_subsets

                # record each scheme that's an improvement
                self.cfg.reporter.write_scheme_summary(
                    self.results.best_scheme, self.results.best_result)

                if len(split_subsets) == 2:
                    log.info(
                        "Splitting subset into %d:%d sites improved the %s score"
                        % (len(split_subsets[0].columns),
                           len(split_subsets[1].columns),
                           self.cfg.model_selection))

                    for s in split_subsets:
                        m = [x % 3 for x in s.columns]
                        l = float(len(s.columns))
                        props = [(float(m.count(1)) / l),
                                 (float(m.count(2)) / l),
                                 (float(m.count(0)) / l)]
                        log.info("%d subset has 1st, 2nd, 3rd props: %s" % (
                        len(s.columns), str(props)))



            else:
                log.info("Splitting this subset did not improve the %s score",
                         self.cfg.model_selection.upper())
                # Move to the next subset in the all_subsets list
                subset_index += 1

        log.info("Best scheme has %s score of %.2f and %d subset(s)"
                 % (self.cfg.model_selection.upper(), best_result.score,
                    len(best_scheme.subsets)))

        if fabricated_subsets:
            log.info("Finalising partitioning scheme")
            log.info("This involves cleaning up small subsets which %s "
                     "can't analyse", self.cfg.phylogeny_program)

        # Now join the fabricated subsets back up with other subsets
        while fabricated_subsets:
            log.info("***Kmeans algorithm step %d***" % step)
            step += 1

            # Take the first subset in the list (to be "popped" off later)
            s = fabricated_subsets[0]
            centroid = s.centroid

            best_match = None

            # Take a list copy of the best scheme
            scheme_list = list(best_scheme)
            scheme_list.remove(s)
            # Loop through the subsets in the best scheme and find the one
            # with the nearest centroid
            for sub in scheme_list:
                centroid_array = [sub.centroid, centroid]
                # euclid_dist = abs(sub.centroid[0] - centroid[0])
                warnings.simplefilter('ignore', DeprecationWarning)
                euclid_dist = spatial.distance.pdist(centroid_array)
                if euclid_dist < best_match or best_match is None:
                    best_match = euclid_dist
                    closest_sub = sub

            # Now merge those subsets
            merged_sub = subset_ops.merge_fabricated_subsets([s, closest_sub])

            # Remove the offending subset from the fabricated subset list
            fabricated_subsets.pop(0)
            # If the closest subset happens to be "fabricated" as well, take
            # it out of the fabricated_subsets list
            if closest_sub in fabricated_subsets:
                fabricated_subsets.remove(closest_sub)

            # Get rid of the two subsets that were merged from the best_scheme
            scheme_list.remove(closest_sub)

            # Now add the new subset to the scheme and see if the new subset
            # can be analyzed
            scheme_list.append(merged_sub)
            merged_scheme = scheme.Scheme(self.cfg, str(step - 1), scheme_list)

            merged_result = self.analyse_scheme(merged_scheme)
            # If it can be analyzed, move the algorithm forward, if it can't
            # be analyzed add it to the list of fabricated_subsets
            for new_subs in merged_scheme:
                if new_subs.fabricated and new_subs not in fabricated_subsets:
                    fabricated_subsets.append(new_subs)
            best_scheme = merged_scheme
            best_result = merged_result

        # Since the AIC will likely be better before we dealt with the
        # fabricated subsets, we need to set the best scheme and best result
        # to those from the last merged_scheme. TODO: add a variable to scheme
        # to take care of this problem so that the best AND analysable scheme
        # is the one that gets automatically flagged as the best scheme
        self.results.best_scheme = best_scheme
        self.results.best_result = best_result

        self.cfg.reporter.write_scheme_summary(
            self.results.best_scheme, self.results.best_result)

        log.info("** Kmeans algorithm finished after %d steps **" % (step - 1))
        log.info("Best scoring scheme is scheme %s, with %s score of %.3f"
                 % (self.results.best_scheme.name, self.cfg.model_selection,
                    self.results.best_score))

        self.cfg.reporter.write_best_scheme(self.results)


class KmeansGreedy(Analysis):
    def do_analysis(self):
        # Copied and pasted from greedy analysis
        partnum = len(self.cfg.user_subsets)
        scheme_count = submodels.count_greedy_schemes(partnum)
        subset_count = submodels.count_greedy_subsets(partnum)

        self.cfg.progress.begin(scheme_count, subset_count)

        # Start with the most partitioned scheme
        start_description = range(partnum)
        start_scheme = scheme.create_scheme(
            self.cfg, "start_scheme", start_description)

        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        old_score = self.analyse_scheme(start_scheme)

        # Get first scheme
        best_scheme = start_scheme
        subset_index = 0
        all_subsets = list(best_scheme.subsets)
        processor = self.cfg.processor
        alignment_path = self.filtered_alignment_path
        tree_path = processor.make_tree_path(alignment_path)

        while subset_index < len(all_subsets):
            current_subset = all_subsets[subset_index]
            split_subsets = kmeans.kmeans_split_subset(self.cfg, self.alignment,
                                                       current_subset,
                                                       tree_path)

            if split_subsets == 1:
                subset_index += 1

            else:
                # Take a copy
                updated_subsets = all_subsets[:]

                # Replace the current one with the split one
                # Google "slice assignments"
                # This list is the key to avoiding recursion. It expands to contain
                # all of the split subsets by replacing them with the split ones
                updated_subsets[subset_index:subset_index + 1] = split_subsets

                test_scheme = scheme.Scheme(self.cfg, "Current Scheme",
                                            updated_subsets)

                try:
                    best_result = self.analyse_scheme(best_scheme)
                    new_result = self.analyse_scheme(test_scheme)

                    log.info("Current best score is: " + str(best_result))
                    log.info("Current new score is: " + str(new_result))
                    if new_result.score < best_result.score:
                        log.info("New score " + str(
                            subset_index) + " is better and will be set to best score")
                        best_scheme = test_scheme

                        # Change this to the one with split subsets in it. Note that
                        # the subset_index now points a NEW subset, one that was split
                        all_subsets = updated_subsets
                    else:
                        # Move to the next subset in the all_subsets list
                        subset_index += 1

                # In PhyML or RAxML, it is likely because of no alignment patterns,
                # catch that and move to the next subset without splitting.
                except PhylogenyProgramError:
                    log.info(
                        "Phylogeny program generated an error so this subset was not split, see error above")
                    subset_index += 1
            # Now start the Greedy Analysis: need to figure out how to make it go through more
        # than one scheme...

        start_scheme = best_scheme
        partnum = len(start_scheme.subsets)
        scheme_count = submodels.count_greedy_schemes(partnum)
        subset_count = submodels.count_greedy_subsets(partnum)
        self.cfg.progress.begin(scheme_count, subset_count)
        start_description = range(partnum)

        step = 1
        cur_s = 2

        # Now we try out all lumpings of the current scheme, to see if we can
        # find a better one and if we do, we just keep going
        while True:
            log.info("***Greedy algorithm step %d***" % step)

            old_best_score = self.results.best_score

            # Get an iterable of all possible pairs of subsets in best_scheme
            lumped_subsets = itertools.combinations(start_scheme.subsets, 2)

            for subset_grouping in lumped_subsets:
                scheme_name = cur_s
                lumped_scheme = neighbour.make_clustered_scheme(
                    start_scheme, scheme_name, subset_grouping, self.cfg)

                new_result = self.analyse_scheme(lumped_scheme)

                log.debug("Difference in %s: %.1f",
                          self.cfg.model_selection,
                          (new_result.score - old_best_score))

                cur_s += 1

            if self.results.best_score != old_best_score:
                log.info("Analysed all schemes for this step. The best "
                         "scheme changed the %s score by %.1f units.",
                         self.cfg.model_selection,
                         (self.results.best_score - old_best_score))

                self.results.best_scheme.name = "step_%d" % step
                self.cfg.reporter.write_scheme_summary(
                    self.results.best_scheme, self.results.best_result)

                # Now we find out which is the best lumping we know of for this step
                start_scheme = self.results.best_scheme
            else:
                log.info(
                    "Analysed all schemes for this step and found no schemes "
                    "that improve the score, stopping")
                break

            # We're done if it's the scheme with everything together
            if len(set(lumped_scheme.subsets)) == 1:
                break

            step += 1

        log.info("Greedy algorithm finished after %d steps" % step)
        log.info("Best scoring scheme is scheme %s, with %s score of %.3f"
                 % (self.results.best_scheme.name, self.cfg.model_selection,
                    self.results.best_score))

        self.cfg.reporter.write_best_scheme(self.results)


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
    elif search == 'kmeans_wss':
        method = KmeansAnalysisWrapper
    elif search == 'kmeans_greedy':
        method = KmeansGreedy
    else:
        log.error("Search algorithm '%s' is not yet implemented", search)
        raise AnalysisError
    return method
