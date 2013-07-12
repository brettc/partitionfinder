# Copyright (C) 2012 Robert Lanfear and Brett Calcott
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


import logging
log = logging.getLogger("method")

import os
import math
import scheme
import algorithm
import submodels
import subset
from analysis import Analysis, AnalysisError
import neighbour
import kmeans
import raxml
import phyml
import subset_ops

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
            log.error("Search set to 'user', but no user schemes detected in .cfg file. Please check.")
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

            # Stop when we've anlaysed the scheme with all subsets combined
            if len(set(clustered_scheme.subsets)) == 1:  # then it's the scheme with everything together
                break
            else:
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
    def do_analysis(self):
        '''A greedy algorithm for heuristic partitioning searches'''

        log.info("Performing greedy analysis")

        partnum = len(self.cfg.user_subsets)
        scheme_count = submodels.count_greedy_schemes(partnum)
        subset_count = submodels.count_greedy_subsets(partnum)

        self.cfg.progress.begin(scheme_count, subset_count)

        # Start with the most partitioned scheme
        start_description = range(partnum)
        start_scheme = scheme.create_scheme(
            self.cfg, "start_scheme", start_description)

        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        self.analyse_scheme(start_scheme)

        step = 1
        cur_s = 2

        # Now we try out all lumpings of the current scheme, to see if we can
        # find a better one and if we do, we just keep going
        while True:
            log.info("***Greedy algorithm step %d***" % step)

            # Get a list of all possible lumpings of the best_scheme
            lumpings = algorithm.lumpings(start_description)

            # Save the current best score we have in results
            old_best_score = self.results.best_score
            for lumped_description in lumpings:
                lumped_scheme = scheme.create_scheme(self.cfg, cur_s, lumped_description)
                cur_s += 1
                # This is just checking to see if a scheme is any good, if it
                # is, we remember and write it later
                self.analyse_scheme(lumped_scheme)

            # Did out best score change (It ONLY gets better -- see in
            # results.py)
            if self.results.best_score == old_best_score:
                # It didn't, so we're done
                break

            # Let's look further. We use the description from our best scheme
            # (which will be the one that just changed in the last lumpings
            # iteration)
            start_description = self.results.best_result.scheme.description

            # Rename and record the best scheme for this step
            self.results.best_scheme.name = "step_%d" % step
            self.cfg.reporter.write_scheme_summary(
                self.results.best_scheme, self.results.best_result)

            # If it's the scheme with everything equal, quit
            if len(set(start_description)) == 1:
                break

            # Go do the next round...
            step += 1

        log.info("Greedy algorithm finished after %d steps" % step)
        log.info("Highest scoring scheme is scheme %s, with %s score of %.3f" %
                 (self.results.best_scheme.name, self.cfg.model_selection,
                  self.results.best_score))

        self.cfg.reporter.write_best_scheme(self.results)


class RelaxedClusteringAnalysis(Analysis):
    '''
    A relaxed clustering algorithm for heuristic partitioning searches

    1. Rank subsets by their similarity (defined by clustering-weights)
    2. Analyse cluster-percent of the most similar schemes
    3. Take the scheme that improves the AIC/BIC score the most
    4. Quit if no improvements.
    '''

    def do_analysis(self):
        log.info("Performing relaxed clustering analysis")

        stop_at = self.cfg.cluster_percent * 0.01

        model_selection = self.cfg.model_selection
        partnum = len(self.cfg.user_subsets)

        scheme_count = submodels.count_relaxed_clustering_schemes(
            partnum, self.cfg.cluster_percent)
        subset_count = submodels.count_relaxed_clustering_subsets(
            partnum, self.cfg.cluster_percent)

        self.cfg.progress.begin(scheme_count, subset_count)

        # Start with the most partitioned scheme, and record it.
        start_description = range(partnum)
        start_scheme = scheme.create_scheme(
            self.cfg, "start_scheme", start_description)
        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        self.analyse_scheme(start_scheme)
        self.cfg.reporter.write_scheme_summary(
            self.results.best_scheme, self.results.best_result)

        # Start by remembering that we analysed the starting scheme
        subset_counter = 1
        step = 1
        while True:

            log.info("***Relaxed clustering algorithm step %d of %d***" % (step, partnum - 1))
            name_prefix = "step_%d" % (step)

            # Get a list of all possible lumpings of the best_scheme, ordered
            # according to the clustering weights
            lumped_subsets = neighbour.get_ranked_clustered_subsets(
                start_scheme, self.cfg)

            # Reduce the size of the lumped subsets to cluster_percent long
            # Round up to stop zeros
            cutoff = int(math.ceil(len(lumped_subsets)*stop_at))
            lumped_subsets = lumped_subsets[:cutoff]

            # Now analyse the lumped schemes
            lumpings_done = 0
            old_best_score = self.results.best_score

            for subset_grouping in lumped_subsets:
                scheme_name = "%s_%d" % (name_prefix, lumpings_done + 1)
                lumped_scheme = neighbour.make_clustered_scheme(
                    start_scheme, scheme_name, subset_grouping, self.cfg)

                new_result = self.analyse_scheme(lumped_scheme)

                log.debug("Difference in %s: %.1f",
                          self.cfg.model_selection,
                          (new_result.score-old_best_score))

                lumpings_done += 1

            if self.results.best_score != old_best_score:
                log.info("Analysed %.1f percent of the schemes for this step. The best "
                         "scheme changed the %s score by %.1f units.",
                         self.cfg.cluster_percent, self.cfg.model_selection,
                         (self.results.best_score - old_best_score))

                self.results.best_scheme.name = "step_%d" % step
                self.cfg.reporter.write_scheme_summary(
                    self.results.best_scheme, self.results.best_result)

                # Now we find out which is the best lumping we know of for this step
                start_scheme = self.results.best_scheme
            else:
                log.info("Analysed %.1f percent of the schemes for this step and found no schemes "
                         "that improve the score, stopping" , self.cfg.cluster_percent)
                break

            # We're done if it's the scheme with everything together
            if len(set(lumped_scheme.subsets)) == 1:
                break

            step += 1

        log.info("Relaxed clustering algorithm finished after %d steps" % step)
        log.info("Best scoring scheme is scheme %s, with %s score of %.3f"
                 % (self.results.best_scheme.name, model_selection, self.results.best_score))

        self.cfg.reporter.write_best_scheme(self.results)

class KmeansAnalysis(Analysis):
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


        while subset_index < len(all_subsets):
            current_subset = all_subsets[subset_index]
            split_subsets = kmeans.kmeans_split_subset(self.cfg, self.alignment, current_subset)

            if split_subsets == 1:
                subset_index += 1

            else:
                # Take a copy
                updated_subsets = all_subsets[:]

                # Replace the current one with the split one
                # Google "slice assignments"
                # This list is the key to avoiding recursion. It expands to contain
                # all of the split subsets by replacing them with the split ones
                updated_subsets[subset_index:subset_index+1] = split_subsets

                test_scheme = scheme.Scheme(self.cfg, "bla", updated_subsets)

                try:
                    best_score = self.analyse_scheme(best_scheme)
                    new_score = self.analyse_scheme(test_scheme)

                    log.info("Current best score is: " + str(best_score))
                    log.info("Current new score is: " + str(new_score))
                    if new_score.score < best_score.score:
                        log.info("New score " + str(subset_index) + " is better and will be set to best score")
                        best_scheme = test_scheme

                        # Change this to the one with split subsets in it. Note that
                        # the subset_index now points a NEW subset, one that was split
                        all_subsets = updated_subsets
                    else:
                        # Move to the next subset in the all_subsets list
                        subset_index += 1

                # In PhyML or RAxML, it is likely because of no alignment patterns,
                # catch that and move to the next subset without splitting.
                except PhylogenyProgramError as e:
                    log.info("Bummer: %s" % e)
                    subset_index += 1
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

        split_subsets = []
        for a_subset in start_scheme:
            how_many = kmeans.kmeans_wrapper(self.cfg, self.alignment, a_subset)
            split_subsets += how_many
        split_scheme = scheme.Scheme(self.cfg, "split_scheme", split_subsets)
        best_score = self.analyse_scheme(best_scheme)
        split_score = self.analyse_scheme(split_scheme)
        if split_score.score < best_score.score:
            best_scheme = split_scheme
            log.info("Initial splits generated superior scheme")
        all_subsets = list(best_scheme.subsets)

        while subset_index < len(all_subsets):
            current_subset = all_subsets[subset_index]
            split_subsets = kmeans.kmeans_split_subset(self.cfg, self.alignment, current_subset)

            if split_subsets == 1:
                subset_index += 1

            else:
                # Take a copy
                updated_subsets = all_subsets[:]

                # Replace the current one with the split one
                # Google "slice assignments"
                # This list is the key to avoiding recursion. It expands to contain
                # all of the split subsets by replacing them with the split ones
                updated_subsets[subset_index:subset_index+1] = split_subsets

                test_scheme = scheme.Scheme(self.cfg, "bla", updated_subsets)

                try:
                    best_score = self.analyse_scheme(best_scheme)
                    new_score = self.analyse_scheme(test_scheme)

                    log.info("Current best score is: " + str(best_score))
                    log.info("Current new score is: " + str(new_score))
                    if new_score.score < best_score.score:
                        log.info("New score " + str(subset_index) + " is better and will be set to best score")
                        best_scheme = test_scheme

                        # Change this to the one with split subsets in it. Note that
                        # the subset_index now points a NEW subset, one that was split
                        all_subsets = updated_subsets
                    else:
                        # Move to the next subset in the all_subsets list
                        subset_index += 1

                # RAxML and PhyML will choke on partitions that have all the
                # same alignment patterns. This will move the analysis along
                # without splitting that subset if that happens.
                except PhylogenyProgramError, e:
                    log.info("Bummer: %s" % e)
                    subset_index += 1
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
    elif search == 'paul':
        method = KmeansAnalysis
    elif search == 'paul2':
        method = KmeansAnalysisWrapper
    else:
        log.error("Search algorithm '%s' is not yet implemented", search)
        raise AnalysisError
    return method
