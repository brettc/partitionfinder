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

class UserAnalysis(Analysis):

    def do_analysis(self):
        log.info("Performing User analysis")
        current_schemes = [s for s in self.cfg.user_schemes]
        scheme_count = len(current_schemes)
        subset_count = subset.count_subsets()

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

        txt = "Best scheme out of all User Schemes, analysed with %s" % self.cfg.model_selection
        self.cfg.reporter.write_best_scheme(txt, self.results)


class ClusteringAnalysis(Analysis):
    """
    This analysis uses model parameters to guess at similar partitions, then
    just joins them together this is much less accurate than other methods, but
    a LOT quicker - it runs in order N time (where N is the number of initial
    datablocks), whereas the greedy algorithm is still N squared.
    """

    def do_analysis(self):
        log.info("Performing clustering analysis")

        partnum = len(self.cfg.partitions)
        subset_count = 2 * partnum - 1
        scheme_count = partnum
        self.cfg.progress.begin(scheme_count, subset_count)

        # Start with the most partitioned scheme
        start_description = range(len(self.cfg.partitions))
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
            log.info("***Clustering algorithm step %d of %d***" %
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

        txt = "Best scheme using Clustering Analysis, analysed with %s" % self.cfg.model_selection
        self.cfg.reporter.write_best_scheme(txt, self.results)


class AllAnalysis(Analysis):

    def do_analysis(self):
        log.info("Performing complete analysis")
        partnum = len(self.cfg.partitions)

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

        txt = "Best scheme out of all Schemes, analysed with %s" % self.cfg.model_selection
        self.cfg.reporter.write_best_scheme(txt, self.results)


class GreedyAnalysis(Analysis):
    def do_analysis(self):
        '''A greedy algorithm for heuristic partitioning searches'''

        log.info("Performing greedy analysis")

        partnum = len(self.cfg.partitions)
        scheme_count = submodels.count_greedy_schemes(partnum)
        subset_count = submodels.count_greedy_subsets(partnum)

        self.cfg.progress.begin(scheme_count, subset_count)

        # Start with the most partitioned scheme
        start_description = range(len(self.cfg.partitions))
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

        txt = "Best scheme according to Greedy algorithm, analysed with %s" % self.cfg.model_selection
        self.cfg.reporter.write_best_scheme(txt, self.results)


class GreediestAnalysis(Analysis):
    '''
    A greediest algorithm for heuristic partitioning searches

    1. Analyse up to greediest-percent of the possible lumpings
    2. If we find greediest-schemes of major improvements (delta score >10)
        then we take the best one*. If we find this many improvements before
        greediest-percent we quit early
    3. If we hit greediest-percent before we find greediest-schemes
        improvements, we take the best improvement*

    *best improvement includes all the imrpovements we find in a given
    scheme, plus all improvements we have found anywhere at all.

    4. If we get to greediest percent, and there are zero improvements we
        can make, then we quit.
    '''

    def do_analysis(self):
        log.info("Performing greediest analysis")

        stop_at = self.cfg.greediest_percent * 0.01

        model_selection = self.cfg.model_selection
        partnum = len(self.cfg.partitions)

        scheme_count = int(math.ceil(submodels.count_greedy_schemes(partnum)*stop_at))
        subset_count = int(math.ceil(submodels.count_greedy_subsets(partnum)*stop_at))

        self.cfg.progress.begin(scheme_count, subset_count)

        # Start with the most partitioned scheme, and record it.
        start_description = range(len(self.cfg.partitions))
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

            log.info("***Greediest algorithm step %d of %d***" % (step, partnum - 1))
            name_prefix = "step_%d" % (step)

            # Get a list of all possible lumpings of the best_scheme, ordered
            # according to the clustering weights
            lumped_subsets = neighbour.get_ranked_clustered_subsets(
                start_scheme, self.cfg)

            # reduce the size of the lumped subsets to greediest_percent long
            cutoff = int(math.ceil(len(lumped_subsets)*stop_at)) #round up to stop zeros            
            lumped_subsets = lumped_subsets[:cutoff]

            # Now analyse the lumped schemes
            lumpings_done = 0
            old_best_score = self.results.best_score

            for subset_grouping in lumped_subsets:
                scheme_name = "%s_%d" % (name_prefix, lumpings_done + 1)
                lumped_scheme = neighbour.make_clustered_scheme(
                    start_scheme, scheme_name, subset_grouping, self.cfg)

                new_result = self.analyse_scheme(lumped_scheme)
                
                log.info("Difference in %s: %.1f", self.cfg.model_selection, (new_result.score-old_best_score))
                
                lumpings_done += 1

            
            if self.results.best_score != old_best_score:
                log.info("Analysed %.1f percent of the schemes for this step. The best "
                         "scheme changed the %s score by %.1f units.", 
                         self.cfg.greediest_percent, self.cfg.model_selection,
                         (self.results.best_score - old_best_score))

                #write out the best scheme
                self.results.best_scheme.name = "step_%d" % step
                self.cfg.reporter.write_scheme_summary(
                    self.results.best_scheme, self.results.best_result)


                # Now we find out which is the best lumping we know of for this step
                start_scheme = self.results.best_scheme
            else:
                log.info("Analysed %.1f percent of the schemes for this step and found no schemes "
                         "that improve the score, stopping" , self.cfg.greediest_percent)
                break

            # We're done if it's the scheme with everything together
            if len(set(lumped_scheme.subsets)) == 1:
                break
    
            step += 1



        log.info("Greediest algorithm finished after %d steps" % step)
        log.info("Best scoring scheme is scheme %s, with %s score of %.3f"
                 % (self.results.best_scheme.name, model_selection, self.results.best_score))

        txt = "Best scheme according to Greediest algorithm, analysed with %s" % self.cfg.model_selection
        self.cfg.reporter.write_best_scheme(txt, self.results)


def choose_method(search):
    if search == 'all':
        method = AllAnalysis
    elif search == 'user':
        method = UserAnalysis
    elif search == 'greedy':
        method = GreedyAnalysis
    elif search == 'clustering':
        method = ClusteringAnalysis
    elif search == 'greediest':
        method = GreediestAnalysis
    else:
        log.error("Search algorithm '%s' is not yet implemented", search)
        raise AnalysisError
    return method
