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


import logging
log = logging.getLogger("analysis_method")

import os
import scheme
import algorithm
import submodels
import subset
from analysis import Analysis, AnalysisError
from neighbour import get_nearest_neighbour_scheme


class UserAnalysis(Analysis):

    def do_analysis(self):
        log.info("Performing User analysis")
        current_schemes = [s for s in self.cfg.schemes]
        scheme_count = len(current_schemes)
        subset_count = subset.count_subsets()

        self.cfg.progress.begin(scheme_count, subset_count)
        if scheme_count > 0:
            for s in current_schemes:
                self.analyse_scheme(s)
        else:
            log.error("Search set to 'user', but no user schemes detected in .cfg file. Please check.")
            raise AnalysisError

        self.cfg.progress.end()


class ClusteringAnalysis(Analysis):
    """This analysis uses model parameters to guess at similar partitions, then just joins them together
        this is much less accurate than other methods, but a LOT quicker - it runs in order
        N time (where N is the number of initial datablocks), whereas the greedy algorithm is
        still N squared.
    """

    def do_analysis(self):
        log.info("Performing clustering analysis")

        partnum = len(self.cfg.partitions)
        subset_count = 2 * partnum - 1
        scheme_count = partnum
        self.cfg.progress.begin(scheme_count, subset_count)

        # Start with the scheme with all subsets separate
        # Analyse that scheme

        # Clear any schemes that are currently loaded
        # TODO Not sure we need this...
        self.cfg.schemes.clear_schemes()

        # Start with the most partitioned scheme
        start_description = range(len(self.cfg.partitions))
        start_scheme = scheme.create_scheme(self.cfg, "start_scheme", start_description)
        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        self.analyse_scheme(start_scheme)

        cur_s = 2

        #now we try out all clusterings of the first scheme, to see if we can find a better one
        while True:
            log.info("***Clustering algorithm step %d of %d***" %
                     (cur_s - 1, partnum - 1))

            #calculate the subsets which are most similar
            #e.g. combined rank ordering of euclidean distances
            #could combine average site-rates, q matrices, and frequencies
            scheme_name ="step_%d" %(cur_s-1)
            clustered_scheme = get_nearest_neighbour_scheme(
                start_scheme, scheme_name, self.cfg)

            #now analyse that new scheme
            cur_s += 1
            self.analyse_scheme(clustered_scheme)

            #stop when we've anlaysed the scheme with all subsets combined
            if len(set(clustered_scheme.subsets)) == 1:  # then it's the scheme with everything together
                break
            else:
                start_scheme = clustered_scheme

        self.cfg.progress.end()


class AllAnalysis(Analysis):

    def do_analysis(self):
        log.info("Performing complete analysis")
        partnum = len(self.cfg.partitions)

        scheme_count = submodels.count_all_schemes(partnum)
        subset_count = submodels.count_all_subsets(partnum)
        log.info("Analysing all possible schemes for %d starting partitions",
                 partnum)
        log.info("This will result in %s schemes being created",
                 scheme_count)
        self.cfg.progress.begin(scheme_count, subset_count)

        log.info("PartitionFinder will have to analyse %d subsets to complete this analysis", subset_count)
        if subset_count > 10000:
            log.warning("%d is a lot of subsets, this might take a long time to analyse", subset_count)
            log.warning("Perhaps consider using a different search scheme instead (see Manual)")

        # Clear any schemes that are currently loaded
        self.cfg.schemes.clear_schemes()

        # Iterate over submodels, which we can turn into schemes afterwards in the loop
        model_iterator = submodels.submodel_iterator([], 1, partnum)

        scheme_name = 1
        for m in model_iterator:
            s = scheme.model_to_scheme(m, scheme_name, self.cfg)
            scheme_name = scheme_name + 1
            self.analyse_scheme(s)


class GreedyAnalysis(Analysis):
    

    def do_analysis(self):
        '''A greedy algorithm for heuristic partitioning searches'''
        log.info("Performing greedy analysis")
        model_selection = self.cfg.model_selection
        partnum = len(self.cfg.partitions)

        scheme_count = submodels.count_greedy_schemes(partnum)
        subset_count = submodels.count_greedy_subsets(partnum)
        log.info("This will result in a maximum of %s schemes being created", scheme_count)
        log.info("PartitionFinder will have to analyse a maximum of %d subsets of sites to complete this analysis", subset_count)
        self.cfg.progress.begin(scheme_count, subset_count)

        if subset_count > 10000:
            log.warning("%d is a lot of subsets, this might take a long time to analyse", subset_count)
            log.warning("Perhaps consider using a different search scheme instead (see Manual)")

        #clear any schemes that are currently loaded
        # TODO Not sure we need this...
        self.cfg.schemes.clear_schemes()

        #start with the most partitioned scheme
        start_description = range(len(self.cfg.partitions))
        start_scheme = scheme.create_scheme(self.cfg, "start_scheme"    , start_description)
        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        result = self.analyse_scheme(start_scheme)

        def get_score(my_result):
            try:
                return getattr(my_result, model_selection)
            except AttributeError:
                log.error("Unrecognised model_selection variable '%s', please check", model_selection)
                raise AnalysisError

        best_result = result
        best_score = get_score(result)

        step = 1
        cur_s = 2

        #now we try out all lumpings of the current scheme, to see if we can find a better one
        #and if we do, we just keep going
        while True:
            log.info("***Greedy algorithm step %d***" % step)

            #get a list of all possible lumpings of the best_scheme
            lumpings = algorithm.lumpings(start_description)

            #we reset the counters as we go, for better user information
            self.total_scheme_num = len(lumpings)
            self.schemes_analysed = 0

            best_lumping_score = None
            for lumped_description in lumpings:
                lumped_scheme = scheme.create_scheme(
                    self.cfg, cur_s, lumped_description)
                cur_s += 1
                #this is just checking to see if a scheme is any good, if it is, we remember and write it later
                result = self.analyse_scheme(lumped_scheme, suppress_writing=True, suppress_memory=True)
                new_score = get_score(result)

                if best_lumping_score is None or new_score < best_lumping_score:
                    best_lumping_score = new_score
                    best_lumping_result = result
                    best_lumping_desc = lumped_description
                    best_scheme = lumped_scheme

            if best_lumping_score < best_score:
                best_score = best_lumping_score
                start_description = best_lumping_desc
                best_scheme = lumped_scheme
                best_result = best_lumping_result
                #now we write out the result of the best scheme for each step...
                fname = os.path.join(self.cfg.schemes_path, "step_%d" %step + '.txt')
                self.cfg.reporter.write_scheme_summary(best_result, open(fname, 'w'))
                self.results.add_scheme_result(best_result)

                if len(set(best_lumping_desc)) == 1:  # then it's the scheme with everything equal, so quit
                    break
                step += 1

            else:
                break

        log.info("Greedy algorithm finished after %d steps" % step)
        log.info("Highest scoring scheme is scheme %s, with %s score of %.3f"
                 % (best_result.scheme.name, model_selection, best_score))

        self.best_result = best_result

    def report(self):
        txt = "Best scheme according to Greedy algorithm, analysed with %s"
        best = [(txt % self.cfg.model_selection, self.best_result)]
        self.cfg.reporter.write_best_schemes(best)
        self.cfg.reporter.write_all_schemes(self.results, info= "Information on the best scheme from each step of the greedy algorithm is here: ")


def choose_method(search):
    if search == 'all':
        method = AllAnalysis
    elif search == 'user':
        method = UserAnalysis
    elif search == 'greedy':
        method = GreedyAnalysis
    elif search == 'clustering':
        method = ClusteringAnalysis
    else:
        log.error("Search algorithm '%s' is not yet implemented", search)
        raise AnalysisError
    return method
