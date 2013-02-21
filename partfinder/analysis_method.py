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
log = logging.getLogger("method")

import os
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

        # Start with the most partitioned scheme
        start_description = range(len(self.cfg.partitions))
        start_scheme = scheme.create_scheme(self.cfg, "start_scheme", start_description)

        # Analyse our first scheme
        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        self.analyse_scheme(start_scheme)

        # Current scheme number
        cur_s = 2

        # Now we try out all clusterings of the first scheme, to see if we can
        # find a better one
        while True:
            log.info("***Clustering algorithm step %d of %d***" % (cur_s - 1, partnum - 1))

            # Calculate the subsets which are most similar
            # e.g. combined rank ordering of euclidean distances
            # Could combine average site-rates, q matrices, and frequencies
            scheme_name = "step_%d" % (cur_s - 1)
            clustered_scheme = neighbour.get_nearest_neighbour_scheme(
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

        partnum = len(self.cfg.partitions)
        scheme_count = submodels.count_greedy_schemes(partnum)
        subset_count = submodels.count_greedy_subsets(partnum)

        log.info("This will result in a maximum of %s schemes being created", scheme_count)
        log.info("PartitionFinder will have to analyse a maximum of %d subsets of sites to complete this analysis", subset_count)
        if subset_count > 10000:
            log.warning("%d is a lot of subsets, this might take a long time to analyse", subset_count)
            log.warning("Perhaps consider using a different search scheme instead (see Manual)")

        self.cfg.progress.begin(scheme_count, subset_count)

        # Start with the most partitioned scheme
        start_description = range(len(self.cfg.partitions))
        start_scheme = scheme.create_scheme(self.cfg, "start_scheme", start_description)

        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        result = self.analyse_scheme(start_scheme)

        step = 1
        cur_s = 2

        # Now we try out all lumpings of the current scheme, to see if we can
        # find a better one and if we do, we just keep going
        while True:
            scheme.tracker.create_snapshot()
            log.info("***Greedy algorithm step %d***" % step)

            # Get a list of all possible lumpings of the best_scheme
            lumpings = algorithm.lumpings(start_description)

            # Save the current best score we have in results
            old_best_score = self.results.best_score
            for lumped_description in lumpings:
                lumped_scheme = scheme.create_scheme(self.cfg, cur_s, lumped_description)
                cur_s += 1
                #this is just checking to see if a scheme is any good, if it is, we remember and write it later
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
            self.cfg.reporter.write_scheme_summary(self.results.best_scheme, self.results.best_result)

            # If it's the scheme with everything equal, quit
            if len(set(start_description)) == 1:
                break

            # Go do the next round...
            step += 1

        log.info("Greedy algorithm finished after %d steps" % step)
        log.info("Highest scoring scheme is scheme %s, with %s score of %.3f" %
                 (self.results.best_scheme.name, self.cfg.model_selection,
                  self.results.best_score))

    def report(self):
        txt = "Best scheme according to Greedy algorithm, analysed with %s" % self.cfg.model_selection
        self.cfg.reporter.write_best_scheme(txt, self.results)


class GreediestAnalysis(Analysis):
    def do_analysis(self):
        '''A greediest algorithm for heuristic partitioning searches
        1. Analyse up to greediest-percent of the possible lumpings
        2. If we find greediest-schemes of major improvements (delta score >10)
           then we take the best one*. If we find this many improvements before greediest-percent
           we quit early
        3. If we hit greediest-percent before we find greediest-schemes improvements, we
           take the best improvement*

        *best improvement includes all the imrpovements we find in a given scheme, plus all
        improvements we have found anywhere at all.

        4. If we get to greediest percent, and there are zero improvements we can make, then
           we quit.


        '''
        log.info("Performing greediest analysis")

        model_selection = self.cfg.model_selection
        partnum = len(self.cfg.partitions)

        scheme_count = submodels.count_greedy_schemes(partnum)
        subset_count = submodels.count_greedy_subsets(partnum)
        log.info("This will result in a maximum of %s schemes being created, probably a lot less", scheme_count)
        log.info("PartitionFinder will have to analyse a maximum of %d subsets of sites to complete this analysis, probably a lot less", subset_count)
        self.cfg.progress.begin(scheme_count, subset_count)

        #start with the most partitioned scheme
        start_description = range(len(self.cfg.partitions))
        start_scheme = scheme.create_scheme(
            self.cfg, "start_scheme", start_description)
        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        result = self.analyse_scheme(start_scheme)

        def get_score(my_result):
            try:
                return getattr(my_result, model_selection)
            except AttributeError:
                log.error("Unrecognised model_selection variable '%s', please check", model_selection)
                raise AnalysisError

        def process_best_scheme(all_improvements, best_score):
            #find the best score we can possibly find
            scores = list(all_improvements.keys())

            new_best_score = min(scores)
            best_improvements = (all_improvements[new_best_score])

            #in case there are >1, we just pick the first one
            r = best_improvements[0]
            best_result = r[0]
            best_lumping = r[1]
            best_name = r[2]
            best_lumped_scheme = r[3]

            write_this_scheme(name_prefix, best_result)
            log.info("Best scheme was number %s. It has %s: %.2f, %s difference: %.2f, and %d subsets." %(best_name.split('_')[-1], model_selection, new_best_score, model_selection, new_best_score-best_score, len(best_result.scheme.subsets)))

            return best_lumped_scheme, best_result

        def write_this_scheme(name_prefix, best_result):
            fname = os.path.join(self.cfg.schemes_path, name_prefix + '.txt')
            fobj = open(fname, 'wb')
            self.cfg.reporter.write_scheme_summary(best_result, fobj)
            self.results.add_scheme_result(best_result)

        best_score = get_score(result)
        best_result = result

        step = 1
        #now we try out all clusterings of the first scheme, to see if we can find a better one
        no_improvements=0
        subset_counter = 1 #start by remembering that we analysed the starting scheme
        while True:
            if no_improvements == 1:
                break
            log.info("***Greediest algorithm step %d of %d***" %
                     (step, partnum - 1))


            all_improvements = {}
            name_prefix ="step_%d" %(step)
            step += 1

            #get a list of all possible lumpings of the best_scheme, ordered according to the clustering weights
            lumped_subsets = neighbour.get_ranked_clustered_subsets(start_scheme, name_prefix, self.cfg)

            #set the scheme counter to run by algorithm step
            self.cfg.progress.schemes_analysed = 0
            self.cfg.progress.scheme_count = len(lumped_subsets)

            #now analyse the lumped schemes
            lumpings_done = 0
            for subset_grouping in lumped_subsets:

                scheme_name = "%s_%d" % (name_prefix, (lumpings_done + 1))
                lumped_scheme = neighbour.make_clustered_scheme(
                    start_scheme, scheme_name, subset_grouping, self.cfg)

                result = self.analyse_scheme(lumped_scheme, suppress_writing=True, suppress_memory=True)
                new_score = get_score(result)

                lumpings_done += 1

                if new_score < best_score:
                    #it's an improvement
                    log.info("Found improved score with delta %s: %.2f" %(model_selection, new_score-best_score))
                    #careful here - there might (might!) be >1 scheme with the same score...
                    results_with_this_score = all_improvements.setdefault(new_score, [])
                    results_with_this_score.append((result, subset_grouping, scheme_name, lumped_scheme))
                    all_improvements[new_score] = results_with_this_score

                #Stopping condition 1 - we found N improvements...
                #if len(all_improvements) >= int(self.cfg.greediest_schemes):
                #   log.info("Found %d schemes that improve the %s score, reached greediest-schemes "
                #             "cutoff condition." %(len(all_improvements), model_selection))
                #
                #   #now we find out which is the best lumping we know of for this step
                #   best_lumped_scheme, best_result = process_best_scheme(all_improvements, best_score)
                #   break

                #Stopping condition 2 - we got to greediest_percent way through...
                if (float(lumpings_done)/float(len(lumped_subsets)) >= self.cfg.greediest_percent*0.01):
                    #if we have any improvements, then we use the best one
                    if (len(all_improvements)>=1):
                        log.info("Analysed %.1f percent of the schemes for this step and found %d schemes "
                                 "that improve %s score, reached greediest-percent cutoff "
                                 "condition" %(self.cfg.greediest_percent, len(all_improvements), model_selection))

                        #now we find out which is the best lumping we know of for this step
                        best_lumped_scheme, best_result = process_best_scheme(all_improvements, best_score)
                        break

                    else: #we have no improvements
                        no_improvements=1
                        break

            #stop when we've anlaysed the scheme with all subsets combined, or if we can't find an improvement
            if len(set(lumped_scheme.subsets)) == 1:  # then it's the scheme with everything together
                break
            elif no_improvements==1:
                log.info("Analysed %.1f percent of the schemes for this step and found no schemes "
                         "that improve %s score, terminating greediest algorithm"
                         %(self.cfg.greediest_percent, model_selection))
                break
            else:
                if step==2:
                    subset_counter = subset_counter + submodels.a_choose_b(len(start_scheme.subsets), 2)
                else:
                    subset_counter = subset_counter + len(start_scheme.subsets)
                #before we loop, let's update the counter as if we'd analysed all the subsets
                self.cfg.progress.subsets_analysed = subset_counter
                start_scheme = best_lumped_scheme
                best_score = get_score(best_result)


        log.info("Greediest algorithm finished after %d steps" % step)
        log.info("Best scoring scheme is scheme %s, with %s score of %.3f"
                 % (best_result.scheme.name, model_selection, best_score))

        self.best_result = best_result

    def report(self):
        txt = "Best scheme according to Greediest algorithm, analysed with %s"
        best = [(txt % self.cfg.model_selection, self.best_result)]
        self.cfg.reporter.write_best_schemes(best)
        self.cfg.reporter.write_all_schemes(self.results, info= "Information on the best scheme from each step of the Greediest algorithm is here: ")


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
