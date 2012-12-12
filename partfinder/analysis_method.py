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
import submodels

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

        percentage_cutoff = float(20)

        model_selection = self.cfg.model_selection
        partnum = len(self.cfg.partitions)

        scheme_count = submodels.count_greedy_schemes(partnum)
        subset_count = submodels.count_greedy_subsets(partnum)
        log.info("This will result in a maximum of %s schemes being created, probably a lot less", scheme_count)
        log.info("PartitionFinder will have to analyse a maximum of %d subsets of sites to complete this analysis, probably a lot less", subset_count)
        self.cfg.progress.begin(scheme_count, subset_count)

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

        def add_result(all_results, this_result, this_name):
            this_score = get_score(this_result)
            results_with_this_score = all_improvements.setdefault(this_score, [])
            results_with_this_score.append((this_result, subset_grouping, this_name))
            all_improvements[this_score] = results_with_this_score

        def process_best_scheme(all_improvements, lumped_subsets):
            #find the best score we can possibly find
            scores = list(all_improvements.keys())
            scores.sort()

            #now we iterate through the scores, and pick the best one that corresponds to any of our subset_groupings
            for score in scores:
                for r in all_improvements[score]: #it's a list of potentially >1 (result, subset_grouping) key
                    for s in lumped_subsets:
                        if r[1]==s: #we have a subset grouping that matches...
                            best_result = r[0]
                            best_lumping = r[1]
                            best_name = r[2]

                            #now we remove that lumping from the all_improvements dictionary
                            all_lumps = all_improvements[score]
                            all_lumps.remove(r)
                            if len(all_lumps)==0:
                                del all_improvements[score]
                            else:
                                all_improvements[score] = all_lumps
                            return all_improvements, best_lumping, best_name


            #if we got to here, then there's no best result
            return all_improvements, None, None

        def write_this_scheme(name_prefix, best_result):
            fname = os.path.join(self.cfg.schemes_path, name_prefix + '.txt')
            fobj = open(fname, 'w')
            self.cfg.reporter.write_scheme_summary(best_result, fobj)
            self.results.add_scheme_result(best_result)
            #before we break, let's update the counters as if we'd analysed all those schemes and subsets
            self.cfg.progress.subsets_analysed = self.cfg.progress.subsets_analysed + len(start_scheme.subsets)


        best_score = get_score(result)

        step = 1
        #now we try out all clusterings of the first scheme, to see if we can find a better one
        no_improvements=0
        while True:
            if no_improvements==1:
                break
            log.info("***Greediest algorithm step %d of %d***" %
                     (step, partnum - 1))

            all_improvements = {}
            schemes_this_step = {}
            name_prefix ="step_%d" %(step)
            step += 1

            #get a list of all possible lumpings of the best_scheme
            lumped_subsets = neighbour.get_ranked_clustered_subsets(start_scheme, name_prefix, self.cfg)

            #set the scheme counter to run by algorithm step
            self.cfg.progress.schemes_analysed = 0
            self.cfg.progress.scheme_count = len(lumped_subsets)

            #now analyse the lumped schemes
            lumpings_done = 0
            major_improvements = 0
            for subset_grouping in lumped_subsets:

                scheme_name = "%s_%d" %(name_prefix, (lumpings_done+1))
                lumped_scheme = neighbour.make_clustered_scheme(start_scheme, scheme_name, subset_grouping, self.cfg)

                result = self.analyse_scheme(lumped_scheme, suppress_writing=True, suppress_memory=True)
                new_score = get_score(result)

                schemes_this_step[scheme_name] = (lumped_scheme, result)
                lumpings_done += 1

                if new_score<best_score:
                    #it's an improvement
                    log.info("Found improved score with delta %s: %.2f" %(model_selection, new_score-best_score))
                    add_result(all_improvements, result, scheme_name) #add it to the all_improvements dictionary
                    if new_score<(best_score-10):
                        major_improvements += 1

                #Stopping condition 1 - we found N major improvements...
                if major_improvements >= int(self.cfg.greediest_schemes):
                    log.info("Found %d schemes that improve the %s score by >10 units, reached greediest-schemes "
                             "cutoff condition." %(major_improvements, model_selection))

                    #now we find out which is the best lumping we know of for this step
                    all_improvements, best_lumping, best_name = process_best_scheme(all_improvements, lumped_subsets)
                    print "best name: ", best_name

                    if best_lumping==None:
                        #we SHOULD be able to find one, since we've had major improvements!
                        log.error("Something has gone wrong with the greediest algorithm")
                        raise AnalysisError

                    #now we check if that scheme already exists
                    if best_name in schemes_this_step:
                        print "EGGER"
                        best_lumped_scheme = schemes_this_step[best_name][0]
                        best_result = schemes_this_step[best_name][1]
                        write_this_scheme(name_prefix, best_result)
                    else:
                        log.error("Something has gone wrong with the greediest algorithm")
                        raise AnalysisError

                    log.info("Best scheme has %s difference: %.2f, and %d subsets" %(model_selection, get_score(best_result)-best_score, len(best_result.scheme.subsets)))
                    #now we move on
                    break

                #Stopping condition 2 - we got to greediest_percent way through...
                if (float(lumpings_done)/float(len(lumped_subsets)) >= self.cfg.greediest_percent*0.01) and (lumpings_done>self.cfg.greediest_schemes):
                    #if we have any improvements, then we use the best one
                    if (len(all_improvements)>=1):
                        log.info("Analysed %.1f percent of the schemes for this step and found only %d schemes "
                                 "that improve %s score by >10 units, reached greediest-percent cutoff "
                                 "condition" %(self.cfg.greediest_percent, major_improvements, model_selection))

                        #now we find out which is the best lumping we know of for this step
                        all_improvements, best_lumping, best_name = process_best_scheme(all_improvements, lumped_subsets)
                        print "best name: ", best_name

                        if best_lumping==None:
                            #we SHOULD be able to find one, since we've had >1 improvement!
                            log.error("Something has gone wrong with the greediest algorithm")
                            raise AnalysisError
                        if best_name in schemes_this_step:
                            print "EGGER2"
                            best_lumped_scheme = schemes_this_step[best_name][0]
                            best_result = schemes_this_step[best_name][1]
                            write_this_scheme(name_prefix.split('_')[0], best_result)
                        else:
                            log.error("Something has gone wrong with the greediest algorithm")
                            raise AnalysisError

                        log.info("Best scheme has %s difference: %.2f, and %d subsets" %(model_selection, get_score(best_result)-best_score, len(best_result.scheme.subsets)))
                        #now we move on
                        break
                    else: #we have no improvements
                        no_improvements=1
                        break


            #stop when we've anlaysed the scheme with all subsets combined
            if len(set(lumped_scheme.subsets)) == 1:  # then it's the scheme with everything together
                break
            elif no_improvements==1:
                log.info("Analysed %.1f percent of the schemes for this step and found no schemes "
                         "that improve %s score by any amount, terminating greediest algorithm"
                         %(self.cfg.greediest_percent, major_improvements, model_selection))
                break
            else:
                start_scheme = best_lumped_scheme
                best_score = get_score(best_result)

        log.info("Greediest algorithm finished after %d steps" % step)
        log.info("Best scoring scheme is scheme %s, with %s score of %.3f"
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
    elif search == 'greediest':
        method = GreediestAnalysis
    else:
        log.error("Search algorithm '%s' is not yet implemented", search)
        raise AnalysisError
    return method
