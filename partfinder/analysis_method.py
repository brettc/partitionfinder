import logging
log = logging.getLogger("analysis_method")

import os
import scheme
import algorithm
import submodels
from analysis import Analysis, AnalysisError

def get_score(my_result, model_selection):
    #TODO: this is bad. Should use self.cfg.model_selection, or write
    #a new model_selection for scheme.py
    if model_selection=="aic":
        score=my_result.aic
    elif model_selection=="aicc":
        score=my_result.aicc
    elif model_selection=="bic":
        score=my_result.bic
    else:
        log.error("Unrecognised model_selection variable '%s', please check" %(score))
        raise AnalysisError
    return score


class UserAnalysis(Analysis):

    def do_analysis(self):
        """Process everything when search=user"""
        models = self.cfg.models

        current_schemes = [s for s in self.cfg.schemes]
        self.total_scheme_num = len(current_schemes)
        if self.total_scheme_num>0:
            for s in current_schemes:
                 self.analyse_scheme(s, models)
        else:
            log.error("Search set to 'user', but no user schemes detected in .cfg file. Please check.")
            raise AnalysisError

class ClusteringAnalysis(Analysis):
    """This analysis uses model parameters to guess at similar partitions, then just joins them together
        this is much less accurate than other methods, but a LOT quicker - it runs in order
        N time (where N is the number of initial datablocks), whereas the greedy algorithm is
        still N squared.
    """
    
    def do_analysis(self):
        log.info("Performing clustering analysis")
    
        models = self.cfg.models
        model_selection = self.cfg.model_selection
        partnum = len(self.cfg.partitions)
        self.total_subset_num = 2*partnum-1
        self.total_scheme_num = partnum

        #start with the scheme with all subsets separate
        #analyse that scheme
        #clear any schemes that are currently loaded
        # TODO Not sure we need this...
        self.cfg.schemes.clear_schemes()        
                
        #start with the most partitioned scheme
        start_description = range(len(self.cfg.partitions))
        start_scheme = scheme.create_scheme(self.cfg, 1, start_description)
        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        result = self.analyse_scheme(start_scheme, models)
        
        best_result = result
        best_score  = get_score(result, self.cfg.model_selection)
                         
        cur_s = 2

        #now we try out all clusterings of the first scheme, to see if we can find a better one
        while True:
            log.info("***Clustering algorithm step %d of %d***" %(cur_s-1, partnum-1))
                                        
            #calculate the subsets which are most similar
            #e.g. combined rank ordering of euclidean distances
            #could combine average site-rates, q matrices, and frequencies
            clustered_scheme = start_scheme.get_clustering(self.cfg, method='hierarchical', scheme_name = cur_s)
            
            #now analyse that new scheme
            cur_s += 1
            result = self.analyse_scheme(clustered_scheme, models)
            
            #stop when we've anlaysed the scheme with all subsets combined
            if len(set(clustered_scheme.subsets)) == 1: #then it's the scheme with everything together
                break
            else:
                start_scheme = clustered_scheme
                                
                
class AllAnalysis(Analysis):

    def do_analysis(self):
        models = self.cfg.models
        partnum = len(self.cfg.partitions)

        self.total_scheme_num = submodels.count_all_schemes(partnum)
        log.info("Analysing all possible schemes for %d starting partitions", partnum)
        log.info("This will result in %s schemes being created", self.total_scheme_num)
        self.total_subset_num = submodels.count_all_subsets(partnum)
        log.info("PartitionFinder will have to analyse %d subsets to complete this analysis" %(self.total_subset_num))
        if self.total_subset_num>10000:
            log.warning("%d is a lot of subsets, this might take a long time to analyse", self.total_subset_num)
            log.warning("Perhaps consider using a different search scheme instead (see Manual)")

        #clear any schemes that are currently loaded
        self.cfg.schemes.clear_schemes()

        #iterate over submodels, which we can turn into schemes afterwards in the loop
        model_iterator = submodels.submodel_iterator([], 1, partnum)

        scheme_name = 1
        for m in model_iterator:
            s = scheme.model_to_scheme(m, scheme_name, self.cfg)
            scheme_name = scheme_name+1
            self.analyse_scheme(s, models)

class GreedyAnalysis(Analysis):

    def do_analysis(self):
        '''A greedy algorithm for heuristic partitioning searches'''
        log.info("Performing greedy analysis")
        models = self.cfg.models
        partnum = len(self.cfg.partitions)

        self.total_scheme_num = submodels.count_greedy_schemes(partnum)
        log.info("This will result in a maximum of %s schemes being created", self.total_scheme_num)

        self.total_subset_num = submodels.count_greedy_subsets(partnum)
        log.info("PartitionFinder will have to analyse a maximum of %d subsets of sites to complete this analysis" %(self.total_subset_num))

        if self.total_subset_num>10000:
            log.warning("%d is a lot of subsets, this might take a long time to analyse", self.total_subset_num)
            log.warning("Perhaps consider using a different search scheme instead (see Manual)")

        #clear any schemes that are currently loaded
        # TODO Not sure we need this...
        self.cfg.schemes.clear_schemes()        
                
        #start with the most partitioned scheme
        start_description = range(len(self.cfg.partitions))
        start_scheme = scheme.create_scheme(self.cfg, 1, start_description)
        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        result = self.analyse_scheme(start_scheme, models)
        
        best_result = result
        best_score  = get_score(result, self.cfg.model_selection)
                         
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
                print lumped_description
                lumped_scheme = scheme.create_scheme(self.cfg, cur_s, lumped_description)
                cur_s += 1
                result = self.analyse_scheme(lumped_scheme, models)
                new_score = get_score(result, self.cfg.model_selection)

                if best_lumping_score==None or new_score < best_lumping_score:
                    best_lumping_score  = new_score
                    best_lumping_result = result
                    best_lumping_scheme = lumped_scheme
                    best_lumping_desc   = lumped_description

            if best_lumping_score < best_score:
                best_scheme = best_lumping_scheme
                best_score  = best_lumping_score
                best_result = best_lumping_result
                start_description = best_lumping_desc               
                if len(set(best_lumping_desc)) == 1: #then it's the scheme with everything equal, so quit
                    break
                step += 1

            else:
                break

        log.info("Greedy algorithm finished after %d steps" % step)
        log.info("Highest scoring scheme is scheme %s, with %s score of %.3f"
                 %(best_result.scheme.name, self.cfg.model_selection, best_score))

        self.best_result = best_result


    def report(self):
        txt = "Best scheme according to Greedy algorithm, analysed with %s"
        best = [(txt % self.cfg.model_selection, self.best_result)]
        self.rpt.write_best_schemes(best)
        self.rpt.write_all_schemes(self.results)

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

