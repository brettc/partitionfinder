import logging
log = logging.getLogger("analysis_method")

import os
import scheme
import algorithm
import submodels
from analysis import Analysis, AnalysisError

class UserAnalysis(Analysis):

    def do_analysis(self):
        """Process everything when search=user"""
        scheme_count = len(self.cfg.schemes)
        self.progress.begin(scheme_count)

        if scheme_count > 0:
            for s in self.cfg.schemes:
                 self.analyse_scheme(s, self.cfg.models)
        else:
            log.error("Search set to 'user', but no user schemes detected in .cfg file. Please check.")
            raise AnalysisError

        self.progress.end()

class AllAnalysis(Analysis):

    def do_analysis(self):
        log.info("Performing a Full Analysis")
        partnum = len(self.cfg.partitions)
        scheme_count = submodels.count_all_schemes(partnum)
        subset_count = submodels.count_all_subsets(partnum)
        self.progress.begin(scheme_count, subset_count)

        # Clear any schemes that are currently loaded
        self.cfg.schemes.clear_schemes()

        # Iterate over submodels, which we can turn into schemes afterwards in the loop
        model_iterator = submodels.submodel_iterator([], 1, partnum)

        scheme_name = 1
        for m in model_iterator:
            s = scheme.model_to_scheme(m, scheme_name, self.cfg)
            scheme_name = scheme_name+1
            self.analyse_scheme(s, self.cfg.models)

        self.progress.end()

class GreedyAnalysis(Analysis):

    def do_analysis(self):
        '''A greedy algorithm for heuristic partitioning searches'''
        log.info("Performing greedy analysis")
        partnum = len(self.cfg.partitions)
        scheme_count = submodels.count_greedy_schemes(partnum)
        subset_count = submodels.count_greedy_subsets(partnum)
        self.progress.begin(scheme_count, subset_count)

        model_selection = self.cfg.model_selection

        #clear any schemes that are currently loaded
        # TODO Not sure we need this...
        self.cfg.schemes.clear_schemes()        
                
        #start with the most partitioned scheme
        start_description = range(len(self.cfg.partitions))
        start_scheme = scheme.create_scheme(self.cfg, 1, start_description)
        log.info("Analysing starting scheme (scheme %s)" % start_scheme.name)
        result = self.analyse_scheme(start_scheme, self.cfg.models)
        
        def get_score(my_result):
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

        best_result = result
        best_score  = get_score(result)
                         
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
                lumped_scheme = scheme.create_scheme(self.cfg, cur_s, lumped_description)
                cur_s += 1
                result = self.analyse_scheme(lumped_scheme, self.cfg.models)
                new_score = get_score(result)

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
                 %(best_result.scheme.name, model_selection, best_score))

        self.best_result = best_result

        self.progress.end()

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
    else:
        log.error("Search algorithm '%s' is not yet implemented", search)
        raise AnalysisError
    return method

