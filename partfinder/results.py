#Copyright (C) 2011 Robert Lanfear and Brett Calcott
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
#program and the PyParsing library both of which are protected by their
#own licenses and conditions, using PartitionFinder implies that you
#agree with those licences and conditions as well.

import logging
log = logging.getLogger("results")

import os
import cPickle as pickle

from util import PartitionFinderError
class ComparisonError(PartitionFinderError):
    pass

class AnalysisResults(object):
    """This should hold all the results
    """

    MAX_ERROR = .01

    def __init__(self):
        self.scheme_results = []

    def add_scheme_result(self, result):
        self.scheme_results.append(result)

    def finalise(self):
        self.scheme_results.sort(key=lambda sch: sch.aic)
        self.best_aic  = self.scheme_results[0]
        self.scheme_results.sort(key=lambda sch: sch.aicc)
        self.best_aicc  = self.scheme_results[0]

        # Lets make the default the one sorted by Bic -- so leave it like this
        self.scheme_results.sort(key=lambda sch: sch.bic)
        self.best_bic  = self.scheme_results[0]

    def get_dump_path(self, cfg):
        return os.path.join(cfg.base_path, 'results.bin')

    def dump(self, cfg):
        pth = self.get_dump_path(cfg) 
        log.info("Dumping all results to '%s'", pth)
        f = open(pth, 'wb')
        pickle.dump(self.scheme_results, f, -1)

    def compare_schemes(self, old, new):

        lnl_err = abs(old.lnl - new.lnl)
        aic_err = abs(old.aic - new.aic)
        aicc_err = abs(old.aicc - new.aicc)
        bic_err = abs(old.bic - new.bic)

        # Keep a list of errors
        self.lnl_errs.append(lnl_err)
        self.aic_errs.append(aic_err)
        self.aicc_errs.append(aicc_err)
        self.bic_errs.append(bic_err)

        if lnl_err > AnalysisResults.MAX_ERROR and\
           aic_err > AnalysisResults.MAX_ERROR and\
           aicc_err > AnalysisResults.MAX_ERROR and\
           bic_err > AnalysisResults.MAX_ERROR:
            return True

        return False 

    def output_error_stats(self):
        def mean_and_var(lst):
            m = sum(lst) / float(len(lst))
            vrs = [(x-m)*(x-m) for x in lst]
            v = sum(vrs)
            return m, v

        log.info("LNL Error, Mean %8.8f, Variance %8.8f",
                    *mean_and_var(self.lnl_errs))
        log.info("AIC Error, Mean %8.8f, Variance %8.8f",
                    *mean_and_var(self.aic_errs))
        log.info("AICc Error, Mean %8.8f, Variance %8.8f",
                    *mean_and_var(self.aicc_errs))
        log.info("BIC Error, Mean %8.8f, Variance %8.8f",
                    *mean_and_var(self.bic_errs))

    def compare(self, cfg):
        pth = self.get_dump_path(cfg) 
        if not os.path.exists(pth):
            log.error("Previous results file not found at '%s'. "
                      "Did you run --dump-results previously?", pth)
            raise ComparisonError

        log.info("Loading old results from '%s'", pth)
        f = open(pth, 'rb')
        old_results = pickle.load(f)
        f.close()

        log.info("Comparing results...")
        # Now do the comparison
        err = False
        old_results = dict([(res.scheme.name, res) for res in old_results])
        new_results = dict([(res.scheme.name, res) for res in self.scheme_results])

        self.lnl_errs = []
        self.aic_errs = []
        self.aicc_errs = []
        self.bic_errs = []

        for sch_name in old_results:
            old = old_results[sch_name]
            try:
                new = new_results[sch_name]
            except KeyError:
                log.error("Results do not have the same schemes")
                raise ComparisonError

            if self.compare_schemes(old, new):
                err = True
                log.error("Scheme Results were more than acceptable value of %s", AnalysisResults.MAX_ERROR)
                log.error("%s to %s", old, new)

        self.output_error_stats()
        if err:
            raise ComparisonError
        else:
            log.info("All results were within an acceptable %s of the dumped results",
                     AnalysisResults.MAX_ERROR)

