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
#program, the RAxML program, and the PyParsing library,
#all of which are protected by their own licenses and conditions, using
#PartitionFinder implies that you agree with those licences and conditions as well.

import logging
log = logging.getLogger("results")

import os
import cPickle as pickle

from util import PartitionFinderError


class ComparisonError(PartitionFinderError):
    pass


class AnalysisResults(object):
    """
    This stores the results, keeping only the winning scheme.
    """

    MAX_ERROR = .01

    def __init__(self, model_selection):
        self.model_selection = model_selection
        self.best_score = None
        self.best_result = None
        self.best_scheme = None

    def add_scheme_result(self, sch, result):
        score = result.score
        if self.best_score is None or score < self.best_score:
            self.best_score = score
            self.best_result = result
            self.best_scheme = sch

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

        # TODO: What's going on here?
        if lnl_err > AnalysisResults.MAX_ERROR and\
            aic_err > AnalysisResults.MAX_ERROR and\
            aicc_err > AnalysisResults.MAX_ERROR and\
                bic_err > AnalysisResults.MAX_ERROR:
            return True

        return False

    def compare(self, cfg):
        """We only compare the best result!"""

        pth = self.get_dump_path(cfg)
        if not os.path.exists(pth):
            log.error("Previous results file not found at '%s'. "
                      "Did you run --dump-results previously?", pth)
            raise ComparisonError

        log.info("Loading old results from '%s'", pth)
        f = open(pth, 'rb')
        old = pickle.load(f)
        f.close()

        log.info("Comparing results...")
        # Now do the comparison

        if self.compare_schemes(self.best_scheme, old.best_scheme):
            log.error("Best Scheme Results were more than acceptable value of %s", AnalysisResults.MAX_ERROR)
            log.error("%s to %s", self, old)
            raise ComparisonError
        else:
            log.info(
                "All results were within an acceptable %s of the dumped results",
                AnalysisResults.MAX_ERROR)
