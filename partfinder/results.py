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
log = logtools.get_logger()

import os
import cPickle as pickle

from util import PartitionFinderError

_check_fields = "lnl aic aicc bic".split()


class ComparisonError(PartitionFinderError):
    pass


class AnalysisResults(object):
    """
    This stores the results, keeping only the winning scheme.
    """

    MAX_ERROR = 10.0

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

    def get_result_fields(self):
        flds = []
        for k in _check_fields:
            flds.append(getattr(self.best_result, k))
        return flds

    def dump(self, cfg):
        pth = self.get_dump_path(cfg)
        log.info("Dumping all results to '%s'", pth)
        f = open(pth, 'wb')
        pickle.dump(self.get_result_fields(), f, -1)

    def compare(self, cfg):
        """We only compare the best result!"""
        pth = self.get_dump_path(cfg)
        if not os.path.exists(pth):
            log.error("Previous results file not found at '%s'. "
                      "Did you run --dump-results previously?", pth)
            raise ComparisonError

        log.info("Loading old results from '%s'", pth)
        f = open(pth, 'rb')
        old_fields = pickle.load(f)
        f.close()

        cur_fields = self.get_result_fields()

        log.info("Comparing results...")
        # Now do the comparison

        errors = 0
        for nm, oldv, curv in zip(_check_fields, old_fields, cur_fields):
            if abs(oldv - curv) > self.MAX_ERROR:
                log.error("Differences were more than acceptable value of %s", AnalysisResults.MAX_ERROR)
                log.error("Old %s value: %s, new %s value %s", nm, oldv, nm, curv)
                errors += 1

        if errors > 0:
            raise ComparisonError
        else:
            log.info(
                "All results were within an acceptable %s of the dumped results",
                AnalysisResults.MAX_ERROR)
