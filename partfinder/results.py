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

    def compare(self, cfg):
        pth = self.get_dump_path(cfg) 
        if not os.path.exists(pth):
            log.error("Previous results file not found at '%s'. "
                      "Did you run --dump-results previously?", pth)
            raise ComparisonError

        log.info("Loading old results from '%s'", pth)
        f = open(pth, 'rb')
        old_results = pickle.load(f)

        log.info("Comparing results...")
        # Now do the comparison
        for old, new in zip(old_results, self.scheme_results):
            if not old.compare(new):
                # TODO Make this better
                log.error("Results were not equal")
                raise ComparisonError

        log.info("Results were equal to dumped results")

