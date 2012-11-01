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
log = logging.getLogger("progress")


class Progress(object):
    def __init__(self, cfg):
        self.cfg = cfg
        self.cfg.progress = self

    def begin(self, scheme_count, subset_count):
        pass

    def next_scheme(self):
        pass

    def subset_begin(self, sub):
        pass

    def subset_done(self, sub):
        pass

    def end(self):
        pass

class NoProgress(Progress):
    pass

class TextProgress(Progress):

    def begin(self, scheme_count, subset_count):
        self.scheme_count = scheme_count
        self.subset_count = subset_count
        self.schemes_analysed = 0
        self.subsets_analysed = 0

    def next_scheme(self):
        self.schemes_analysed += 1
        log.info("Analysing scheme %d/%d", self.schemes_analysed,
                 self.scheme_count)

    def subset_begin(self, sub):
        log.info("Begin analysing subset %d/%d",
                 self.subsets_analysed, self.subset_count)
        pass

    def subset_done(self, sub):
        self.subsets_analysed += 1
        log.info("Finished analysing subset %d/%d",
                 self.subsets_analysed, self.subset_count)
        # percent_done = float(self.subsets_analysed)*100.0/float(self.total_subset_num)
        # log.info("Analysing subset %d/%d: %.2f%s done" %
                     # (self.subsets_analysed,self.total_subset_num,
                      # percent_done, r"%"))

    def end(self):
        pass
