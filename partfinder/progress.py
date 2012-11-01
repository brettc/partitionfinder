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

class TextProgress(Progress):

    def update_subsets(self, sub):
        # Keep people informed about what's going on
        log.info("Analysing subset %s", sub)

    def blarg(self):
        if self.total_subset_num is None:
            self.obj = len(sub._cache)
        old_num_analysed = self.subsets_analysed
        self.subsets_analysed_set.add(sub.name)
        self.subsets_analysed = len(self.subsets_analysed_set)
        if self.subsets_analysed > old_num_analysed: # We've just analysed a subset we haven't seen yet
            percent_done = float(self.subsets_analysed)*100.0/float(self.total_subset_num)
            log.info("Analysing subset %d/%d: %.2f%s done" %
                     (self.subsets_analysed,self.total_subset_num,
                      percent_done, r"%"))
