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
        self.subsets_analysed = set()

        if "kmeans" not in self.cfg.search:
            if subset_count > 10000:
                log.warning("""%d is a lot of subsets, this might take a
                long time to analyse. Perhaps consider using a different
                search scheme instead (see the Manual)"""
                % subset_count)

    def next_scheme(self):
        self.schemes_analysed += 1
        #log.info("Analysing scheme %d/%d", self.schemes_analysed,self.scheme_count)

    def subset_begin(self, sub):
        #log.info("Begin analysing subset %s", sub)
        pass

    def subset_done(self, sub):
        old_num_done = len(self.subsets_analysed)
        self.subsets_analysed.add(sub.subset_id)
        num_subs_done = len(self.subsets_analysed)
        if old_num_done != num_subs_done:

            if self.cfg.search in ["kmeans", "krmeans", "user"]:
                # we don't know the total number of possible subsets
                # log.info("Finished subset %d" %(num_subs_done))
                pass
            else:    
                percent_done = (
                    float(num_subs_done) * 100.0) / float(self.subset_count)
                log.info("Finished subset %d/%d, %.2f percent done" %
                         (num_subs_done, self.subset_count, percent_done))

    def end(self):
        pass
