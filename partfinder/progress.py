import logging
log = logging.getLogger("progress")
import scheme
import subset

class Progress(object):
    """Provide progress reports"""
    def __init__(self):
        self.subsets_done_set = set()
        self.subsets_done = 0
        self.schemes_done = 0

    def begin(self, scheme_count, subset_count=None, approx=False):
        if subset_count == None:
            # TODO: make this nicer 
            subset_count = len(subset.Subset._cache)

        self.scheme_count = scheme_count
        self.subset_count = subset_count

        if approx:
            log.info("This will result in a maximum of "
                     "%s schemes being created", self.scheme_count)
            log.info("PartitionFinder will have to analyse a maximum of "
                     "%d subsets of sites to complete this analysis",
                     self.subset_count)
        else:
            log.info("This will result in %s schemes being created", scheme_count)
            log.info("PartitionFinder will have to analyse %d subsets "
                    "to complete this analysis", subset_count)

        if subset_count > 10000:
            log.warning("%d is a lot of subsets, this might take a long time to analyse", self.total_subset_num)
            log.warning("Perhaps consider using a different search scheme instead (see Manual)")
            
    def another_subset(self, sub):
        if sub in self.subsets_done_set:
            return

        self.subsets_done_set.add(sub)
        self.subsets_done = len(self.subsets_done_set)

        done = float(self.subsets_done) / float(self.subset_count)
        self.percent_done = done * 100.0

        log.info("Analysing subset %d/%d: %%%.2f done", self.subsets_done,
                 self.subset_count,
                 self.percent_done)

    def another_scheme(self, sch):
        self.schemes_done = self.schemes_done + 1        
        log.info("Analysing scheme %d/%d",
                 self.schemes_done, 
                 self.scheme_count)

    def end(self):
        pass

        
