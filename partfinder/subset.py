import logging
log = logging.getLogger("subset")
import os

import alignment

class SubsetError(Exception):
    pass

results_cache = {}

class Subset(object):
    """A Subset of Partitions
    """

    def __init__(self, *parts):

        partset = set()
        for p in parts:
            if p.partition_set is None:
                log.error("You cannot add a Partition to a Subset until \
                          the Partition belongs to a PartitionSet")
                raise SubsetError

            if p in partset:
                log.error("%s is duplicated in a Subset", p)
                raise SubsetError

            partset.add(p)

        self.partitions = frozenset(partset)
        
        # Append all of the columns in the partition
        self.columns = []
        self.columnset = set()
        for p in parts:
            self.columns += p.columns
            self.columnset |= p.columnset
        self.columns.sort()

        log.debug("Created %s", self)

    def __str__(self):
        return "Subset(%s)" % ", ".join([str(p) for p in self.partitions])

    def make_filename(self):
        s = sorted([p.name for p in self.partitions])
        return '-'.join(s)

    def analyse(self):
        # Check first to see if we've got the results, otherwise calculate and
        # cache them.
        if self.partitions in results_cache:
            log.debug("Returning cached result for %s", self)
            return results_cache[self.partitions]

        log.debug("Calculating result for %s", self)
        result = self._really_analyse()
        results_cache[self.partitions] = result
        return result

    def _really_analyse(self):
        fname = self.make_filename()
        sa = alignment.SubsetAlignment(
            fname, config.settings.source_alignment, self)
        return sa.analyse()

    def __iter__(self):
        return iter(self.partitions)

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.DEBUG)
    import config
    from partition import Partition
    config.initialise("~/tmp", True)
    a = alignment.TestAlignment('test', r"""
>spp1
CTTGAGGTTCAGAATGGTAATGAA------GTGCTGG
>spp2
CTTGAGGTACAAAATGGTAATGAG------AGCCTGG
>spp3
CTTGAGGTACAGAATAACAGCGAG------AAGCTGG
>spp4
CTCGAGGTGAAAAATGGTGATGCT------CGTCTGG
    """)
    config.settings.source_alignment = a

    pa = Partition('a', (1, 10, 3))
    pb = Partition('b', (2, 10, 3))
    pc = Partition('c', (3, 10, 3))

    s1 = Subset(pa, pb)
    s2 = Subset(pa, pb)
    s3 = Subset(pc)

    s1.analyse()
    s2.analyse()
    s3.analyse()


