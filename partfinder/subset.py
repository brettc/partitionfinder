import logging
log = logging.getLogger("subset")
import os
import weakref

import alignment

class SubsetError(Exception):
    pass

class Subset(object):
    """A Subset of Partitions
    """

    # TODO: changes this to AllSubsets?
    _cache = weakref.WeakValueDictionary()

    # CLEVER BIT:
    # Return the SAME subset if the partitions are identical
    # This is basically a pythonized factory
    # See here: http://codesnipers.com/?q=python-flyweights
    def __new__(cls, *parts):
        cacheid = frozenset(parts)
        obj = Subset._cache.get(cacheid, None)
        if not obj:
            obj = object.__new__(cls)
            Subset._cache[cacheid] = obj

            # Error checking....
            tempparts = set()
            for p in parts:
                if p.partition_set is None:
                    log.error("You cannot add a Partition to a Subset until "
                              "the Partition belongs to a PartitionSet")
                    raise SubsetError

                if p in tempparts:
                    log.error("%s is duplicated in a Subset", p)
                    raise SubsetError

                tempparts.add(p)

            obj.partitions = cacheid
            
            # Append all of the columns in the partition -- I think these are
            # useful...
            obj.columns = []
            obj.columnset = set()
            for p in parts:
                obj.columns += p.columns
                obj.columnset |= p.columnset
            obj.columns.sort()

            obj.has_analysis = False
            log.debug("Created %s", obj)
        # else:
            # log.debug("Reused %s", obj)

        return obj

    # def __init__(self, *parts):
        # Everything is relegated to above...

    def __str__(self):
        return "Subset(%s)" % ", ".join([str(p) for p in self.partitions])

    @property
    def name(self):
        s = sorted([p.name for p in self.partitions])
        return '-'.join(s)

    def __iter__(self):
        return iter(self.partitions)

if __name__ == '__main__':
    import logging
    import tempfile
    logging.basicConfig(level=logging.DEBUG)
    import config
    from partition import Partition

    tmp = tempfile.mkdtemp()
    config.initialise(tmp, True)

    pa = Partition('a', (1, 10, 3))
    pb = Partition('b', (2, 10, 3))
    pc = Partition('c', (3, 10, 3))

    s1 = Subset(pa, pb)
    s2 = Subset(pa, pb)
    s3 = Subset(pa, pc)
    s4 = Subset(pa, pb)
    print s1 is s2
    print s1
    # s2 = Subset(pa, pb, pc)
    # s3 = Subset(pc)

    # print s1.name
    # print s2.name


