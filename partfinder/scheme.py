import logging
log = logging.getLogger("scheme")

from subset import Subset

class SchemeError(Exception):
    pass

class Scheme(object):
    def __init__(self, name, *subsets):
        """A set of subsets of partitions"""
        self.name = name
        self.subsets = set()

        partitions = set()
        duplicates = []
        for s in subsets:
            if s in self.subsets:
                duplicates.append(s)
            else:
                self.subsets.add(s)
                partitions |= s.partitions

        if duplicates:
            log.error("Scheme '%s' contains duplicate partitions: %s", 
                      name, ', '.join(duplicates))
            raise SchemeError

        # Hm. It seems this is the only way to get just one item out of a set
        # as pop would remove one...
        pset = iter(partitions).next().partition_set
        # Set difference
        missing = pset.partitions - partitions
        if missing:
            log.error("Scheme '%s' is missing partitions: %s", 
                      name, ', '.join([str(p) for p in missing]))
            raise SchemeError

        log.debug("Created %s", self)

    def __str__(self):
        ss = ', '.join([str(s) for s in self.subsets])
        return "Scheme(%s, %s)" % (self.name, ss)


class SchemeSet(object):
    """All the schemes added, and also a list of all unique subsets"""
    def __init__(self, partitions):
        """A collection of schemes"""
        self.partitions = partitions
        self.schemes = {}
        self.subsets = {}

    def add_scheme(self, scheme):
        if scheme.name in self.schemes:
            log.error("Cannot add two schemes with same name: '%s'" %
                      scheme.name)
            raise SchemeError
        self.schemes[scheme.name] = scheme

    def analyse(self, config):
        # Analyse all of the subsets...
        for s in self.subsets.values():
            print s.create_alignment(config)

    # Easy iteration
    def __iter__(self):
        return iter(self.schemes.itervalues())

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.DEBUG)
    from partition import Partition, PartitionSet
    from subset import Subset

    pa = Partition('a', (1, 10, 3))
    pb = Partition('b', (2, 10, 3))
    pc = Partition('c', (3, 10, 3))
    ps = PartitionSet(pa, pb, pc)
    s = Scheme('x', Subset(pa), Subset(pc), Subset(pb))
    
