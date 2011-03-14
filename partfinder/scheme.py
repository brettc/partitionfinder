import logging
log = logging.getLogger("scheme")

from partition import all_partitions

class SchemeError(Exception):
    pass

class Scheme(object):
    def __init__(self, name, *subsets):
        """A set of subsets of partitions"""
        global all_schemes
        self.name = name
        self.subsets = set()

        # This one is a set of frozensets of partitions...
        part_subsets = set()

        partitions = set()
        duplicates = []
        for s in subsets:
            for p in s:
                if p in partitions:
                    # This is an error -- we'll collect them up
                    duplicates.append(str(p))
                else:
                    partitions.add(p)
            self.subsets.add(s)
            part_subsets.add(s.partitions)

        self.part_subsets = frozenset(part_subsets)

        # Report the errors 
        if duplicates:
            log.error("Scheme '%s' contains duplicate partitions: %s", 
                      name, ', '.join(duplicates))
            raise SchemeError

        # Hm. It seems this is the only way to get just one item out of a set
        # as pop would remove one...
        pset = iter(partitions).next().partition_set

        # Do a set-difference to see what is missing...
        missing = pset.partitions - partitions
        if missing:
            log.error("Scheme '%s' is missing partitions: %s", 
                      name, ', '.join([str(p) for p in missing]))
            raise SchemeError

        # This locks down whether new partitions can be created.
        if not all_partitions.finalised:
            all_partitions.finalise()
        
        # Now add to all_schemes -- more possibility of errors, see below
        all_schemes.add_scheme(self)
        log.debug("Created %s", self)

    def __iter__(self):
        return iter(self.subsets)

    def __str__(self):
        ss = ', '.join([str(s) for s in self.subsets])
        return "Scheme(%s, %s)" % (self.name, ss)

class AllSchemes(object):
    """All the schemes added, and also a list of all unique subsets"""
    def __init__(self):
        """A collection of schemes"""
        self.schemes_by_name = {}
        self.schemes_by_subsets = {}

    def add_scheme(self, scheme):
        if scheme.name in self.schemes_by_name:
            log.error("Cannot add two schemes with same name: '%s'" %
                      scheme.name)
            raise SchemeError

        if scheme.part_subsets in self.schemes_by_subsets:
            existing_scheme = \
                    self.schemes_by_subsets[scheme.part_subsets]
            log.error(
                "Scheme named %s being added is identical to existing %s",
                scheme.name, existing_scheme)
            raise SchemeError

        self.schemes_by_name[scheme.name] = scheme
        self.schemes_by_subsets[scheme.part_subsets] = scheme

    def __len__(self):
        return len(self.schemes_by_name)
    # Easy iteration
    def __iter__(self):
        return iter(self.schemes_by_name.itervalues())

# Container for all schemes that are created
all_schemes = AllSchemes()

def generate_all_schemes():
    """Convert the abstract schema given by the algorithm into subsets"""
    import subset
    import submodels

    # Make sure that no schemes have been defined!
    if len(all_schemes) > 0:
        log.error("Cannot generate schemes if some already exist!")
        raise SchemeError
    
    partnum = len(all_partitions)
    # Now generate the pattern for this many partitions
    mods = submodels.get_submodels(partnum)
    scheme_name = 1
    for m in mods:
        subs = {}
        # We use the numbers returned to group the different subsets
        for sub_index, grouping in enumerate(m):
            insub = subs.setdefault(grouping, [])
            insub.append(sub_index)
        # We now have what we need to create a subset. Each entry will have a
        # set of values which are the index for the partition
        created_subsets = []
        for sub_indexes in subs.values():
            sub = subset.Subset(*tuple([all_partitions[i] for i in sub_indexes]))
            created_subsets.append(sub)

        Scheme(str(scheme_name), *tuple(created_subsets))
        scheme_name += 1

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.DEBUG)
    from partition import Partition
    from subset import Subset

    pa = Partition('a', (1, 10, 3))
    pb = Partition('b', (2, 10, 3))
    pc = Partition('c', (3, 10, 3))
    pd = Partition('d', (11, 20))
    # s = Scheme('x', Subset(pa, pc), Subset(pb))

    generate_all_schemes()
    # This should give us an error!
    # s = Scheme('y', Subset(pa, pc), Subset(pb))
    
