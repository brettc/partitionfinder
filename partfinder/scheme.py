import logging
log = logging.getLogger("scheme")

from subset import Subset

class SchemeError(Exception):
    pass

class Scheme(object):
    def __init__(self, name, *subsets):
        """A set of sets of partitions"""
        # self.schemeset = schemeset
        self.name = name

        for s in subsets:
            self.subsets.add((subset)

        # Now check that the scheme is complete without duplications...
        # Gather all the names ...
        all_names = []
        for subset in self.subsets:
            all_names += list(subset.subset_id)

        # ... check for duplicates
        all_names_set = set()
        duplicates = []
        for nm in all_names:
            if nm in all_names_set:
                duplicates.append(nm)
            else:
                all_names_set.add(nm)

        if duplicates:
            log.error("Scheme '%s' contains duplicate partitions: %s", 
                      name,
                      ', '.join(duplicates)
                     )
            raise SchemeError

        # Check for missing
        all_partitions = set(self.schemeset.partitions.names())
        missing = all_partitions - all_names_set
        if missing:
            log.error("Scheme '%s' is missing partitions: %s", 
                      name,
                      ', '.join(list(missing))
                     )
            raise SchemeError

        # Finally, add it to the schemeset
        log.debug("Creating Scheme '%s'", name)
        schemeset.add_scheme(self)

    def make_subset_id(self, subset_def):
        """Check to make sure the partitions exist"""
        parts = self.schemeset.partitions
        for partname in subset_def:
            # Check it is a valid partition
            if partname not in parts:
                log.error(
                    "Creating scheme '%s': '%s' is not a defined partition",
                    name, partname)
                raise SchemeError
        return frozenset(subset_def)

class SchemeSet(object):
    """All the schemes added, and also a list of all unique subsets"""
    def __init__(self, partitions):
        """A collection of schemes"""
        self.partitions = partitions
        self.schemes = {}
        self.subsets = {}

    def get_subset(self, schemename, partlist):
        """Return an existing subset or make a new one"""

        partset = set()
        # This is laborious, but we can error check
        for p in partlist:
            if p in partset:
                log.error(
                    "Partitions '%s' is duplicated within a subset in Scheme '%s'",
                    p.name, schemename)
                raise PartitionError
            partset.add(p)

        # Freeze it. We can use the set of partitions as a key
        partset = frozenset(partset)
        if partset in self.subsets:
            return self.subsets[partset]

        # Create a new subset
        sub = Subset(partset)
        self.subsets[partset] = sub
        return sub

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
