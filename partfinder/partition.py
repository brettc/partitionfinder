import logging
log = logging.getLogger("partition")

class PartitionError(Exception):
    pass

def columnset_to_string(colset):
    s = list(colset)
    s.sort()
    # Add one, cos we converted to zero base...
    return ', '.join([str(x+1) for x in s])

class Partition(object):
    """A set of columns from an alignment"""
    def __init__(self, name, *partlist):
        """A named partition

        """
        self.name = name
        parts = []

        # This will get set later, when they are added to PartitionSet
        self.partition_set = None

        # We now need to convert to column definitions. Note that these are
        # zero based, which is not how they are specified in the config. So we
        # must do some fiddling to make sure they are right. In addition, we
        # use range(...) which excludes the final column, whereas the
        # definitions assume inclusive...
        columns = []
        for p in partlist:

            # Make sure it is sensible
            if len(p) < 2 or len(p) > 3:
                log.error("Partition definition '%s' should contain a list of start, a stop, and an optional step", self.name)
                raise PartitionError
            if len(p) == 2:
                start, stop = p
                step = 1
            else:
                start, stop, step = p
            if start >= stop:
                log.error("Partition '%s' has beginning after end (%s > %s)",
                          name, start, stop)
                raise PartitionError

            # Actually, subtracting 1 deals with both issues...
            columns.extend(range(start-1, stop, step))
            parts.append((start, stop, step))

        self.parts = tuple(parts)

        # Normalise it all
        columns.sort()
        columnset = set(columns)

        # If there was any overlap then these will differ...
        if len(columns) != len(columnset):
            log.error("Partition '%s' has internal overlap", name)
            raise PartitionError

        # Both of these are useful?
        self.columns = columns
        self.columnset = columnset


    def __repr__(self):
        outlist = ", ".join(["%s-%s\\%s" % tuple(p) for p in self.parts])
        return "Partition<%s: %s>" % (self.name, outlist)

    def __str__(self):
        outlist = ", ".join(["%s-%s\\%s" % tuple(p) for p in self.parts])
        return "Partition('%s', %s)" % (self.name, outlist)

class PartitionSet(object):
    """The set of all partitions loaded from a configuration file"""
    def __init__(self, *partitions):
        """A set of Partitions"""
        self.parts = {}

        # All of the columns
        self.columns = []
        self.columnset = set()

        # Freeze it afterwards
        self.finalised = False

        for p in partitions:
            self.add_partition(p)
        self.finalise()

    def __str__(self):
        return "PartitionSet(" + ", ".join([str(p) for p in self.parts.values()]) + ")"

    def add_partition(self, p):
        """Check for overlap (= intersection)"""
        if self.finalised:
            log.error("Cannot add partitions when to set that is finalised")
            raise PartitionError

        if p.name in self.parts:
            log.error( "Attempt to add %s which already exists", p)
            raise PartitionError

        overlap = self.columnset & p.columnset
        if overlap:
            log.error(
                "%s overlaps with previous partitions at columns %s",
                p, columnset_to_string(overlap))
            raise PartitionError

        p.partition_set = self
        self.parts[p.name] = p
        self.columns.extend(p.columns)
        self.columns.sort()
        self.columnset |= p.columnset

        log.debug("Created %s", p)

    def finalise(self):
        """Internal check -- just for gaps now"""
        # It is sorted -- so the last one is the biggest
        colmax = self.columns[-1]
        fullset = set(range(colmax))
        leftout = fullset - self.columnset
        if leftout:
            # This does not raise an error
            log.warn(
                "These columns are missing from the partition definitions: %s", 
                columnset_to_string(leftout))

        self.finalised = True

    # We can treat this like a bit like a dictionary
    def __getitem__(self, k):
        return self.parts[k]

    def __contains__(self, k):
        return k in self.parts

    def names(self):
        return self.parts.keys()


# def test_partition():
    # p = Partition("A", [[1, 10, 3]])
    # assert p.columns == [0, 3, 6, 9]
    # >>> p = Partition("A", [[1, 10, 1],[1, 10, 2]])
    # PartitionError
    # """
if __name__ == '__main__':
    import logging
    logging.basicConfig()
    p1 = Partition('one', (1, 10))
    p2 = Partition('two', (11, 20))
    ps = PartitionSet(p1, p2)
    print ps

    # print ps
