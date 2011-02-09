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
    def __init__(self, name, partlist):
        """A named partition

        The partlist should be a list of column definitions
        e.g. [[1, 100, 3],[100, 155, 2]]
        """
        self.name = name
        self.parts = partlist

        # We now need to convert to column definitions. Note that these are
        # zero based, which is not how they are specified in the config. So we
        # must do some fiddling to make sure they are right. In addition, we
        # use range(...) which excludes the final column, whereas the
        # definitions assume inclusive...
        columns = []
        for p in partlist:
            start, stop, step = p
            # Actually, this is all we need do to deal with both issues...
            start -= 1 
            columns.extend(range(start, stop, step))

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
        outlist = " and ".join(["%s-%s\\%s" % tuple(p) for p in self.parts])
        return "Partition '%s' made of %s" % (self.name, outlist)

class PartitionSet(object):
    """A set of partitions loaded from a configuration file"""
    def __init__(self):
        self.parts = {}

        # All of the columns
        self.columns = []
        self.columnset = set()

    def add_partition(self, p):
        """Check for overlap (= intersection)"""
        if p.name in self.parts:
            log.error(
                "Attempt to add partition '%s' which already exists", 
                p.name)
            raise PartitionError

        overlap = self.columnset & p.columnset
        if overlap:
            log.error(
                "Partition '%s' overlaps with previous partitions at columns %s",
                p.name, columnset_to_string(overlap))
            raise PartitionError

        self.parts[p.name] = p
        self.columns.extend(p.columns)
        self.columns.sort()
        self.columnset |= p.columnset

    def validate(self):
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

    # We can treat this like a bit like a dictionary
    def __getitem__(self, k):
        return self.parts[k]

    def __contains__(self, k):
        return k in self.parts

    def __len__(self):
        return len(self.columnset)

def test_partition():
    p = Partition("A", [[1, 10, 3]])
    assert p.columns == [1, 4, 7]

    # >>> p = Partition("A", [[1, 10, 1],[1, 10, 2]])
    # PartitionError
    # """

