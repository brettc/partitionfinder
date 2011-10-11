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
log = logging.getLogger("partition")

from util import PartitionFinderError
class PartitionError(PartitionFinderError):
    pass

def columnset_to_string(colset):
    s = list(colset)
    s.sort()
    # Add one, cos we converted to zero base...
    return ', '.join([str(x+1) for x in s])

class PartitionSet(object):
    """The set of all partitions loaded from a configuration file"""
    def __init__(self):
        """A set of Partitions"""
        self.sequence = 0
        self.parts_by_name = {}
        self.parts_by_number = {}
        self.partitions = set()

        # All of the columns
        self.columns = []
        self.columnset = set()

        self.finalised = False

    def __str__(self):
        return "PartitionSet(%s)" % ", ".join([str(p) for p in self.partitions])

    def add_partition(self, p):
        """Check for overlap (= intersection)"""
        if self.finalised:
            log.error("Cannot add partitions after a Scheme has been created")
            raise PartitionError

        if p.name in self.parts_by_name:
            log.error("Attempt to add %s when that name already exists", p)
            raise PartitionError

        overlap = []
        for otherp in self.partitions:
            if p.columnset & otherp.columnset:
                overlap.append(str(otherp))
        if overlap:
            log.error("%s overlaps with previously defined "
                      "partitions: %s",
                      p, ", ".join(overlap))
            raise PartitionError

        # Assign the partition to this set
        p.partition_set = self

        # Make sure we can look up by name
        self.parts_by_name[p.name] = p
        self.parts_by_number[self.sequence] = p
        p.sequence = self.sequence
        self.sequence += 1
        self.partitions.add(p)

        # Merge all the columns
        self.columns.extend(p.columns)
        self.columns.sort()
        self.columnset |= p.columnset

    def finalise(self):
        """Ensure that no more partitions can be added"""
        self.finalised = True

    def check_against_alignment(self, alignment):
        """Check the partition definitions against the alignment"""

        # TODO: pbly should check the converse too -- stuff defined that is
        # missing??
        self.fullset = set(range(0, alignment.sequence_len))
        leftout = self.fullset - self.columnset
        if leftout:
            # This does not raise an error, just a warning
            log.warn(
                "Columns defined in partitions range from %s to %s, "
                "but these columns in the alignment are missing: %s", 
                self.columns[0]+1, self.columns[-1]+1,
                columnset_to_string(leftout))
        
    # We can treat this like a bit like a dictionary
    def __iter__(self):
        return iter(self.partitions)

    def __len__(self):
        return len(self.partitions)

    def __getitem__(self, k):
        if type(k) is int:
            return self.parts_by_number[k]
        return self.parts_by_name[k]

    def __contains__(self, k):
        return k in self.parts_by_name

    def names(self):
        return self.parts_by_name.keys()

class Partition(object):
    """A set of columns from an alignment"""
    def __init__(self, cfg, name=None, *partlist):
        """A named partition

        """
        self.name = name
        description = []

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
                log.error("The Partition '%s' should contain\
                          a list of start, a stop, and an optional step",
                          self.name)
                raise PartitionError
            if len(p) == 2:
                start, stop = p
                step = 1
            else:
                start, stop, step = p
            if start > stop:
                log.error("Partition '%s' has beginning after end (%s > %s)",
                          name, start, stop)
                raise PartitionError

            # Actually, subtracting 1 deals with both issues...
            columns.extend(range(start-1, stop, step))
            description.append((start, stop, step))

        self.description = tuple(description)

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

        cfg.partitions.add_partition(self)
        log.debug("Created %s", self)

    def __repr__(self):
        outlist = ", ".join(["%s-%s\\%s" % tuple(p) for p in self.description])
        return "Partition<%s: %s>" % (self.name, outlist)

    def __str__(self):
        outlist = ", ".join(["%s-%s\\%s" % tuple(p) for p in self.description])
        return "Partition(%s, %s)" % (self.name, outlist)


if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.DEBUG)
    p1 = Partition('one', (1, 10))
    p2 = Partition('two', (11, 20))
    p3 = Partition('three', (21, 21))
    print p3.columnset

    # ps = PartitionSet(p1, p2)

    # print ps
