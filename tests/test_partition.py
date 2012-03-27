from basetest import *
from partfinder.partition import Partition, PartitionError
from partfinder.config import Configuration

from nose.tools import raises

def test_construction():
    c = Configuration()
    p1 = Partition(c, 'one', (1, 10))
    p2 = Partition(c, 'two', (11, 20))
    p3 = Partition(c, 'three', (21, 21))

    assert p1.columnset == set(range(10))
    assert c.partitions.columnset == set(range(0, 21))

@raises(PartitionError)
def test_overlap():
    c = Configuration()
    p1 = Partition(c, 'one', (1, 10))
    p2 = Partition(c, 'two', (10, 20))

if __name__ == '__main__':
    import nose
    nose.runmodule()
