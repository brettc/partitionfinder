import pytest
from partfinder.partition import Partition, PartitionError
from partfinder.config import Configuration


def test_construction():
    c = Configuration()
    p1 = Partition(c, 'one', (1, 10))
    Partition(c, 'two', (11, 20))
    Partition(c, 'three', (21, 21))

    assert p1.columnset == set(range(10))
    assert c.partitions.columnset == set(range(0, 21))


def test_overlap():
    c = Configuration()
    with pytest.raises(PartitionError):
        Partition(c, 'one', (1, 10))
        Partition(c, 'two', (10, 20))
