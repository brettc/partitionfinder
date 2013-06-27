from partfinder.partition import Partition
from partfinder.subset import Subset
from partfinder.config import Configuration

# TODO: put all of the subset splitting joining tests in here

def test_identity():
    c = Configuration()

    s1 = Subset(c, set(range(10)))
    s2 = Subset(c, set(range(20)))
    s3 = Subset(c, set(range(10, 20)))
    s4 = Subset(c, set(range(20)))

    # Not just equal BUT THE SAME
    assert s1 is not s2
    assert s1 is not s3
    assert s2 is s4
