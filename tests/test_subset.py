import pytest
from partfinder.subset import Subset
from partfinder.scheme import Scheme, SchemeError
from partfinder.config import Configuration

# TODO: put all of the subset splitting joining tests in here

def test_identity():
    c = Configuration()
    c.init()

    s1 = Subset(c, set(range(10)))
    s2 = Subset(c, set(range(20)))
    s3 = Subset(c, set(range(10, 20)))
    s4 = Subset(c, set(range(20)))

    # Not just equal BUT THE SAME (see the __new__ member of the class Subset)
    assert s1 is not s2
    assert s1 is not s3
    assert s2 is s4


def test_overlap(caplog):
    c = Configuration()
    c.init()

    s1 = Subset(c, set(range(10)))
    s2 = Subset(c, set(range(10, 20)))
    s3 = Subset(c, set(range(9, 20)))

    # This should be okay...
    Scheme(c, 'a', [s1, s2])

    # This isn't
    with pytest.raises(SchemeError):
        Scheme(c, 'a', [s1, s3])
    assert "contains overlapping" in caplog.text()




