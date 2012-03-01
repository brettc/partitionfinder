from basetest import *
from partfinder.partition import Partition
from partfinder.subset import Subset
from partfinder.config import Configuration

class TestSubset(PartitionFinderTestCase):

    def test_identity(self):
        c = Configuration()

        pa = Partition(c, 'a', (1, 10, 3))
        pb = Partition(c, 'b', (2, 10, 3))
        pc = Partition(c, 'c', (3, 10, 3))

        s1 = Subset(pa, pb)
        s2 = Subset(pa, pb)
        s3 = Subset(pa, pc)
        s4 = Subset(pa, pb)

        assert s1 is s2

if __name__ == '__main__':
    unittest.main()
