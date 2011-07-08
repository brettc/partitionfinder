from basetest import *
from partfinder.partition import Partition, PartitionError
from partfinder.config import Configuration

class TestPartition(PartitionFinderTestCase):

    def test_construction(self):
        c = Configuration()
        p1 = Partition(c, 'one', (1, 10))
        p2 = Partition(c, 'two', (11, 20))
        p3 = Partition(c, 'three', (21, 21))

        self.assertEqual(p1.columnset, set(range(10)))
        self.assertEqual(c.partitions.columnset, set(range(0, 21)))

    def test_overlap(self):
        c = Configuration()
        with self.assertRaises(PartitionError):
            p1 = Partition(c, 'one', (1, 10))
            p2 = Partition(c, 'two', (10, 20))


if __name__ == '__main__':
    unittest.main()
