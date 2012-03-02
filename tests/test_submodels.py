from basetest import *
from partfinder.submodels import get_submodels, count_all_schemes 

class TestSubmodels(PartitionFinderTestCase):

    def test_consistency(self):
        known_results = [
            [0, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0],
            [0, 0, 1, 1],
            [0, 0, 1, 2],
            [0, 1, 0, 0],
            [0, 1, 0, 1],
            [0, 1, 0, 2],
            [0, 1, 1, 0],
            [0, 1, 1, 1],
            [0, 1, 1, 2],
            [0, 1, 2, 0],
            [0, 1, 2, 1],
            [0, 1, 2, 2],
            [0, 1, 2, 3],
        ]
        submodels = get_submodels(4)
        self.assertEqual(submodels, known_results)

    def test_scheme_lengths(self):
        self.assertEqual(count_all_schemes(1), 1)
        self.assertEqual(count_all_schemes(5), 52)
        self.assertEqual(count_all_schemes(10), 115975)

if __name__ == '__main__':
    unittest.main()
