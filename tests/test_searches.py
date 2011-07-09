from basetest import *
import os, shutil

class TestTheLot(PartitionFinderTestCase):

    def test_search_all(self):
        '''Run a full analysis with search=all'''
        os.system("python PartitionFinder.py tests/test_searches/all")
        shutil.rmtree("tests/test_searches/all/analysis")

    def test_search_greedy(self):
        '''Run a full analysis with search=greedy'''
        os.system("python PartitionFinder.py tests/test_searches/greedy")
        shutil.rmtree("tests/test_searches/greedy/analysis")

    def test_search_user(self):
        '''Run a full analysis with search=user'''
        os.system("python PartitionFinder.py tests/test_searches/user")
        shutil.rmtree("tests/test_searches/user/analysis")

if __name__ == '__main__':
    unittest.main()
