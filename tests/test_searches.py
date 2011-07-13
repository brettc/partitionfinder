from basetest import *
import os, shutil

from partfinder import analysis, config

class TestTheLot(PartitionFinderTestCase):

    def load_cfg_and_run(self, name):
        try:
            pth = os.path.join(self.test_path, "test_searches", name)
            cfg = config.Configuration()
            cfg.load_base_path(pth)
            anal = analysis.Analysis(cfg, True)
            anal.do_analysis()
        finally:
            # Always do this
            shutil.rmtree(cfg.output_path)

    def test_search_all(self):
        '''Run a full analysis with search=all'''
        self.load_cfg_and_run("all")

    def test_search_greedy(self):
        '''Run a full analysis with search=greedy'''
        self.load_cfg_and_run("greedy")

    def test_search_user(self):
        '''Run a full analysis with search=user'''
        self.load_cfg_and_run("user")

    def test_search_maclinebreaks(self):
        '''Load and run with mac linebreaks in input'''
        self.load_cfg_and_run("maclinebreaks")

    def test_search_maclinebreaks(self):
        '''Load and run with windows linebreaks in input'''
        self.load_cfg_and_run("maclinebreaks")

if __name__ == '__main__':
    unittest.main()
