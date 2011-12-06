from basetest import *
import os, shutil

from partfinder import analysis, config

class TestTheLot(PartitionFinderTestCase):

    def load_cfg_and_run(self, name, user_tree=None):
        try:
            pth = os.path.join(ANALYSIS_PATH, name)
            cfg = config.Configuration()
            cfg.load_base_path(pth)
            anal = analysis.Analysis(cfg, True, False)
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

    def test_search_windowslinebreaks(self):
        '''Load and run with windows linebreaks in input'''
        self.load_cfg_and_run("windowslinebreaks")

    def test_search_interleaved(self):
        '''try an interleaved phylip alignment'''
        self.load_cfg_and_run("aln_interleaved")

    def test_search_RY(self):
        '''try an RY-coded phylip alignment'''
        self.load_cfg_and_run("aln_interleaved")


if __name__ == '__main__':
    unittest.main()
