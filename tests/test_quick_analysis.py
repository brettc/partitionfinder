from basetest import *
import os, shutil

from partfinder import config, analysis_method, reporter

QUICK_PATH = os.path.join(TEST_PATH, 'quick_analysis')

def load_cfg_and_run(name):
    try:
        pth = os.path.join(QUICK_PATH, name)
        cfg = config.Configuration()
        cfg.load_base_path(pth)
        method = analysis_method.choose_method(cfg.search)
        rpt = reporter.TextReporter(cfg)
        meth = method(cfg, rpt, True, False)
        results = meth.analyse()
    finally:
        # Always do this
        shutil.rmtree(cfg.full_output_path)
        cfg.reset()

def test_all_analyses():
    analysis_dirs = os.listdir(QUICK_PATH)
    for f in analysis_dirs:
        yield load_cfg_and_run, f

if __name__ == '__main__':
    nose.runmodule()
