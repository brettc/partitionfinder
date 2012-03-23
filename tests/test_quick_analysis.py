from basetest import *
import os, shutil

from partfinder import config, analysis_method, reporter

def load_cfg_and_run(name):
    try:
        pth = os.path.join(TEST_PATH, 'quick_analysis', name)
        cfg = config.Configuration()
        cfg.load_base_path(pth)
        method = analysis_method.choose_method(cfg.search)
        rpt = reporter.TextReporter(cfg)
        meth = method(cfg, rpt, True, False)
        results = meth.analyse()
        results.compare(cfg)
    finally:
        # Always do this
        shutil.rmtree(cfg.output_path)

def test_all_analyses():
    analysis_dirs = os.listdir(ANALYSIS_PATH)
    for f in analysis_dirs:
        yield load_cfg_and_run, f

if __name__ == '__main__':
    nose.runmodule()
