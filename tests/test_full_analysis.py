from basetest import *
import os, shutil

from partfinder import config, analysis_method, reporter
from zipfile import ZipFile

from nose.plugins.attrib import attr
from nose.tools import raises

FULL_PATH = os.path.join(TEST_PATH, 'full_analysis')

def load_cfg_and_run(pth, compare=True, restart=True):
    # try:
    # pth = os.path.join(FULL_PATH, name)
    cfg = config.Configuration()
    cfg.load_base_path(pth)
    method = analysis_method.choose_method(cfg.search)
    rpt = reporter.TextReporter(cfg)
    meth = method(cfg, rpt, force_restart=restart, threads=-1)
    results = meth.analyse()
    if compare:
        results.compare(cfg)
    # finally:
        # Always do this
        # shutil.rmtree(cfg.output_path)

# @attr('slow')
# def test_DNA1():
    # load_cfg_and_run('DNA1')

def load_rerun(pth):
    dna3 = ZipFile(os.path.join(FULL_PATH, 'DNA3-analysis.zip'))
    dna3.extractall(pth)
    load_cfg_and_run(pth, compare=False, restart=False)
    # Now run it....

# @attr('rerun')
@raises(config.ConfigurationError)
def load_rerun_configuration_failure(pth):
    load_rerun(pth)

RERUN_SUCCESS = "01 02 03 04 05 06 07 08 17".split()
RERUN_CONFIGURATION_ERROR = "09 10 11 12 13 14 15 16 17 18 19 20 21".split()
# RERUN_FAIL = "14 15".split()

# @attr('slow')
# @raises(config.ConfigurationError)
def test_failed_rerun():
    for t in RERUN_CONFIGURATION_ERROR:
        pth = os.path.join(FULL_PATH, 'rerun%s' % t)
        yield load_failed_rerun, pth

if __name__ == '__main__':
    nose.runmodule()
