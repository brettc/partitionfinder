from basetest import *
import os, shutil, sys
from zipfile import ZipFile

from nose.plugins.attrib import attr
from nose.tools import raises

from partfinder import config, analysis_method, reporter
from partfinder import analysis, scheme, util

FULL_PATH = os.path.join(TEST_PATH, 'full_analysis')

# We grab the name of the path from the function name. A bit sneaky, but it
# avoids screwing up function name / folder correspondence. See here:
# http://code.activestate.com/recipes/66062-determining-current-function-name/
def path_from_function():
    # get name of the caller function (up 1 level in stack frame)
    funcname = sys._getframe(1).f_code.co_name
    # Remove 'test_'
    dirname = funcname[5:]
    return os.path.join(FULL_PATH, dirname)

def cleanup(cfg):
    shutil.rmtree(cfg.full_output_path)
    # TODO CAUSE PROBLEMS ON WINDOWS -- maybe call logging.shutdown()??
    # os.remove(os.path.join(cfg.base_path, 'log.txt'))

def load_cfg_and_run(pth, compare=True, fails=False):

    try:
        # This is a crappy way to do it, but it works
        if "DNA" in pth or "rerun" in pth:
            dt = "DNA"
        else:
            dt = "protein"
        cfg = config.Configuration(dt)
        cfg.load_base_path(pth)
        method = analysis_method.choose_method(cfg.search)
        reporter.TextReporter(cfg)
        meth = method(cfg, force_restart=False, threads=-1)
        results = meth.analyse()
        if compare:
            results.compare(cfg)

        # Note -- we only get here if everything completes
        # Failed tests leave shit around!
    finally:
        # If it fails, but we expect it to, then clean up
        if fails:
            cleanup(cfg)
        cfg.reset()

    if not fails: # Might already be done
        cleanup(cfg)

def load_rerun(pth, fails=False):
    dna3 = ZipFile(os.path.join(FULL_PATH, 'DNA3-analysis.zip'))
    dna3.extractall(pth)
    load_cfg_and_run(pth, compare=False, fails=fails)

# See ./full_analysis/tests.txt for details
# NOTE We could get all of these automatically, but declaring makes it easier
# to play with individual details
# DNA -----------------------------------

@attr('slow', 'DNA')
def test_DNA1():
    load_cfg_and_run(path_from_function())

@attr('slow', 'DNA')
def test_DNA2():
    load_cfg_and_run(path_from_function())

@attr('slow', 'DNA')
def test_DNA3():
    load_cfg_and_run(path_from_function())

@attr('slow', 'DNA')
def test_DNA4():
    load_cfg_and_run(path_from_function())

@attr('slow', 'DNA')
def test_DNA5():
    load_cfg_and_run(path_from_function())

@attr('slow', 'DNA')
def test_DNA6():
    load_cfg_and_run(path_from_function())

@attr('slow', 'DNA')
def test_DNA7():
    load_cfg_and_run(path_from_function())

@attr('slow', 'DNA')
def test_DNA8():
    load_cfg_and_run(path_from_function())

# Protein -----------------------------------
@attr('slow', 'prot', 'ug')
def test_prot1():
    load_cfg_and_run(path_from_function())

@attr('slow', 'prot')
def test_prot2():
    load_cfg_and_run(path_from_function())

@attr('slow', 'prot')
def test_prot3():
    load_cfg_and_run(path_from_function())

@attr('slow', 'prot')
def test_prot4():
    load_cfg_and_run(path_from_function())

@attr('prot')
def test_prot5():
    load_cfg_and_run(path_from_function())

@attr('slow', 'prot')
def test_prot6():
    load_cfg_and_run(path_from_function())

@attr('slow', 'prot')
def test_prot7():
    load_cfg_and_run(path_from_function())

@attr('slow', 'prot')
def test_prot8():
    load_cfg_and_run(path_from_function())

# Re-running ------------------------------
@attr('rerun')
def test_rerun01():
    load_rerun(path_from_function())

@attr('rerun')
def test_rerun02():
    load_rerun(path_from_function())

@attr('rerun')
def test_rerun03():
    load_rerun(path_from_function())

@attr('rerun')
def test_rerun04():
    load_rerun(path_from_function())

@attr('rerun')
def test_rerun05():
    load_rerun(path_from_function())

@attr('rerun', 'slow')
def test_rerun06():
    load_rerun(path_from_function())

@attr('rerun', 'slow')
def test_rerun07():
    load_rerun(path_from_function())

@attr('rerun')
def test_rerun08():
    load_rerun(path_from_function())

# FAILS
@attr('rerun', 'fails')
@raises(config.ConfigurationError)
def test_rerun09():
    load_rerun(path_from_function(), fails=True)

@attr('rerun', 'fails')
@raises(config.ConfigurationError)
def test_rerun10():
    load_rerun(path_from_function(), fails=True)

@attr('rerun', 'fails')
@raises(config.ConfigurationError)
def test_rerun11():
    load_rerun(path_from_function(), fails=True)

@attr('rerun', 'fails')
@raises(config.ConfigurationError)
def test_rerun12():
    load_rerun(path_from_function(), fails=True)

@attr('rerun', 'fails')
@raises(config.ConfigurationError)
def test_rerun13():
    load_rerun(path_from_function(), fails=True)

@attr('rerun', 'fails')
@raises(analysis.AnalysisError)
def test_rerun14():
    load_rerun(path_from_function(), fails=True)

@attr('rerun', 'fails')
@raises(analysis.AnalysisError)
def test_rerun15():
    load_rerun(path_from_function(), fails=True)

@attr('rerun', 'fails')
@raises(config.ConfigurationError)
def test_rerun16():
    load_rerun(path_from_function(), fails=True)

@attr('rerun', 'slow')
def test_rerun17():
    load_rerun(path_from_function())

@attr('rerun', 'fails')
@raises(util.PartitionFinderError)
def test_rerun18():
    load_rerun(path_from_function(), fails=True)

@attr('rerun', 'fails')
@raises(util.PartitionFinderError)
def test_rerun19():
    load_rerun(path_from_function(), fails=True)

@attr('rerun', 'fails')
@raises(util.PartitionFinderError)
def test_rerun20():
    load_rerun(path_from_function(), fails=True)

@attr('rerun', 'fails')
@raises(util.PartitionFinderError)
def test_rerun21():
    load_rerun(path_from_function(), fails=True)

if __name__ == '__main__':
    nose.runmodule()
