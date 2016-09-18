import os
import pytest
from partfinder import main, util, analysis, config
from zipfile import ZipFile, ZIP_DEFLATED

HERE = os.path.abspath(os.path.dirname(__file__))

organise_tests = {
    "dna"               : ["DNA%d" % n for n in range(1, 9)],
    "prot"              : ["prot%d" % n for n in range(1, 9)],
    "rerun_success"     : ["rerun%02d" % int(n) for n in "1 2 3 4 5 6 7 8".split()],
    "rerun_pf_error"    : ["rerun%02d" % int(n) for n in "9 10 11 12 13 20 21".split()],
    "rerun_ana_error"   : ["rerun%02d" % int(n) for n in "14 15".split()],
}


def pytest_generate_tests(metafunc):
    # This function feeds the output of the above function into the tests below
    for test_type, folders in organise_tests.items():
        if test_type in metafunc.fixturenames:
            metafunc.parametrize(test_type, folders)


def test_dna(dna):
    full_path = os.path.join(HERE, dna)
    main.call_main("DNA", '--no-ml-tree --compare "%s"' % full_path)


def test_prot(prot):
    full_path = os.path.join(HERE, prot)
    main.call_main("protein", '--no-ml-tree --compare "%s"' % full_path)


def load_rerun(pth):
    dna3 = ZipFile(os.path.join(HERE, 'DNA3-analysis.zip'))
    dna3.extractall(pth)


def test_rerun_success(rerun_success):
    full_path = os.path.join(HERE, rerun_success)
    load_rerun(full_path)
    main.call_main("DNA", '--no-ml-tree "%s"' % full_path)


def test_rerun_pf_error(rerun_pf_error):
    full_path = os.path.join(HERE, rerun_pf_error)
    load_rerun(full_path)
    with pytest.raises(util.PartitionFinderError):
        main.call_main("DNA", '--no-ml-tree "%s"' % full_path)


def test_rerun_analysis_error(rerun_ana_error):
    full_path = os.path.join(HERE, rerun_ana_error)
    load_rerun(full_path)
    with pytest.raises(analysis.AnalysisError):
        main.call_main("DNA", '--no-ml-tree "%s"' % full_path)

