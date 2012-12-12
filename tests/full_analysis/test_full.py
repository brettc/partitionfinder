import os
import pytest
from partfinder import main, util, analysis, config
from zipfile import ZipFile

HERE = os.path.abspath(os.path.dirname(__file__))


DNA_names = ["DNA%d" % n for n in range(1, 9)]
prot_names = ["prot%d" % n for n in range(1, 9)]
rerun_success_names = ["rerun%02d" % int(n) for n in "1 2 3 4 5 6 7 8 17".split()]
rerun_config_error = ["rerun%02d" % int(n) for n in "9 10 11 12 13 16".split()]
rerun_analysis_error = ["rerun%02d" % int(n) for n in "14 15".split()]
rerun_pf_error = ["rerun%02d" % int(n) for n in "18 19 20 21".split()]

def pytest_generate_tests(metafunc):
    # This function feeds the output of the above function into the tests below
    if "dna_folder" in metafunc.fixturenames:
        metafunc.parametrize("dna_folder", DNA_names)
    if "prot_folder" in metafunc.fixturenames:
        metafunc.parametrize("prot_folder", prot_names)
    if "rerun_success_folder" in metafunc.fixturenames:
        metafunc.parametrize("rerun_success_folder", rerun_success_names)
    if "rerun_pf_error_folder" in metafunc.fixturenames:
        metafunc.parametrize("rerun_pf_error_folder", rerun_pf_error)
    if "rerun_analysis_error_folder" in metafunc.fixturenames:
        metafunc.parametrize("rerun_analysis_error_folder", rerun_analysis_error)
    if "rerun_config_error_folder" in metafunc.fixturenames:
        metafunc.parametrize("rerun_config_error_folder", rerun_config_error)


def test_dna(dna_folder):
    full_path = os.path.join(HERE, dna_folder)
    main.call_main("DNA", "%s" % full_path)


def test_prot(prot_folder):
    full_path = os.path.join(HERE, prot_folder)
    main.call_main("protein", "%s" % full_path)


def load_rerun(pth):
    dna3 = ZipFile(os.path.join(HERE, 'DNA3-analysis.zip'))
    dna3.extractall(pth)


def test_rerun_success(rerun_success_folder):
    full_path = os.path.join(HERE, rerun_success_folder)
    load_rerun(full_path)
    main.call_main("DNA", "%s" % full_path)


def test_rerun_pf_error(rerun_pf_error_folder):
    full_path = os.path.join(HERE, rerun_pf_error_folder)
    load_rerun(full_path)
    with pytest.raises(util.PartitionFinderError):
        main.call_main("DNA", "%s" % full_path)


def test_rerun_analysis_error(rerun_analysis_error_folder):
    full_path = os.path.join(HERE, rerun_analysis_error_folder)
    load_rerun(full_path)
    with pytest.raises(analysis.AnalysisError):
        main.call_main("DNA", "%s" % full_path)


def test_rerun_config_error(rerun_config_error_folder):
    full_path = os.path.join(HERE, rerun_config_error_folder)
    load_rerun(full_path)
    with pytest.raises(config.ConfigurationError):
        main.call_main("DNA", "%s" % full_path)
