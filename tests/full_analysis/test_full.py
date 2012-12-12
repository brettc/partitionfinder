import os
from partfinder import main


HERE = os.path.abspath(os.path.dirname(__file__))


# def load_rerun(pth, fails=False):
    # dna3 = ZipFile(os.path.join(FULL_PATH, 'DNA3-analysis.zip'))
    # dna3.extractall(pth)
    # load_cfg_and_run(pth, compare=False, fails=fails)


def get_dna_names():
    return ["DNA%d" % n for n in range(1, 9)]


def get_prot_names():
    return ["prot%d" % n for n in range(1, 9)]


def pytest_generate_tests(metafunc):
    # This function feeds the output of the above function into the tests below
    if "dna_folder" in metafunc.fixturenames:
        metafunc.parametrize("dna_folder", get_dna_names())
    if "prot_folder" in metafunc.fixturenames:
        metafunc.parametrize("prot_folder", get_prot_names())


def test_dna(dna_folder):
    full_path = os.path.join(HERE, dna_folder)
    main.call_main("DNA", "%s --force-restart" % full_path)


def test_prot(prot_folder):
    full_path = os.path.join(HERE, prot_folder)
    main.call_main("protein", "%s --force-restart" % full_path)
