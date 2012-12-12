import os
from partfinder import main
from zipfile import ZipFile


HERE = os.path.abspath(os.path.dirname(__file__))



DNA_names = ["DNA%d" % n for n in range(1, 9)]
prot_names = ["prot%d" % n for n in range(1, 9)]
rerun_ok_names = ["rerun%02d" % int(n) for n in "1 2 3 4 5 6 7 8 17".split()]
rerun_part_error = ["rerun%02d" % int(n) for n in "9 10 11 12 13 16 18 19 20 21".split()]
rerun_anal_error = ["rerun%02d" % int(n) for n in "14 15".split()]


def pytest_generate_tests(metafunc):
    # This function feeds the output of the above function into the tests below
    if "dna_folder" in metafunc.fixturenames:
        metafunc.parametrize("dna_folder", DNA_names)
    if "prot_folder" in metafunc.fixturenames:
        metafunc.parametrize("prot_folder", prot_names)
    if "rerun_ok_folder" in metafunc.fixturenames:
        metafunc.parametrize("rerun_ok_folder", rerun_ok_names)
    if "rerun_part_folder" in metafunc.fixturenames:
        metafunc.parametrize("rerun_part_folder", rerun_part_error)
    if "rerun_anal_folder" in metafunc.fixturenames:
        metafunc.parametrize("rerun_anal_folder", rerun_anal_error)


def test_dna(dna_folder):
    full_path = os.path.join(HERE, dna_folder)
    main.call_main("DNA", "%s --force-restart" % full_path)


def test_prot(prot_folder):
    full_path = os.path.join(HERE, prot_folder)
    main.call_main("protein", "%s --force-restart" % full_path)


def load_rerun(pth):
    dna3 = ZipFile(os.path.join(HERE, 'DNA3-analysis.zip'))
    dna3.extractall(pth)


def test_rerun_ok(rerun_ok_folder):
    full_path = os.path.join(HERE, rerun_ok_folder)
    load_rerun(full_path)
    main.call_main("DNA", "%s --force-restart" % full_path)



