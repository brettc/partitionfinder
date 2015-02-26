import os
import pytest
from partfinder import util, config

HERE = os.path.abspath(os.path.dirname(__file__))


def get_file_list(strlist):
    return ["test%02d.cfg" % int(n) for n in strlist.split()]

all_tests = {
    "DNA_success": get_file_list("1 2 3 4 5 6 7 8 9 10 11"),
    "DNA_failure": get_file_list("14"),
    "prot_success": get_file_list("13"),
    "prot_failure": get_file_list("12 15 16"),
}


def pytest_generate_tests(metafunc):
    # This function feeds the output of the above function into the tests below
    for c in all_tests:
        if c in metafunc.fixturenames:
            metafunc.parametrize(c, all_tests[c])


def do_success(ftype, fname):
    pth = os.path.join(HERE, fname)
    c = config.init(datatype=ftype)
    c.load(pth)


def do_failure(ftype, fname):
    pth = os.path.join(HERE, fname)
    c = config.init(datatype=ftype)
    with pytest.raises(util.PartitionFinderError):
        c.load(pth)


def test_DNA_success(DNA_success):
    do_success("DNA", DNA_success)


def test_DNA_failure(DNA_failure):
    do_failure("DNA", DNA_failure)


def test_protein_success(prot_success):
    do_success("protein", prot_success)


def test_protein_failure(prot_failure):
    do_failure("protein", prot_failure)
