from partfinder import config
import os
import fnmatch

HERE = os.path.abspath(os.path.dirname(__file__))


def get_file_names():
    names = os.listdir(HERE)
    names = fnmatch.filter(names, "*.cfg")
    names.sort()
    return names


def pytest_generate_tests(metafunc):
    # This function feeds the output of the above function into the tests below
    if "config_file" in metafunc.fixturenames:
        metafunc.parametrize("config_file", get_file_names())


def test_config(config_file):
    pth = os.path.join(HERE, config_file)
    c = config.Configuration()
    c.load(pth)
