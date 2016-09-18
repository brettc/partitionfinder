import os
from partfinder import main

HERE = os.path.abspath(os.path.dirname(__file__))


def get_test_folders():
    testdirs = []
    # Get everything that is here...
    for p in os.listdir(HERE):
        if p.startswith('__'):
            continue
        d = os.path.join(HERE, p)
        if os.path.isdir(d):
            testdirs.append(p)
    return testdirs


def pytest_generate_tests(metafunc):
    # This function feeds the output of the above function into the tests below
    if 'test_folder' in metafunc.fixturenames:
        metafunc.parametrize("test_folder", get_test_folders())


def test_quick_analysis(test_folder):
    full_path = os.path.join(HERE, test_folder)
    main.call_main("DNA", '--no-ml-tree "%s" --force-restart' % full_path)
