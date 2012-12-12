import os
from partfinder import main

HERE = os.path.abspath(os.path.dirname(__file__))


def pytest_generate_tests(metafunc):
    if 'test_folder' in metafunc.fixturenames:
        testdirs = []
        for p in os.listdir(HERE):
            if p.startswith('__'):
                continue
            d = os.path.join(HERE, p)
            if os.path.isdir(d):
                testdirs.append(p)
        metafunc.parametrize("test_folder", testdirs)


def test_quick_asss(test_folder):
    "ntoheu oeu oeu"
    full_path = os.path.join(HERE, test_folder)
    main.call_main("DNA", "%s --force-restart" % full_path)




