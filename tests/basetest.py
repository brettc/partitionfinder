import unittest
import sys
import os

# Magically append the root folder so we can run stuff directly from the tests
# folder if we want
def append_path():
    global _root_path, _test_path 
    _test_path = os.path.dirname(os.path.abspath(__file__))
    _root_path, here = os.path.split(_test_path)
    sys.path.append(_root_path)
append_path()

# Setup the logging. We don't want anything going to the output. Just send it
# out to a testing.log file in the tests folder
import logging
logging.getLogger("").addHandler(
    logging.FileHandler(
        os.path.join(_test_path, 'testing.log'),
        mode='w', encoding=None, delay=False))

# Use a base class for all our tests. We can stick useful stuff in it at the
# class level...
class PartitionFinderTestCase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Set up some useful variables
        cls.test_path = _test_path
        cls.cfg_path = os.path.join(cls.test_path, 'cfg')
        cls.phyml_path = os.path.join(cls.test_path, 'phyml')

