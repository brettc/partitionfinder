import unittest
import sys
import os

# Magically append the root folder so we can run stuff directly from the tests
# folder if we want
TEST_PATH = os.path.dirname(os.path.abspath(__file__))
ROOT_PATH, here = os.path.split(TEST_PATH)
sys.path.append(ROOT_PATH)
CFG_PATH = os.path.join(TEST_PATH, 'cfg')
PHYML_PATH = os.path.join(TEST_PATH, 'phyml')
ANALYSIS_PATH = os.path.join(TEST_PATH, 'analysis')
del here

# Setup the logging. We don't want anything going to the output. Just send it
# out to a testing.log file in the tests folder
import logging
logging.getLogger("").addHandler(
    logging.FileHandler(
        os.path.join(TEST_PATH, 'testing.log'),
        mode='w', encoding=None, delay=False))
logging.getLogger("").setLevel(logging.INFO)

# Use a base class for all our tests.
class PartitionFinderTestCase(unittest.TestCase):
    pass

