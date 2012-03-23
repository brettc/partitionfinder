import sys
import os
import nose

# Magically append the root folder so we can run stuff directly from the tests
# folder if we want
TEST_PATH = os.path.dirname(os.path.abspath(__file__))
ROOT_PATH, here = os.path.split(TEST_PATH)
ROOT_PATH = os.path.relpath(ROOT_PATH)
TEST_PATH = os.path.relpath(TEST_PATH)
sys.path.append(ROOT_PATH)
CFG_PATH = os.path.join(TEST_PATH, 'cfg')
ANALYSIS_PATH = os.path.join(TEST_PATH, 'analysis')
del here

