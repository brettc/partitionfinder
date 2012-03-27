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

# Setup the logging. We don't want anything going to the output. Just send it
# out to a testing.log file in the tests folder
import logging
handler = logging.FileHandler(os.path.join(TEST_PATH, 'testing.log'),
        mode='w', encoding=None, delay=False)
formatter = logging.Formatter(
        "%(levelname)-8s | %(asctime)s | %(name)-10s | %(message)s")
handler.setFormatter(formatter)
logging.getLogger("").addHandler(handler)
logging.getLogger("").setLevel(logging.INFO)
del handler
del formatter
