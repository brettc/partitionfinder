"""New tests for PF2"""
import os
import inspect
from partfinder import main
import pytest
from partfinder.util import PartitionFinderError

HERE = os.path.abspath(os.path.dirname(__file__))

# Far too clever function
def path_from_function():
    funname = inspect.stack()[1][3]
    # remove test_
    funname = funname[5:]
    funname = funname.replace('_', '-')
    pth = os.path.join(HERE, funname)
    return pth

# ---------------- SUCCESS ---------------------


# ---------------- ERRORS ---------------------
# Check the exception, and then look in the log output for specific
# details of the exception.

def test_missing_sites_warning(caplog):
    main.call_main("protein", '--no-ml-tree "%s"' % path_from_function())
    assert "These columns are missing from the block definitions" in caplog.text()


def test_overlapping_blocks(caplog):
    with pytest.raises(PartitionFinderError):
        main.call_main("protein", '--no-ml-tree "%s"' % path_from_function())
    assert "following sites overlap" in caplog.text()


