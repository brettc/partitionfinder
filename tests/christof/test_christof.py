"""Incorporate some tests from Christoph Mayer (originally in Perl)"""
# TODO This could be automated in the same way "full_analysis'

import os
import inspect
import pytest
from partfinder import main, util, analysis, config, alignment

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

def test_greedy_phyml_dna():
    main.call_main("DNA", '--no-ml-tree "%s"' % path_from_function())

def test_greedy_raxml_dna():
    main.call_main("DNA", '--no-ml-tree "%s" --raxml' % path_from_function())

def test_greedy_phyml_protein():
    main.call_main("protein", '--no-ml-tree "%s"' % path_from_function())

def test_greedy_raxml_protein():
    main.call_main("protein", '--no-ml-tree "%s" --raxml' % path_from_function())

def test_clustering_raxml_dna():
    main.call_main("DNA", '--no-ml-tree "%s" --raxml' % path_from_function())


# ---------------- ERRORS ---------------------
# Check the exception, and then look in the log output for specific
# details of the exception.

def test_alignment_error(caplog):
    with pytest.raises(alignment.AlignmentError):
        main.call_main("protein", '--no-ml-tree "%s"' % path_from_function())
    assert "Site 1000 is specified in [data_blocks], but the alignment only has 949 sites." in caplog.text()

def test_overlap_error(caplog):
    with pytest.raises(util.PartitionFinderError):
        main.call_main("protein", '--no-ml-tree "%s"' % path_from_function())
    assert "sites overlap in your block definitions" in caplog.text()

def test_clustering_phyml_dna(caplog):
    with pytest.raises(util.PartitionFinderError):
        main.call_main("DNA", '--no-ml-tree "%s"' % path_from_function())
    assert "Clustering methods are only available when using raxml" in caplog.text()


