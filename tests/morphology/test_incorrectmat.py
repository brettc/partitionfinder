import pytest
import os
from partfinder import main, util, scheme, analysis_method, config, alignment, results, reporter


def test_incmat():
	'''This test should fail due to incorrect line length'''
	HERE = os.path.abspath(os.path.dirname(__file__))
	full_path = os.path.join(HERE, "incmat")
	with pytest.raises(util.PartitionFinderError):
		main.call_main("morphology", '--no-ml-tree --min-subset-size 1 --raxml "%s"' % full_path)

    
