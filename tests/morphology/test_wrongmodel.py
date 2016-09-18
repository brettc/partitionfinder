import pytest
import os
from partfinder import main, util



def test_mixed():
	'''This test should fail, as we've given PFinder a non-morph model'''
	HERE = os.path.abspath(os.path.dirname(__file__))
	full_path = os.path.join(HERE, "wrongmodel")
	with pytest.raises(util.PartitionFinderError):
		main.call_main("morphology", '--no-ml-tree --min-subset-size 1 --raxml "%s"' % full_path)

	
		
		