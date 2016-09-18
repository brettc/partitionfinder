import pytest
import shlex
import os
from partfinder import main, util, scheme, analysis_method, config, alignment, results, reporter



def test_multiASC():
	'''This test should pass'''
	HERE = os.path.abspath(os.path.dirname(__file__))
	full_path = os.path.join(HERE, "multiASC")
	main.call_main("morphology", '--no-ml-tree --min-subset-size 1 --raxml "%s"' % full_path)




		