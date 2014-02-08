import sys
import numpy as np
rescale_morphology = sys.path.append('/home/april/project_files/partitions/partitionfinder/partfinder/') 

import rescale_morphology as rem

def test_type_of():
	"""Test to make sure matrix being imported is actually a numpy array, otherwise the next step does not work. Tests the functionality of the conversion between the user's data and a numpy array."""
	a = rem.make_mat()
	try:
		isinstance(a, np.ndarray)
		print 'Morphology data passes initial format check.'
	except	TypeError():
		print "This matrix is not a numpy array."


