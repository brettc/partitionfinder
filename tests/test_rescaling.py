import sys
import numpy as np
rescale_morphology = sys.path.append('/home/april/project_files/partitions/partitionfinder/partfinder/') 

import rescale_morphology as rem

def test_shape():
	"""check that the size and shape (number of characters and taxa) of the rescaled matrix with appropriate missing data symbols matches the input matrix. Tests that no data has been lost in the functions rescale_mat and fix_missing."""
	a = rem.make_mat()
	b = rem.fix_missing()
	b = rem.np.array(b)
	assert a.shape == b.shape, "Input and output matrices are the same size. All data present and accounted for, captain!"
	if not a.shape == b.shape:
		raise AssertionError, "Data has been lost between input and output matrix. This is most likely due to non-numeric characters in your data, please check that no non-numeric or missing data ('-', '?') characters exist in input file."

def test_rescaling():
	"""check that the rescaling worked: calling unique on the rescaled matrix should yield a sequential series of numbers beginning from zero. Tests function of rescale_mat."""
	b = rem.rescale_mat()
	c = np.unique(b[1])
	try:
		len(c)-1 == max(c)-min(c)
		min(c) == 0
		print "Input matrix appropriately scaled; continuing to write out matrices."
	except:
		raise TypeError, "These values are not appropriately scaled, please check that your input matrix is comprised of numbers and not letters."

def test_nans_replaced():
	"""Check that all NoneType characters have been replaced with a missing data symbol ('?') that RAxML can actually read. Checks functionality of fix_missing."""
	a = rem.fix_missing()
	new_list = [x for sublist in a for x in sublist]
	unique_vals = set(new_list)
	try:
		None not in unique_vals
		'?' in unique_vals
		print "NaNs removed, matrix has appropriate missing data symbols."
	except:
		raise TypeError, "There are NoneType characters in this array, which RAxML cannot handle. Please check that there are no non-numeric characters in your input matrix."

