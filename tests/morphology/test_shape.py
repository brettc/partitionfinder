import sys
import numpy as np
from partfinder import rescale_morphology as rem

def test_shape():
	"""check that the size and shape (number of characters and taxa) of the 
    rescaled matrix with appropriate missing data symbols matches the input 
    matrix. Tests that no data has been lost in the functions rescale_mat and 
    fix_missing."""
    a = rem.make_mat()
	b = rem.fix_missing()
	b = rem.np.array(b)
	assert a.shape == b.shape, """Input and output matrices are the same size. 
                                All data present and accounted for, captain!"""
	if not a.shape == b.shape:
		raise AssertionError, """Data has been lost between input and output 
                               matrix. This is most likely due to non-numeric 
                               characters in your data, please check that no 
                               non-numeric or missing data ('-', '?') 
                               characters exist in input file."""

