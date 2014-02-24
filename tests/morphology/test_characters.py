import sys
import numpy as np
from partfinder import rescale_morphology as rem

def test_nans_replaced():
	"""Check that all NoneType characters have been replaced with a missing 
    data symbol ('?') that RAxML can actually read. Checks functionality of 
    fix_missing.
"""
	a = rem.fix_missing()
	new_list = [x for sublist in a for x in sublist]
	unique_vals = set(new_list)
	try:
		None not in unique_vals
		'?' in unique_vals
		print "NaNs removed, matrix has appropriate missing data symbols."
	except:
		raise TypeError, """There are NoneType characters in this array, which 
                         RAxML cannot handle. Please check that there are no 
                         non-numeric characters in your input matrix."""
test_nans_replaced()
