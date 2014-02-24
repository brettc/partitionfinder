import sys
import numpy as np
from partfinder import rescale_morphology as rem

def test_rescaling():
	"""check that the rescaling worked: calling unique on the rescaled matrix 
    should yield a sequential series of numbers beginning from zero. 
    Tests function of rescale_mat."""
	b = rem.rescale_mat()
	c = np.unique(b[1])
	try:
		len(c)-1 == max(c)-min(c)
		min(c) == 0
		print "Input matrix appropriately scaled; continuing to write out 
              matrices."
	except:
		raise TypeError, """These values are not appropriately scaled, please 
                         check that your input matrix is comprised of numbers 
                         and not letters."""

