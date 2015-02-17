import sys
import numpy as np
import os
from partfinder import rescale_morphology as rem

def test_shape():
    """check that the size and shape (number of characters and taxa) of the 
    rescaled matrix with appropriate missing data symbols matches the input 
    matrix. Tests that no data has been lost in the functions rescale_mat and 
    fix_missing.
"""
    a = rem.read_in_array()
    b = rem.fix_missing()
    b = rem.np.array(b)
    c = a.shape == b.shape
    assert c == True

test_shape()
