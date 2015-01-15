import numpy as np
import os
import random
import glob

HERE = os.path.abspath(os.path.join(os.path.dirname(__file__),"../examples/morphology"))

def rescale_mat(alignment):
    """RAxML expects symbols to be constant across columns, so we rescale
    the matrices to meet this requirement

"""
    align_array = np.ma.masked_array(alignment,np.isnan(alignment))
    uniq = np.unique(align_array)
    rescaled_matrix = uniq.searchsorted(aalign_array)
    return rescaled_matrix]

def fix_missing():
    """Rescaled matrices don't have missing data symbols, as numpy does 
    not use them, so we need to reinsert those
	"""
    mat_list = rescale_mat()
    a = mat_list[1]
    a = a.tolist()
#   timeit solutions, including all the dependent functions: 1000x1000 size: 
#   365ms per loop. 1000x10000: 3.67s per loop. 1000x100000: 36.8s per loop. 
#   Seems good.
    a=[['?' if x is None else x for x in i] for i in a]
    print a
    return(a)

fix_missing()

if __name__ == '__main__':
    rescale_mat()
    fix_missing()
