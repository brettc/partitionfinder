import unittest
import numpy as np
import random
import scipy
from scipy import stats

def make_mat():
	a = np.random.randint(5, size=(85, 300))
	return a

def read_in_array():
	mat = np.genfromtxt("mat.nex", skip_header=1, delimiter=',')
	return mat

def rescale_mat():
	"""RAxML expects symbols to be constant across columns, so we rescale the matrices to meet this reuirement"""
	a = read_in_array()
	a = np.ma.masked_array(a,np.isnan(a))
#Solution with scipy: Times at 249ms per loop with 100 x 10000 matrix. 2.73s per loop for 1000 x 10000 and runs out of memory at 1000 x 10000.
#	rescaled_matrix = scipy.stats.rankdata(a,'dense').reshape(a.shape)-1
#Solution 2: This runs into a memory error on my computer with a 1000 x 100000 or 1000 x 10000 data set. Times at 410ms per loop with 100 x 10000 matrix. I think this is not our solution.
#	uniq, inv = np.unique(a, return_inverse=True)
#	rescaled_matrix = inv.reshape(a.shape)
#Solution 3: 192ms per loop with a 100 x 10000 matrix, 19.2s per loop with 1000 x 100000. 
	uniq = np.unique(a)
	rescaled_matrix = uniq.searchsorted(a)
	return a, rescaled_matrix

def fix_missing():
	"""Rescaled matrices are don't have missing data symbols, as numpy does not use them, so we need to reinsert those"""
	mat_list = rescale_mat()
	a = mat_list[1]
	a = a.tolist()
#timeit solutions, including all the dependent functions: 1000x1000 size: 365ms per loop. 1000x10000: 3.67s per loop. 1000x100000: 36.8s per loop. Seems good.
	a=[['?' if x is None else x for x in i] for i in a]
	return a
	
def test_type_of():
	"""Test to mkae sure matrix being imported is actually a numpy array, otherwise the next step does not work"""
	a = read_in_array()
	assert isinstance(a, np.ndarray)
	print 'Morphology data passes initial format check.'

def test_shape():
	'''check that the rescaled matrix with appropriate missing data symbols matches the input matrix'''
	a = read_in_array()
	b = fix_missing()
	b = np.array(b)
#	print b, b.shape
#	print a, a.shape
	try:
		a.shape == b.shape
		print "Input and output matrices are the same size. All data present and accounted for, captain!"
	except:
		a.shape != b.shape
		print "Something has gone wrong. This is most likely due to non-numeric characters in your data"

def test_rescaling():
	"""check that the rescaling worked: calling unique on the rescaled matrix should yield a sequential series of numbers beginning from zero"""
	b = rescale_mat()
	c = np.unique(b[1])
	try:
		len(c)-1 == max(c)-min(c)
		min(c) == 0
		print "Input matrix appropriately scaled; continuing to write out matrices."
	except:
		len(c)-1 != max(c)-min(c)
		print "These values are not appropriately scaled, please check that your input matrix is comprised of numbers and not letters."
		min(c) != 0
		print "These values are not appropriately scaled, please check that your input matrix is comprised of numbers and not letters."

def test_nans_replaced():
	a = fix_missing()
	new_list = [x for sublist in a for x in sublist]
	unique_vals = set(new_list)
	try:
		None not in unique_vals
		'?' in unique_vals
		print "NaNs removed, matrix has appropriate missing data symbols."
	except:
		None in unique_vals
		print 'There are NoneTypes in this array, which RAxML cannot handle. Please check that there are no non-numeric characters in your input matrix.'


test_nans_replaced()



