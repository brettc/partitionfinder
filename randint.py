import unittest
import numpy as np
import random
import scipy
from scipy import stats


def make_mat():
	a = np.random.randint(5, size=(1000, 100000))
	return a

def read_in_array():
	mat = np.genfromtxt("mat.nex", skip_header=1, delimiter=',')
	return mat

def rescale_mat():
	a = read_in_array()
#	a = np.array([[1, 2, 10],[1, 2, 99],[1,2,30]])
	print a
#Solution with scipy: Times at 249ms per loop with 100 x 10000 matrix. 2.73s per loop for 1000 x 10000 and runs out of memory at 1000 x 10000.
#	rescaled_matrix = scipy.stats.rankdata(a,'dense').reshape(a.shape)-1
#Solution 2: This runs into a memory error on my computer with a 1000 x 100000 or 1000 x 10000 data set. Times at 410ms per loop with 100 x 10000 matrix. I think this is not our solution.
#	uniq, inv = np.unique(a, return_inverse=True)
#	rescaled_matrix = inv.reshape(a.shape)
#Solution 3: 192ms per loop with a 100 x 10000 matrix, 19.2s per loop with 1000 x 100000. 
	uniq = np.unique(a)
	rescaled_matrix = uniq.searchsorted(a)
	print rescaled_matrix
	return a, rescaled_matrix

def test_type_of():
	"""Test to mkae sure matrix being imported is actually a numpy array, otherwise the next step does not work"""
	a = read_in_array()
	assert isinstance(a, np.ndarray)
	print 'Morphology data passes initial format check.'

def test_rescaling():
	"""check that the rescaling worked: the input matrix should not be the same as the rescaled, unless input contained zero"""
	a = rescale_mat()
#	print a[0], "\n", a[1]
	rescaled = a[1]
	try:
		 np.array_equal(a[:1],a[1:])
		 np.any(rescaled[:, 0] == 0)
		 print "Input matrix appropriately scaled; continuing to write out matrices."
	except:
		print "These values are not appropriately scaled, please check that your input matrix is comprised of numbers and not letters."
	

test_rescaling()



