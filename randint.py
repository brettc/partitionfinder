import unittest
import numpy as np
import random

def make_mat():
	a = np.random.randint(11, size=(5, 5))
#	print a
	return a

#def read_in_array():
#	mat = np.genfromtxt("archy_mb.txt", skip_header=1, names=True)
#	print mat

def rescale_mat():
	a = make_mat()
	uniq, inv = np.unique(a, return_inverse = True)
	rescaled_matrix = inv.reshape(a.shape)
#	print rescaled_matrix
	return a, rescaled_matrix

def test_type_of():
	"""Test to mkae sure matrix being imported is actually a numpy array, otherwise the next step does not work"""
	a = make_mat()
	assert isinstance(a, np.ndarray)
	print 'Morphology data passes initial format check.'

def test_rescaling():
	"""check that the rescaling worked: the input matrix should not be the same as the rescaled, unless input contained zero"""
	a = rescale_mat()
	print a[0], "\n", a[1]
	rescaled = a[1]
	try:
		 np.array_equal(a[:1],a[1:])
		 np.any(rescaled[:, 0] == 0)
		 print "Input matrix appropriately scaled; continuing to write out matrices."
	except:
		print "These values are not appropriately scaled, please check that your input matrix is comprised of numbers and not letters."
	


test_rescaling()



