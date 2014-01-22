import unittest
import numpy as np
import random

class ZeroException(Exception):
    pass

def make_mat():
	a = np.random.randint(11, size=(5, 5))
	print a
	return a

def read_in_array():
	mat = np.genfromtxt("archy_mb.txt", skip_header=1, names=True)
	print mat

def rescale_mat():
	a = make_mat()
	a = a - a.min(axis=0)
	print a
	return a

def test_type_of():
	"""Test to mkae sure matrix being imported is actually a numpy array, otherwise the next step does not work"""
	a = make_mat()
	assert isinstance(a, np.ndarray)
	print 'Morphology data passes initial format check.'

def test_lower_lim():
	"""check that the rescaling worked: each column must begin from zero"""
	a = rescale_mat()
	b = a.min(axis = 0)
	print b
	for x in b:
		if (x != 0):
			raise ZeroException("These values are not appropriately scaled, please check that your input matrix is comprised of numbers and not letters.")
	print "Input matrix appropriately scaled; continuing to write out matrices."


test_lower_lim()



