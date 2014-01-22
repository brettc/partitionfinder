import unittest
import numpy as np
import random

def make_mat():
	a = np.random.randint(11, size=(5, 5))
	print a
	return a

def rescale_mat():
	a = make_mat()
	a = a - a.min(axis=0)
	print a
	return a


def test_type_of():
	a = make_mat()
	assert isinstance(a, np.ndarray)
	print 'Morphology data passes initial format check.'

def test_lower_lim():
	a = rescale_mat()
	b = a.min(axis = 0)
	print b
	for x in b:
		assert (x == 0), "These values are not appropriately scaled, please check that all characters are numerical."



rescale_mat()



