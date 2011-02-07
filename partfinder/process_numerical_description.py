from os import path
from MinimalSubsets import process_partition_file


def process_numerical_description(description, all_minimal_subsets, debug = 0):
	"""turn a numerical description like this: 012012012 into a nice string like this (part1, part4, part7) (part2, part5, part8) (part3, part6, part9)"""

	if debug>0:
		print description

	#make a dict of minimal subsets keyed by number
	s_dict = {}
	for key in all_minimal_subsets:
		s = all_minimal_subsets[key]
		s_dict[s.line_number] = s

	if debug>0:
		print "s_dict keys", s_dict.keys()
		
	d_list = list(description)
	d_set = list(set(d_list))
	d_set.sort()

	if debug>0:
		print "d_set", d_set
		print "d_list", d_list
	
	subsets = []

	for index in d_set:
		subset = ['(']
		for j in range(len(d_list)): #run through the description e.g. 0 1 2 0 1 2 0 1 2 
			k = d_list[j]
			if debug>0: print index, j, k
			if index==k:
				subset.append(s_dict[int(j)].name)
				subset.append(',')
		subset.append(')')
		subset = ''.join(subset)
		subset = subset.replace(",)", ")")
		subsets.append(subset)
		
	subsets = ','.join(subsets)
		
	if debug>0:
		print subsets
	return subsets

if __name__ == '__main__':
	curdir = path.dirname(path.abspath(__file__))
	partitions_filename = "%s/testfiles/test_partitions" %(curdir)

	all_minimal_subsets = process_partition_file(partitions_filename)

	description = "000111222"
	
	subsets = process_numerical_description(description, all_minimal_subsets)
	
	print subsets
