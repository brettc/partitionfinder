from os import path
from string import strip

class MinimalSubset(object):
	'''an object for a given minimal subset of sites'''
	def __init__(self, from_file, line_number, line):
		self.from_file  	= from_file		#the file that this partition was loaded from
		self.line_number	= line_number	#the order that this partition was in the input file (i.e. 0 would be first in the file)
		self.line			= line			#a string of the line of the partition file e.g. "pos1 = 1-500\3;"
		self.name 			= None			#the name of the subset as listed in the partitions file supplied by the user
		self.description	= None			#the text description of this partition, as supplied by the user in the file
		self.columns 		= None			#a list of alignment columns for this partition NB: these are zero-indexed
		self.errors			= []			#any errrors loading?
		self.process_subset_line()
		self.process_description()

	def __repr__(self): #what gets done when you call print on this class
		return "Minimal_Subset<%s %s %s>" %(self.name, self.description, self.columns)

	def process_subset_line(self):
		"""process a line of a subsets file, to figure out which columns of the alignment are what"""		
		words = self.line.split("=")
		self.name = strip(words[0]) #that bit's easy
		description = strip(words[1]) #the description might have whitespace on the ends
		description = description.strip(";") #some people put semicolons too
		description = description.strip() # and there may be whitespace between the semicolon and the rest of it  
		self.description = description

	def process_description(self):
		'''descriptions can look like this: 1-200\3   2-200\3   210-500 
		   this function turns that into a list of columns in an alignment'''
		description_units = self.description.split()
		all_columns = []
		for unit in description_units:	
			split_unit = unit.split("\\") #split by the \ character, if it's there
			start_stop = split_unit[0].split("-") #extract the start and stop values
			start = int(start_stop[0])-1 # remember we're zero indexed here...
			stop = int(start_stop[1])
			if start>stop:
				self.errors.append("Partition %s seems to have a problem - at least one of the start positions is less than the stop position, please check.\nThe partition description looks like this: %s" %(self.name, self.line))
			if len(split_unit)==2:
				step = int(split_unit[1])
			else:
				step = 1
			unit_columns = range(start, stop, step)
			all_columns = all_columns + unit_columns
		all_columns.sort()
		if len(all_columns)>len(set(all_columns)):
			self.errors.append("Partition %s seems to have duplicate alignment columns specified in it, these will be ignored in the analysis, please check\nThe partition description looks like this: %s" %(self.name, self.line) )
			all_columns = list(set(all_columns))
			all_columns.sort()
		self.columns = all_columns

def process_partition_file(filepath):
	"""input a filepath for a partitions file, and output a dict of MinimumSubset objects"""
	minimum_subsets = {}
	partitions_file = open(filepath, 'r')
	lines = partitions_file.readlines()
	order = 0
	for line in lines:
		if line.count("="): #a simple check to make sure that empty lines aren't a problem
			subset = MinimalSubset(filepath, order, line)
			order = order+1
			minimum_subsets[subset.name] = subset
	return minimum_subsets


if __name__ == "__main__":
	
	curdir = path.dirname(path.abspath(__file__))
	partitions_filename = "%s/testfiles/test_difficult_partition_file" %(curdir)

	subsets = process_partition_file(partitions_filename)
	
	for key in subsets:
		subset = subsets[key]
		print "\n****** subset number %d ******" % subset.line_number
		print "subset.from_file:\t\t%s" % subset.from_file
		print "subset.line_number:\t\t%d" % subset.line_number
		print "subset.line:\t\t%s" % subset.line
		print "subset.name:\t\t'%s'" % subset.name
		print "subset.description:\t\t'%s'" % subset.description
		print "subset.columns:"
		print subset.columns
		print "self.errors:"
		for error in subset.errors:
				print error

	print "\n\nTesting the MinimalSubsetFactory"
	
	test_name = subsets.keys()[0]
	print test_name
	
