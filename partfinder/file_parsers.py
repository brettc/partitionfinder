from os import path	
from string import strip

def parse_RaxML_partitions_file(partitions_filename):
	#Parse a partitions file, to get a dictionary keyed by name,
	#the entries in this dictionary have: start, stop, and step values to get data from alignments
	print "Reading partitions data from %s" %(partitions_filename) 
	partitions = open(partitions_filename)
	partitions_dict = {}
	partitions_list = [] #an ordered list of the keys, useful for some things
	for line in partitions:
		words = line.split("=")
		name = strip(words[0])
		description = strip(words[1])
		columns = description.split("\\") #split by the \ character, if it's there
		start_stop = columns[0].split("-") #extract the start and stop values
		start = int(start_stop[0])-1 # remember we're zero indexed here...
		stop = int(start_stop[1])
		step = 1
		if len(columns)==2:
			step = int(columns[1])
		partitions_dict[name] = start, stop, step #this lets us slice the sequence appropriately
		partitions_list.append(name)
	print "The following partitions have been defined"
	for partition in partitions_list:
		print partition,
		print "start: %d\t\tstop: %d\t\tgapsize: %d" %(partitions_dict[partition][0], partitions_dict[partition][1], partitions_dict[partition][2])
	return partitions_dict, partitions_list

def parse_combinations_file(combinations_filename):
	#Parse a combinations file, to get a dictionary keyed by name,
	#the entries in this dictionary are list representations of the scheme
	combinations = open(combinations_filename)
	combinations_dict = {}
	for line in combinations:
		words = line.split("=")
		name = strip(words[0])
		scheme = list(strip(words[1]))
		combinations_dict[name] = scheme #this lets us slice the sequence appropriately
	return combinations_dict

if __name__ == "__main__":

	curdir = path.dirname(path.abspath(__file__))
	test_partitionfile 		= "%s/testfiles/test_partitions" %(curdir)
	test_combinationfile 	= "%s/testfiles/test_combinations" %(curdir)

	parse_RaxML_partitions_file(test_partitionfile)
