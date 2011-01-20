class Partitioning_Scheme(object):
	"""This object is a way of lumping together the minimal_subsets into larger subsets
	   this is the essence of what PartitionFinder is trying to do: find the best lumping, or partitioning scheme.
	"""
	def __init__(self, method, line_number, line):
		self.method  		= method		#the method used to generate this scheme - either "file" "all" or "search"
		self.line_number	= line_number	#the position that this scheme was in that file (i.e. 0 would be first in the file)
		self.line			= line			#a string of the line of the schemes file e.g. "[pos1,pos2] [pos3]"
		self.name 			= None			#the name of the file as listed in the partitions file supplied by the user
		self.description	= None			#the text description of this partition, as supplied by the user in the file
		self.from_file		= None			#if method == file, then this records the input file
		self.n_subsets		= None			#the number of subsets in this scheme
		self.subsets		= None			#a list of the subsets that make up this scheme
		self.lnL			= None			#the log-likelihood of this scheme
		self.params			= None			#the number of parameters of this scheme
		self.AIC			= None			#the AIC score of this scheme
		self.AICc			= None			#the AICc score of this scheme
		self.BIC			= None			#the BIC score of this scheme