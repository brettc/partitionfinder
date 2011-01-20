from os import path
from MinimalSubsets import process_partition_file
from string import strip
from process_numerical_description import process_numerical_description
from Subset import SubsetFactory

class Scheme(object):
	"""This object is a way of lumping together the minimal_subsets into larger subsets
	   this is the essence of what PartitionFinder is trying to do: find the best lumping, or partitioning scheme.
	"""
	def __init__(self, allminimalsubsets, from_file, line_number, line, subset_factory, alignment_filename):
		self.allminimalsubsets	= allminimalsubsets	#a dict of all of the MinimalSubset objects, which identify the partitions
		self.from_file  		= from_file			#the file that this scheme was loaded from
		self.line_number		= line_number		#i.e., this was the nth scheme created (e.g. useful for remembering where in the input file the scheme was)
		self.line				= line				#a string of the line of the partition file e.g. "pos1 = 1-500\3;"
		self.subset_factory		= subset_factory
		self.alignment_filename = alignment_filename
		self.name 				= None				#the name of the subset as listed in the partitions file supplied by the user
		self.description		= None				#the text description of this partition, as supplied by the user in the file
		self.subset_list		= None				#a list of the subsets that make up this scheme, each of these is a string
		self.subset_objs		= None				#a list of subset objects from the subset_list
		self.n_subsets			= None				#the number of subsets in this scheme
		self.lnL				= None				#the log-likelihood of this scheme
		self.n_params			= None				#the number of parameters of this scheme
		self.AIC				= None				#the AIC score of this scheme
		self.AICc				= None				#the AICc score of this scheme
		self.BIC				= None				#the BIC score of this scheme
		self.process_scheme_line()
		self.process_description()
		self.get_subset_objs()
		self.process_subsets()
		
	def process_scheme_line(self):
		"""turn a line from a scheme into a name and a description"""
		words = self.line.split("=")
		self.name = strip(words[0]) #that bit's easy
		description = strip(words[1]) #the description might have whitespace on the ends
		description = description.strip(";") #some people put semicolons too
		description = description.strip() # and there may be whitespace between the semicolon and the rest of it  
		self.description = description

	def process_description(self):
		"""turn a description into a list of subsets, the entries should all match up to the partitions in the partitions file"""
		if self.description.isdigit(): #it's the numbery kind of description
			self.description = process_numerical_description(self.description, self.allminimalsubsets, debug=0)
		sublist = self.description.split('(')
		self.subset_list = []
		for sub in sublist:
			if len(sub)>0:
				sub = sub.strip()
				sub = sub.strip(',') #no commas at the start or end please
				sub = sub.strip()
				sub = ''.join(['(', sub]) #put that bracket back on
				self.subset_list.append(sub)
		self.n_subsets = len(self.subset_list)
	
	def get_subset_objs(self):
		"""get a list of subset objects"""
		self.subset_objs = []
		for description in self.subset_list:
			self.subset_objs.append(self.subset_factory.get_subset(self.allminimalsubsets, self.alignment_filename, description))

	def process_subsets(self):
		"""process_subsets to get likelihood, and various other metrics for a given Scheme"""
		sum_lnL = 0
		sum_AIC = 0
		sum_AICc = 0
		sum_BIC = 0
		for subset in self.subset_objs:
			sum_lnL = sum_lnL + subset.lnL[subset.lnL_best]
			sum_AIC = sum_AIC + subset.AIC[subset.AIC_best]
			sum_AICc = sum_AICc + subset.AICc[subset.AICc_best]
			sum_BIC = sum_BIC + subset.BIC[subset.BIC_best]
		self.lnL = sum_lnL
		self.AIC = sum_AIC
		self.BIC = sum_BIC
		self.AICc = sum_AICc
			

if __name__ == '__main__':
	curdir = path.dirname(path.abspath(__file__))
	partitions_filename = "%s/testfiles/test_partitions" %(curdir)
	alignment_filename = "%s/testfiles/test.fas" %(curdir)
	schemes_filename = "%s/testfiles/test_schemes1" %(curdir)

	subset_factory = SubsetFactory()
	 

	all_minimal_subsets = process_partition_file(partitions_filename)

	line = "12_3_by_gene 	= 		001223445"
	
	s = Scheme(all_minimal_subsets, schemes_filename, 0, line, subset_factory, alignment_filename)

	print "s.from_file:", s.from_file
	print "s.line_number:", s.line_number
	print "s.line:", s.line
	print "s.name:", s.name
	print "s.description:", s.description
	print "s.subsetlist:\n", 
	for thing in s.subset_list:
		print thing
	print "s.n_subsets:", s.n_subsets
	print "s.lnL:", s.lnL
	print "s.AIC:", s.AIC
