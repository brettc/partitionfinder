from os import path
from partition_finder.MinimalSubsets import process_partition_file
from string import strip
from partition_finder.process_numerical_description import process_numerical_description
from partition_finder.Subset import SubsetFactory
from partition_finder.Scheme import Scheme
from partition_finder.submodels import get_submodels

curdir = path.dirname(path.abspath(__file__))
subset_factory = SubsetFactory()


partitions_filename = "%s/testfiles/test_partitions" %(curdir)
alignment_filename = "%s/testfiles/test.fas" %(curdir)
schemes_filename = "%s/testfiles/test_schemes2" %(curdir)

all_minimal_subsets = process_partition_file(partitions_filename)

N = len(all_minimal_subsets)
results2 = {}
all_submodels = get_submodels(N)
for model in all_submodels:
	m = ''
	for thing in model:
		m = ''.join([m,str(thing)])
	line = "%s = %s" %(m, m)
	print m
	s = Scheme(all_minimal_subsets, schemes_filename, 1, line, subset_factory, alignment_filename)

	results2[s.description] = s

print "\n\n\nResults\n"

scoresA = {}
scoresAc = {}
scoresB = {}
scoresL = {}

for k in results2:
	s = results2[k]
#	print s.name
#	print "lnL:", s.lnL
#	print "AIC:", s.AIC
	scoresA[s.AIC] = s.description
	scoresAc[s.AICc] = s.description
	scoresB[s.BIC] = s.description
	scoresL[s.lnL] = s.description
best_scoreA = min(scoresA.keys())
best_modelA = scoresA[best_scoreA]
best_scoreAc = min(scoresAc.keys())
best_modelAc = scoresAc[best_scoreAc]
best_scoreB = min(scoresB.keys())
best_modelB = scoresB[best_scoreB]
best_scoreL = max(scoresL.keys())
best_modelL = scoresL[best_scoreL]

print "The best AIC model is:"
print best_modelA
print "which has AIC score:", best_scoreA
print "\nThe best AICc model is:"
print best_modelAc
print "which has AICc score:", best_scoreAc
print "\nThe best BIC model is:"
print best_modelB
print "which has BIC score:", best_scoreB
print "\nThe best lnL model is:"
print best_modelL
print "which has lnL score:", best_scoreL


for i in range(N):
	n = i+1
	best = 9999999999999999
	for k in results2:
		s = results2[k]
		if s.n_subsets==n:
			if s.BIC<best: best=s.BIC
	print n, best
	
