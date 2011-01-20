def submodel_generator(result, pat, current, maxn):
	if pat:
		curmax = max(pat)
	else: 
		curmax = 0
	for i in range(current):
		if i-1<=curmax:
			newpat = pat[:]
			newpat.append(i)
			if current == maxn:
				result.append(newpat)
			else:
				submodel_generator(result, newpat, current+1, maxn)

def get_submodels(N):
	result = []
	submodel_generator(result, [], 1, N)
	return result


def submodel_counter(result, pat, current, maxn):
	if pat:
		curmax = max(pat)
	else: 
		curmax = 0
	for i in range(current):
		if i-1<=curmax:
			newpat = pat[:]
			newpat.append(i)
			if current == maxn:
				result[0] = result[0]+1
			else:
				submodel_counter(result, newpat, current+1, maxn)

def count_submodels(N):
	result = [0]
	submodel_counter(result, [], 1, N)
	return result[0]


if __name__ == "__main__":
	
	N = range(10)
	
	for i in N:
		n = i+1
	
		print "Number of Submodels for N = %d\t" %(n),
		print count_submodels(n)
		
	print "\nThe 15 submodels for N=4 look like this:"
	
	submodels = get_submodels(4)
	
	for thing in submodels:
		print thing