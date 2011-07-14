#Copyright (C) 2011 Robert Lanfear and Brett Calcott
#
#This program is free software: you can redistribute it and/or modify it
#under the terms of the GNU General Public License as published by the
#Free Software Foundation, either version 3 of the License, or (at your
#option) any later version.
#
#This program is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#General Public License for more details. You should have received a copy
#of the GNU General Public License along with this program.  If not, see
#<http://www.gnu.org/licenses/>. PartitionFinder also includes the PhyML
#program and the PyParsing library both of which are protected by their
#own licenses and conditions, using PartitionFinder implies that you
#agree with those licences and conditions as well.

import logging
log = logging.getLogger("submodels")
import algorithm

def submodel_generator(result, pat, current, maxn, countonly=False):
    if pat:
        curmax = max(pat)
    else: 
        curmax = 0
    for i in range(current):
        if i-1 <= curmax:
            newpat = pat[:]
            newpat.append(i)
            if current == maxn:
                if not countonly:
                    result.append(newpat)
                else:
                    result[0] += 1
            else:
                submodel_generator(result, newpat, current+1, maxn, countonly)


def a_choose_b(n,k):
    return reduce(lambda a,b: a*(n-b)/(b+1),xrange(k),1)

def count_greedy_schemes(N):
    '''oeis.org reveals this is 1+(N*(N+1)*(N-1))/6'''
    count = 1+(N*(N+1)*(N-1))/6
    return count

def count_greedy_subsets(N):
    '''oeis.org says thes are Central polygonal numbers: n^2 - n + 1. '''
    count = (N*N) - N + 1
    return count


def bell_numbers(N):
    ## Return the bell number for N subsets
    # script modified from Wikipedia: http://en.wikipedia.org/wiki/Bell_number
    N = N+1                          ## Bell numbers are indexed from zero 
    t = [[1]]                        ## Initialize the triangle as a two-dimensional array
    c = 1                            ## Bell numbers count
    while c <= N:
        if c >= N:
            return t[-1][0]          ## Yield the Bell number of the previous row
        row = [t[-1][-1]]            ## Initialize a new row
        for b in t[-1]:
            row.append(row[-1] + b)  ## Populate the new row
        c += 1                       ## We have found another Bell number
        t.append(row)                ## Append the row to the triangle
 

def get_submodels(N):
    """Return all the submodels

    Let's do a test case:
    >>> submodels = get_submodels(4)
    >>> for s in submodels:
    ...     print s
    [0, 0, 0, 0]
    [0, 0, 0, 1]
    [0, 0, 1, 0]
    [0, 0, 1, 1]
    [0, 0, 1, 2]
    [0, 1, 0, 0]
    [0, 1, 0, 1]
    [0, 1, 0, 2]
    [0, 1, 1, 0]
    [0, 1, 1, 1]
    [0, 1, 1, 2]
    [0, 1, 2, 0]
    [0, 1, 2, 1]
    [0, 1, 2, 2]
    [0, 1, 2, 3]
    """
    log.debug("Generating submodels for %s partitions", N)
    result = []
    submodel_generator(result, [], 1, N)
    log.debug("Resulting number of partitions is %s", len(result))
    return result

def count_all_schemes(N):
    """Count the number of submodels we've got
    These are the right numbers...
    >>> print count_all_schemes(1)
    1
    >>> print count_all_schemes(5)
    52
    >>> print count_all_schemes(10)
    115975
    """
    count = bell_numbers(N)
    return count

def count_all_subsets(N):
    '''Count the number of subses we'll have to look at given a certain number of starting partitions'''
    count = (2**N) - 1
    return count

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
    print "Getting All Submodels for N=5"
    result = get_submodels(5)
    for a in result:
        print a
    
    print "A table of number of partitions versus number of schemes and number of subsets for greedy analyses"
    print "Parts\tGreedySchemes\tGreedySubsets\tAllSchemes\tAllSubsets"
 
    for i in range(1,201):
        print "%d, %g, %g, %g, %g" %(i, count_greedy_schemes(i), count_greedy_subsets(i), count_all_schemes(i), count_all_subsets(i))