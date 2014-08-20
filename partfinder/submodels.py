# Copyright (C) 2012 Robert Lanfear and Brett Calcott
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details. You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# PartitionFinder also includes the PhyML program, the RAxML program, and the
# PyParsing library, all of which are protected by their own licenses and
# conditions, using PartitionFinder implies that you agree with those licences
# and conditions as well.

import logtools
log = logtools.get_logger()

import math

def submodel_generator(result, pat, current, maxn):
    """ result is a list to append to
        pat is the current pattern (starts as empty list)
        current is the current number of the pattern
        maxn is the number of items in the pattern
    """
    if pat:
        curmax = max(pat)
    else: 
        curmax = 0
    for i in range(current):
        if i-1 <= curmax:
            newpat = pat[:]
            newpat.append(i)
            if current == maxn:
                result.append(newpat)
            else:
                submodel_generator(result, newpat, current+1, maxn)

def submodel_iterator(pat, current, maxn):
    """same as generator but yields instead"""
    if pat:
        curmax = max(pat)
    else: 
        curmax = 0
    for i in range(current):
        if i-1 <= curmax:
            newpat = pat[:]
            newpat.append(i)
            if current == maxn:
                yield newpat
            else:
                for b in submodel_iterator(newpat, current+1, maxn):
                    yield b

def a_choose_b(n,k):
    return reduce(lambda a,b: a*(n-b)/(b+1),xrange(k),1)

def count_relaxed_clustering_subsets(N, cluster_percent, cluster_max):
    #startscheme    
    start_scheme = N
    #firstbatch is just cluster_percent of N choose 2
    step_1 = int(math.ceil(a_choose_b(N, 2)*cluster_percent*0.01))
    if cluster_max != None and step_1>cluster_max:
        step_1 = cluster_max    
    previous = step_1
    cumsum = start_scheme+step_1

    #now for the rest
    for i in reversed(xrange(N)):

        # once we get to the all combined scheme we can stop  
        if i == 1:
            break
        num_new_schemes = int(math.ceil((a_choose_b(i, 2))*cluster_percent*0.01))
        # but those new schemes include a lot we will have already analysed
        # so we want to subtract that many. We could have already seen up to i-1 choose 2
        # the worst case is that the scheme we chose knocked out the maximum number of 
        # previously analysed schemes, which is just 2(i)-1, so:
        worst_case = 2*i - 1
        num_already_analysed = previous - worst_case
        # now we transfer over the 'previous' for the next round of the loop
        previous = num_new_schemes
        # now we calculate the final number of new schemes
        num_new_schemes -= num_already_analysed
        if cluster_max != None and num_new_schemes>cluster_max:
            num_new_schemes = cluster_max
        cumsum += num_new_schemes

    return cumsum

def count_relaxed_clustering_schemes(N, cluster_percent, cluster_max):
    #startscheme    
    start_scheme = 1
    #firstbatch is just cluster_percent of N choose 2
    step_1 = int(math.ceil(a_choose_b(N, 2)*cluster_percent*0.01))

    if cluster_max != None and step_1>cluster_max:
        step_1 = cluster_max

    cumsum = start_scheme+step_1
    #now for the rest
    for i in reversed(xrange(N)):
        # each subsequent step is cluster_percent of i choose 2  
        if i == 1:
            break
        num_new_schemes = int(math.ceil((a_choose_b(i, 2))*cluster_percent*0.01))
        if cluster_max != None and num_new_schemes>cluster_max:
            num_new_schemes = cluster_max
        cumsum += num_new_schemes
    return cumsum

def count_greedy_schemes(N):
    """oeis.org reveals this is 1+(N*(N+1)*(N-1))/6"""
    count = 1+(N*(N+1)*(N-1))/6
    return count

def count_greedy_subsets(N):
    """oeis.org says thes are Central polygonal numbers: n^2 - n + 1. """
    count = (N*N) - N + 1
    return count

def bell_numbers(N):
    ## Return the bell number for N subsets
    # script modified from Wikipedia: http://en.wikipedia.org/wiki/Bell_number
    N += 1## Bell numbers are indexed from zero
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
    """
    log.debug("Generating submodels for %s partitions", N)
    result = []
    submodel_generator(result, [], 1, N)
    log.debug("Resulting number of partitions is %s", len(result))
    return result

def count_all_schemes(N):
    """Count the number of submodels we"ve got"""
    count = bell_numbers(N)
    return count

def count_all_subsets(N):
    """Count the number of subses we'll have to look at given a certain number of starting partitions"""
    count = (2**N) - 1
    return count

