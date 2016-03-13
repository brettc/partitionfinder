# array lookups

import numpy as np
import random
import sklearn.metrics.pairwise
import scipy.spatial.distance
import heapq
from bisect    import insort
from itertools import islice
from scipy.spatial import KDTree



N = 1000
r = np.array([random.randrange(1, 1000) for _ in range(0, 10000)])
c = r[:, None]

dists = scipy.spatial.distance.pdist(c, 'cityblock')

def method1(dists, N):
    d_sq = scipy.spatial.distance.squareform(dists)
    unique_dists = np.unique(dists)[:N]
    l = []
    for d in unique_dists:
        locations = np.where(d_sq==d)

        for i in range(len(locations[0])):
            pair = [locations[0][i], locations[1][i]]
            l.append(pair)
    return l



def condensed_to_square_index(n, c):
    # converts an index in a condensed array to the 
    # pair of observations it represents
    # modified from here: http://stackoverflow.com/questions/5323818/condensed-matrix-function-to-find-pairs
    ti = np.triu_indices(n, 1)
    return ti[0][c]+ 1, ti[1][c]+ 1

def condensed_to_square_index2(ti, c):
    return ti[0][c]+ 1, ti[1][c]+ 1


def method2(dists, N):
    closest = dists.argsort()[:N]
    n = np.ceil(np.sqrt(2* len(dists)))
    r = []
    for i in closest:
        pair = condensed_to_square_index(n, i)
        r.append(pair)
    return r


def method3(dists, N):
    closest = dists.argsort()[:N]
    n = np.ceil(np.sqrt(2* len(dists)))
    ti = np.triu_indices(n, 1)
    r = np.vstack(ti)[:, closest] + 1
    return r

def method4(dists, N):
    closest = dists.argsort()[:N]
    n = np.ceil(np.sqrt(2* len(dists)))
    ti = np.triu_indices(n, 1)
    r  = zip(ti[0][closest] + 1, ti[1][closest] + 1)
    return r

def method5(dists, N):
    closest = dists.argsort()[:N]
    n = np.ceil(np.sqrt(2* len(dists)))
    ti = np.triu_indices(n, 1)
    r = np.column_stack((np.take(ti[0], closest), np.take(ti[1], closest))) + 1
    return r


def method6(dists, N):
    closest = dists.argsort()[:N]
    n = np.ceil(np.sqrt(2* len(dists)))
    ti = np.triu_indices(n, 1)
    r = []
    for i in closest:
        pair = condensed_to_square_index2(ti, i)
        r.append(pair)
    return r

def argsmallest_n(a, n):
    ret = np.argpartition(a, n)[:n]
    b = np.take(a, ret)
    return np.take(ret, np.argsort(b))

# NB method7 requires NumPy 1.8
def method7(dists, N):
    closest = argsmallest_n(dists, N)
    n = np.ceil(np.sqrt(2* len(dists)))
    tu = np.triu_indices(n, 1)
    pairs = np.column_stack((np.take(tu[0], closest), np.take(tu[1], closest))) + 1    
    return pairs

def method8(dists, N):
    closest = heapq.nsmallest(N, dists)
    n = np.ceil(np.sqrt(2* len(dists)))
    ti = np.triu_indices(n, 1)
    r  = zip(ti[0][closest] + 1, ti[1][closest] + 1)
    return r


def method9(dists, N, insort=insort):
    # http://stackoverflow.com/questions/350519/getting-the-lesser-n-elements-of-a-list-in-python
    it   = iter(dists)
    closest = sorted(islice(it, N))
    for el in it:
        if el <= closest[-1]: #NOTE: equal sign is to preserve duplicates
            insort(closest, el)
            closest.pop()

    n = np.ceil(np.sqrt(2* len(dists)))
    ti = np.triu_indices(n, 1)
    r  = zip(ti[0][closest] + 1, ti[1][closest] + 1)
    return r

def method10(dists, N):
    # http://stackoverflow.com/questions/350519/getting-the-lesser-n-elements-of-a-list-in-python    
    mins = dists[:N]
    mins.sort()
    for i in dists[N:]:
        if i <= mins[-1]: 
            np.append(mins, i)
            mins.sort()
            mins = mins[:N]
    closest = [np.where(dists==x) for x in mins]
    n = np.ceil(np.sqrt(2* len(dists)))
    ti = np.triu_indices(n, 1)
    r  = zip(ti[0][closest] + 1, ti[1][closest] + 1)
    return r


def method11(dists, N):
    # http://stackoverflow.com/questions/350519/getting-the-lesser-n-elements-of-a-list-in-python    
    mins = list(dists[:N])
    mins.sort()
    for i in dists[N:]:
        if i <= mins[-1]: 
            mins.append(i)
            insort(mins, i)
            mins = mins[:N]

    closest = [np.where(dists==x) for x in mins]
    n = np.ceil(np.sqrt(2* len(dists)))
    ti = np.triu_indices(n, 1)
    r  = zip(ti[0][closest] + 1, ti[1][closest] + 1)
    return r


def method12(dists, N):
    tree = KDTree(dists, leafsize=dists.shape[0]+1)
    distances, ndx = tree.query(0, k=N)

def method4_2(dists, N):
    closest = dists.argsort()[:N]

# timing in ipython, best is method 4 when the array is 1000 by 1000, and N = 1000

# best method is 9 when array gets > ~6000 by 6000, and 9 is only marginally slower for smaller arrays

# but note that either using the bottleneck module, or method7 (both of which use 
# efficient sorting) would be a lot quicker.


# NB arrays of 10k by 10k still take ~20s with method 4. This is not great.









