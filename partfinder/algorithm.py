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

from math import sqrt
from itertools import izip


def k_subsets_i(n, k):
    '''
    http://code.activestate.com/recipes/500268-all-k-subsets-from-an-n-set/
    Yield each subset of size k from the set of intergers 0 .. n - 1
    n -- an integer > 0
    k -- an integer > 0
    '''
    # Validate args
    if n < 0:
        raise ValueError('n must be > 0, got n=%d' % n)
    if k < 0:
        raise ValueError('k must be > 0, got k=%d' % k)
    # check base cases
    if k == 0 or n < k:
        yield set()
    elif n == k:
        yield set(range(n))

    else:
        # Use recursive formula based on binomial coeffecients:
        # choose(n, k) = choose(n - 1, k - 1) + choose(n - 1, k)
        for s in k_subsets_i(n - 1, k - 1):
            s.add(n - 1)
            yield s
        for s in k_subsets_i(n - 1, k):
            yield s


def k_subsets(s, k):
    '''
    http://code.activestate.com/recipes/500268-all-k-subsets-from-an-n-set/
    Yield all subsets of size k from set (or list) s
    s -- a set or list (any iterable will suffice)
    k -- an integer > 0
    '''
    s = list(s)
    n = len(s)
    for k_set in k_subsets_i(n, k):
        yield set([s[i] for i in k_set])


def lumpings(scheme):
    """
    generate all possible lumpings of a given scheme, where a lumping involves
    joining two partitions together scheme has to be a list of digits
    """
    # Get the numbers involved in the scheme
    nums = set(scheme)
    lumpings = []
    for sub in k_subsets(nums, 2):
        lump = list(scheme)
        sub = list(sub)
        sub.sort()
        # Now replace all the instance of one number in lump with the other in
        # sub
        while lump.count(sub[1]) > 0:
            lump[lump.index(sub[1])] = sub[0]
        lumpings.append(lump)

    return lumpings


def euclidean_distance(x, y):
    sum = 0
    for xval, yval in izip(x, y):
        sum += (xval - yval) ** 2
    dist = sqrt(sum)
    return dist


def levels_to_scheme(levels, namedict):
    """
    take the return from Cluster.getlevel
    and return it as a list of partition names description
    """

    levels = str(levels)

    for key in namedict.keys():
        old = str(namedict[key])
        new = '"%s"' % key
        levels = levels.replace(old, new)

    levels = eval(levels)
    return levels
