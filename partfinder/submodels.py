import logging
log = logging.getLogger("submodels")

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

def count_submodels(N):
    """Count the number of submodels we've got

    These are the right numbers...
    >>> print count_submodels(1)
    1
    >>> print count_submodels(5)
    52
    >>> print count_submodels(10)
    115975
    """
    result = [0]
    submodel_generator(result, [], 1, N, countonly=True)
    return result[0]

if __name__ == "__main__":
    import doctest
    doctest.testmod()
    
