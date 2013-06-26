"""Operation on subsets: Joining / Splitting etc

    We keep some of the more arcane subset manipulations and operations in this
    file.

"""

import hashlib
import cPickle as pickle
import subset

def subset_unique_name(subset):
    """Return a unique string based on the subsets columns (which are unique)"""

    # Use pickle to give us a string from the object
    pickled_columns = pickle.dumps(subset.column_set, -1)

    # Now get an md5 hash from this. There is some vanishingly small chance that
    # we'll get the same thing. Google "MD5 Hash Collision"
    return hashlib.md5(pickled_columns).hexdigest()


def merge_subsets(subset_list):
    """Take a set of subsets and merge them together"""

    columns = set()

    # We just need the columns
    for sub in subset_list:
        columns |= sub.column_set

    return subset.Subset(sub.cfg, columns)


# TODO: FIX THESE
def has_overlap(subset_list):
    return False


def has_missing(subset_list):
    return False

