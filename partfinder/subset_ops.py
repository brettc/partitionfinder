"""Operation on subsets: Joining / Splitting etc

    We keep some of the more arcane subset manipulations and operations in this
    file.

"""

import hashlib
import cPickle as pickle

def subset_unique_name(subset):
    """Return a unique string based on the subsets columns (which are unique)"""

    # Use pickle to give us a string from the object
    pickled_columns = pickle.dumps(subset.column_set, -1)

    # Now get an md5 hash from this. There is some vanishingly small chance that
    # we'll get the same thing. Google "MD5 Hash Collision"
    return hashlib.md5(pickled_columns).hexdigest()



