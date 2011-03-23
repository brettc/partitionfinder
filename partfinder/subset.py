import logging
log = logging.getLogger("subset")
import os
import weakref
import cPickle as pickle

import alignment
import phyml_models

class SubsetError(Exception):
    pass

class Subset(object):
    """A Subset of Partitions
    """

    # TODO: changes this to AllSubsets?
    _cache = weakref.WeakValueDictionary()

    # CLEVER BIT:
    # Return the SAME subset if the partitions are identical
    # This is basically a pythonized factory
    # See here: http://codesnipers.com/?q=python-flyweights
    def __new__(cls, *parts):
        cacheid = frozenset(parts)
        obj = Subset._cache.get(cacheid, None)
        if not obj:
            obj = object.__new__(cls)
            Subset._cache[cacheid] = obj

            # Error checking....
            tempparts = set()
            for p in parts:
                if p.partition_set is None:
                    log.error("You cannot add a Partition to a Subset until "
                              "the Partition belongs to a PartitionSet")
                    raise SubsetError

                if p in tempparts:
                    log.error("%s is duplicated in a Subset", p)
                    raise SubsetError

                tempparts.add(p)

            obj.partitions = cacheid
            
            # Append all of the columns in the partition -- I think these are
            # useful...
            obj.columns = []
            obj.columnset = set()
            for p in parts:
                obj.columns += p.columns
                obj.columnset |= p.columnset
            obj.columns.sort()

            obj.results = {}
            obj.best_aic = None
            obj.best_model = None
            log.debug("Created %s", obj)
        # else:
            # log.debug("Reused %s", obj)

        return obj

    # def __init__(self, *parts):
        # Everything is relegated to above...

    def __str__(self):
        return "Subset(%s)" % ", ".join([str(p) for p in self.partitions])

    @property
    def name(self):
        s = sorted([p.name for p in self.partitions])
        return '-'.join(s)

    def __iter__(self):
        return iter(self.partitions)

    def add_model_result(self, model, result):
        result.model = model
        result.params = phyml_models.get_num_params(model)
        result.aic = 2 * (result.params - result.lnl)
        if model in self.results:
            log.error("Can't add model result %s, it already exists in %s",
                    model, self)
        self.results[model] = result

        if self.best_aic is None or result.aic < self.best_aic:
            self.best_lnl = result.lnl
            self.best_aic = result.aic
            self.best_model = result.model
            self.best_params = result.params

    _template = "%-15s | %-15s | %-15s\n"
    def write_summary(self, path):
        # Sort everything
        model_results = [(r.aic, r) for r in self.results.values()]
        model_results.sort()
        f = open(path, 'w')
        f.write("Results for %s\n\n" % self)
        f.write(Subset._template % ("Model", "lNL", "AIC"))
        for aic, r in model_results:
            f.write(Subset._template % (r.model, r.lnl, r.aic))

    def write_binary_summary(self, path):
        """Write out the results we've collected to a binary file"""
        f = open(path, 'rb')
        pickle.dump(self.results, f, -1)

    def read_binary_summary(self, path):
        if not os.path.exists(path):
            return False

        f = open(path, 'wb')
        self.results = pickle.load(path)


if __name__ == '__main__':
    import logging
    import tempfile
    logging.basicConfig(level=logging.DEBUG)
    import config
    from partition import Partition

    tmp = tempfile.mkdtemp()
    config.initialise(tmp, True)

    pa = Partition('a', (1, 10, 3))
    pb = Partition('b', (2, 10, 3))
    pc = Partition('c', (3, 10, 3))

    s1 = Subset(pa, pb)
    s2 = Subset(pa, pb)
    s3 = Subset(pa, pc)
    s4 = Subset(pa, pb)
    print s1 is s2
    print s1
    # s2 = Subset(pa, pb, pc)
    # s3 = Subset(pc)

    # print s1.name
    # print s2.name


