import logging
log = logging.getLogger("scheme")

from math import log as logarithm

from util import PartitionFinderError
class SchemeError(PartitionFinderError):
    pass

class Scheme(object):
    def __init__(self, cfg, name, *subsets):
        """A set of subsets of partitions"""
        self.name = name
        self.subsets = set()

        # This one is a set of frozensets of partitions...
        part_subsets = set()

        # This is really long-winded, but it is mainly for error-checking
        partitions = set()
        duplicates = []
        for s in subsets:
            for p in s:
                if p in partitions:
                    # This is an error -- we'll collect them up
                    duplicates.append(str(p))
                else:
                    partitions.add(p)
            self.subsets.add(s)
            part_subsets.add(s.partitions)

        self.part_subsets = frozenset(part_subsets)

        # Report the errors 
        if duplicates:
            log.error("Scheme '%s' contains duplicate partitions: %s", 
                      name, ', '.join(duplicates))
            raise SchemeError

        # Hm. It seems this is the only way to get just one item out of a set
        # as pop would remove one...
        pset = cfg.partitions

        # Do a set-difference to see what is missing...
        missing = pset.partitions - partitions
        if missing:
            log.error("Scheme '%s' is missing partitions: %s", 
                      name, ', '.join([str(p) for p in missing]))
            raise SchemeError

        # This locks down whether new partitions can be created.
        if not cfg.partitions.finalised:
            cfg.partitions.finalise()
        
        # Now add to all_schemes -- more possibility of errors, see below
        cfg.schemes.add_scheme(self)
        log.debug("Created %s", self)

    def __iter__(self):
        return iter(self.subsets)

    def __str__(self):
        ss = ', '.join([str(s) for s in self.subsets])
        return "Scheme(%s, %s)" % (self.name, ss)

    def assemble_results(self, nseq, branchlengths):
		#calculate AIC, BIC, AICc for each scheme.
		#how you do this depends on whether brlens are linked or not.
        self.nsubs = len(self.subsets) #number of subsets
        sum_subset_k = sum([s.best_params for s in self]) #sum of number of parameters in the best model of each subset
        
        if branchlengths == 'linked': #linked brlens - only one extra parameter per subset
            self.sum_k = sum_subset_k + (self.nsubs-1) + ((2*nseq)-3) #number of parameters in a scheme
        elif branchlengths == 'unlinked': #unlinked brlens - every subset has its own set of brlens
            self.sum_k = sum_subset_k + (self.nsubs*((2*nseq)-3)) #number of parameters in a scheme

        else:
            # WTF?
            log.error("Unknown option for branchlengths: %s", branchlengths)
            raise AnalysisError
        
        self.lnl = sum([s.best_lnl for s in self])
        self.nsites = sum([len(s.columnset) for s in self])

        K = float(self.sum_k)
        n = float(self.nsites)
        lnL = float(self.lnl)		
   
        self.aic  = (-2.0*lnL) + (2.0*K)
        self.bic  = (-2.0*lnL) + (K * logarithm(n))
        self.aicc = self.aic + (((2.0*K)*(K+1.0))/(n-K-1.0))


    _header_template = "%-15s: %s\n"
    _subset_template = "%-6s | %-10s | %-50s | %-40s\n"
    def write_summary(self, path, write_type = 'wb', extra_line=None):
        f = open(path, write_type)
        if extra_line:
		    f.write(extra_line)
        f.write(Scheme._header_template % ("Scheme Name", self.name))
        f.write(Scheme._header_template % ("Scheme lnL", self.lnl))
        f.write(Scheme._header_template % ("Scheme AIC", self.aic))
        f.write(Scheme._header_template % ("Scheme AICc", self.aicc))
        f.write(Scheme._header_template % ("Scheme BIC", self.bic))
        f.write(Scheme._header_template % ("Num params", self.sum_k))
        f.write(Scheme._header_template % ("Num sites", self.nsites))
        f.write(Scheme._header_template % ("Num subsets", self.nsubs))
        f.write("\n")
        f.write(Scheme._subset_template % (
            "Subset", "Best Model", "Subset Partitions",  "Alignment"))
        number = 1
        
        sorted_subsets = [sub for sub in self]
        sorted_subsets.sort(key=lambda sub: min(sub.columns), reverse=False)
                
        for sub in sorted_subsets:
            desc = []
            names= []
            for part in sub:
                desc.extend(part.description)
                names.append(part.name)
            names.sort()
            names = ', '.join(names)
			
            f.write(Scheme._subset_template % (
                number, sub.best_model, names, sub.alignment_path))
            number = number + 1

class SchemeSet(object):
    """All the schemes added, and also a list of all unique subsets"""
    def __init__(self):
        """A collection of schemes"""
        self.clear_schemes()

    def clear_schemes(self):
        self.schemes_by_name = {}
        self.schemes_by_subsets = {}

    def add_scheme(self, scheme):
        if scheme.name in self.schemes_by_name:
            log.error("Cannot add two schemes with same name: '%s'" %
                      scheme.name)
            raise SchemeError

        if scheme.part_subsets in self.schemes_by_subsets:
            existing_scheme = \
                    self.schemes_by_subsets[scheme.part_subsets]
            log.warning(
                "Scheme named %s being added is identical to existing %s",
                scheme.name, existing_scheme)
            # raise SchemeError

        self.schemes_by_name[scheme.name] = scheme
        self.schemes_by_subsets[scheme.part_subsets] = scheme

    def __len__(self):
        return len(self.schemes_by_name)
    # Easy iteration
    def __iter__(self):
        return iter(self.schemes_by_name.itervalues())

def create_scheme(cfg, scheme_name, scheme_description):
    """Generate a single scheme given a list of numbers e.g. [0,1,2,3,4,5,6,7]"""
    import subset
    import submodels
    
    partnum = len(cfg.partitions) #total number of partitions defined by user

    #check that the correct number of items are in the list
    if len(scheme_description)!=partnum:
        log.error("There's a problem with the description of scheme %s" % scheme_name)
        raise SchemeError

    # Now generate the pattern
    subs = {}
    # We use the numbers returned to group the different subsets
    for sub_index, grouping in enumerate(scheme_description):
        insub = subs.setdefault(grouping, [])
        insub.append(sub_index)
    # We now have what we need to create a subset. Each entry will have a
    # set of values which are the index for the partition
    created_subsets = []
    for sub_indexes in subs.values():
        sub = subset.Subset(*tuple([cfg.partitions[i] for i in sub_indexes]))
        created_subsets.append(sub)

    new_scheme = Scheme(cfg, str(scheme_name), *tuple(created_subsets))
		
    return new_scheme

def generate_all_schemes(cfg):
    """Convert the abstract schema given by the algorithm into subsets"""
    import subset
    import submodels

    log.info("Generating all possible schemes for the partitions...")

    # Make sure that no schemes have been defined!
    # if len(all_schemes) > 0:
        # log.error("Cannot generate schemes if some already exist!")
        # raise SchemeError
    
    partnum = len(cfg.partitions) #total number of partitions defined by user

    # Now generate the pattern for this many partitions
    all_schemes = submodels.get_submodels(partnum)
    scheme_name = 1
    scheme_list = []
    for scheme in all_schemes:
        subs = {}
        # We use the numbers returned to group the different subsets
        for sub_index, grouping in enumerate(scheme):
            insub = subs.setdefault(grouping, [])
            insub.append(sub_index)
        # We now have what we need to create a subset. Each entry will have a
        # set of values which are the index for the partition
        created_subsets = []
        for sub_indexes in subs.values():
            sub = subset.Subset(*tuple([cfg.partitions[i] for i in sub_indexes]))
            created_subsets.append(sub)

        scheme_list.append(
            Scheme(cfg, str(scheme_name), *tuple(created_subsets)))
        #log.info("Created scheme %d of %d" %(scheme_name, len(all_schemes)))

        scheme_name += 1
		
    return scheme_list

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.DEBUG)
    from partition import Partition
    from subset import Subset

    pa = Partition('a', (1, 10, 3))
    pb = Partition('b', (2, 10, 3))
    pc = Partition('c', (3, 10, 3))
    pd = Partition('d', (11, 20))
    # s = Scheme('x', Subset(pa, pc), Subset(pb))

    generate_all_schemes()
    # This should give us an error!
    # s = Scheme('y', Subset(pa, pc), Subset(pb))
    
