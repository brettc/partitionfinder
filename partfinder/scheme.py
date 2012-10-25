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
log = logging.getLogger("scheme")
from cluster import *
from algorithm import *

from math import log as logarithm

from util import PartitionFinderError
class SchemeError(PartitionFinderError):
    pass

class SchemeResult(object):
    def __init__(self, sch, nseq, branchlengths):
        self.scheme = sch

		# Calculate AIC, BIC, AICc for each scheme.
		# How you do this depends on whether brlens are linked or not.
        self.nsubs = len(sch.subsets) #number of subsets
        sum_subset_k = sum([s.best_params for s in sch]) #sum of number of parameters in the best model of each subset

        log.debug("Calculating number of parameters in scheme:")
        log.debug("Total parameters from subset models: %d" %(sum_subset_k))
        
        if branchlengths == 'linked': #linked brlens - only one extra parameter per subset
            self.sum_k = sum_subset_k + (self.nsubs-1) + ((2*nseq)-3) #number of parameters in a scheme
            log.debug("Total parameters from brlens: %d" %((2*nseq) -3))
            log.debug("Parameters from subset multipliers: %d" %(self.nsubs-1))        

        elif branchlengths == 'unlinked': #unlinked brlens - every subset has its own set of brlens
            self.sum_k = sum_subset_k + (self.nsubs*((2*nseq)-3)) #number of parameters in a scheme
            log.debug("Total parameters from brlens: %d" %((2*nseq) -3)*self.nsubs)

        else:
            # WTF?
            log.error("Unknown option for branchlengths: %s", branchlengths)
            raise AnalysisError
        
        log.debug("Grand total parameters: %d" %(self.sum_k))
        
        self.lnl = sum([s.best_lnl for s in sch])
        self.nsites = sum([len(s.columnset) for s in sch])

        K = float(self.sum_k)
        n = float(self.nsites)
        lnL = float(self.lnl)		

        log.debug("n: %d\tK: %d" %(n, K))
   
        #here we put in a catch for small subsets, where n<K+2
        #if this happens, the AICc actually starts rewarding very small datasets, which is wrong
        #a simple but crude catch for this is just to never allow n to go below k+2
        if n<(K+2): 
            log.warning("Scheme '%s' has a very small"
                        " number of sites (%d) compared to the number of parameters"
                        " in the models that make up the subsets"
                        " This may give misleading AICc results, so please check carefully"
                        " if you are using the AICc for your analyses."
                        " The results for this scheme are in the following file:" 
                        " /analysis/schemes/%s.txt\n" % (sch.name, n, sch.name))
            n = K+2 

        self.aic  = (-2.0*lnL) + (2.0*K)
        self.bic  = (-2.0*lnL) + (K * logarithm(n))
        self.aicc = (-2.0*lnL) + ((2.0*K)*(n/(n-K-1.0)))
        
    def __repr__(self):
        return "SchemeResult<aic:%f, aicc:%f, bic:%f>" % (self.aic, self.aicc,
                                                          self.bic)


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

    def get_clustering(self, cfg, method, scheme_name):
        """Ask a scheme which two subsets are most similar
        we return a scheme that has the closest subsets joined together
        """
        import subset
    
        #first we get a dictionary of subset parameters, keyed by subset name
        param_dict = {}
        for s in self.subsets:
            param_dict[s.full_name] = s.get_param_values()

        if method=='hierarchical':                
            #perform hierarchical clustering with means of the parameters
            cl_H = HierarchicalClustering(param_dict.values(), euclidean_distance)
            cl_H.setLinkageMethod("uclus")
            cl_H.cluster()
            
            #now we get all the levels and order them
            d = cl_H.data[0]
            levs = getLevels(d, [])
            levs.sort()
        
            #and now we use the first level to get the two most similar subsets
            newscheme_description = levels_to_scheme(cl_H.getlevel(levs[0]), param_dict)

        print "yo"
        print newscheme_description
        created_subsets = []
        for subset_names in newscheme_description:
            #some names look like this "Gene1_pos1-Gene2_pos2", so we need to extract the part names
            partition_names = []
            for x in subset_names:            
                [partition_names.append(y) for y in x.split('-')]
            print partition_names
            sub = subset.Subset(*tuple([cfg.partitions.parts_by_name[i] for i in partition_names]))            
            print sub
            created_subsets.append(sub)
    
        scheme = (Scheme(cfg, str(scheme_name), *tuple(created_subsets)))
        print scheme
        return scheme
        
        
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
        
    print scheme_description
    print subs
    print cfg.partitions.parts_by_number
    print cfg.partitions.parts_by_name
    
    # We now have what we need to create a subset. Each entry will have a
    # set of values which are the index for the partition
    created_subsets = []
    for sub_indexes in subs.values():
        sub = subset.Subset(*tuple([cfg.partitions[i] for i in sub_indexes]))
        print sub
        created_subsets.append(sub)

    new_scheme = Scheme(cfg, str(scheme_name), *tuple(created_subsets))
	
    print new_scheme
	
    return new_scheme

def model_to_scheme(model, scheme_name, cfg):
	"""Turn a model definition e.g. [0, 1, 2, 3, 4] into a scheme"""
	import subset
	
	subs = {}
	# We use the numbers returned to group the different subsets
	for sub_index, grouping in enumerate(model):
		insub = subs.setdefault(grouping, [])
		insub.append(sub_index)
	# We now have what we need to create a subset. Each entry will have a
	# set of values which are the index for the partition
	created_subsets = []
	for sub_indexes in subs.values():
		sub = subset.Subset(*tuple([cfg.partitions[i] for i in sub_indexes]))
		created_subsets.append(sub)

	scheme = (Scheme(cfg, str(scheme_name), *tuple(created_subsets)))
	#log.info("Created scheme %d of %d" %(scheme_name, len(all_schemes)))	
	return scheme

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

