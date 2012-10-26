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
log = logging.getLogger("analysis")

import os
import weakref

from hashlib import md5

# import base64
# from zlib import compress

import cPickle as pickle

import alignment
import phyml_models

from math import log as logarithm

from util import PartitionFinderError
class SubsetError(PartitionFinderError):
    pass

class Subset(object):
    """A Subset of Partitions
    """

    # TODO: changes this to AllSubsets?
    _cache = weakref.WeakValueDictionary()

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
            
            # a list of columns in the subset
            obj.columns = []
            obj.columnset = set()
            for p in parts:
                obj.columns += p.columns
                obj.columnset |= p.columnset
            obj.columns.sort()

            obj.results = {}
            obj.best_info_score = None #e.g. AIC, BIC, AICc
            obj.best_model = None
            obj.best_params = None
            obj.best_lnl = None
            obj.alignment_path = None
            log.debug("Created %s", obj)
        # else:
            # log.debug("Reused %s", obj)

        return obj

    # def __init__(self, *parts):
        # Everything is relegated to above...

    def __str__(self):
        return "(%s)" % ", ".join([str(p.name) for p in self.partitions])

    @property
    def full_name(self):
        if hasattr(self, '_full_name'):
            nm = self._full_name
        else:
            s = sorted([p.name for p in self.partitions])
            nm = '-'.join(s)
            self._full_name = nm
        return nm

    @property
    def name(self):
        # Cache this
        if hasattr(self, '_name'):
            nm = self._name
        else:
            nm = self.full_name
            # This gets super long -- we can shorten it like this...  This is
            # a slightly lazy solution. There is some vanishingly small chance
            # that we'll get the same thing. Google "MD5 Hash Collision"
            nm = md5(nm).hexdigest()
            self._name = nm
        return nm

    def __iter__(self):
        return iter(self.partitions)

    def add_model_result(self, model, result):
        result.model = model
        result.params = phyml_models.get_num_params(model)
    
        K = float(result.params)
        n = float(len(self.columnset))
        lnL = float(result.lnl)
        #here we put in a catch for small subsets, where n<K+2
        #if this happens, the AICc actually starts rewarding very small datasets, which is wrong
        #a simple but crude catch for this is just to never allow n to go below k+2
        if n<(K+2): 
            log.warning("The subset containing the following data_blocks: %s, has a very small"
                        " number of sites (%d) compared to the number of parameters"
                        " in the model being estimated (the %s model which has %d parameters)."
                        " This may give misleading AICc results, so please check carefully"
                        " if you are using the AICc for your analyses."
                        " The model selection results for this subset are in the following file:" 
                        " /analysis/subsets/%s.txt\n" % (self, n, model, K, self.name))
            n = K+2 

        result.aic  = (-2.0*lnL) + (2.0*K)
        result.bic  = (-2.0*lnL) + (K * logarithm(n))
        result.aicc = (-2.0*lnL) + ((2.0*K)*(n/(n-K-1.0)))

        #this is the rate per site of the model - used in some clustering analyses
        result.site_rate = float(result.tree_size)

        log.debug("Adding model to subset. Model: %s, params %d, site_rate %f" %(model, K, result.site_rate))

        if model in self.results:
            log.error("Can't add model result %s, it already exists in %s",
                    model, self)
        self.results[model] = result

    def model_selection(self, method, models):
        #model selection is done after we've added all the models
        self.best_info_score = None #reset this before model selection

        for model in models:
            result=self.results[model]
            if method=="AIC" or method=="aic":
                info_score = result.aic
            elif method=="AICc" or method=="AICC" or method=="aicc":
                info_score = result.aicc
            elif method=="BIC" or method=="bic":
                info_score = result.bic
            else:
                log.error("Model selection option %s not recognised, please check", method)
                raise SubsetError
		
            if self.best_info_score is None or info_score < self.best_info_score:
                self.best_lnl = result.lnl
                self.best_info_score = info_score
                self.best_model = result.model
                self.best_params = result.params
                self.best_site_rate = result.site_rate
        log.debug("Model Selection. best model: %s, params: %d, site_rate: %f" %(self.best_model, self.best_params, self.best_site_rate))

    def get_param_values(self):
        param_values =[]
        
        #add any parameters you want to this list
        param_values.append(self.best_site_rate)
        
        return param_values
        

    # These are the fields that get stored for quick loading
    _cache_fields = "alignment_path results".split()
    def write_cache(self, path):
        """Write out the results we've collected to a binary file"""
        f = open(path, 'wb')
        store = dict([(x, getattr(self, x)) for x in Subset._cache_fields])
        pickle.dump(store, f, -1)
        f.close()

    def read_cache(self, path):
        if not os.path.exists(path):
            return False

        log.debug("Reading binary cached results for %s", self)
        f = open(path, 'rb')
        self.__dict__.update(pickle.load(f))
        f.close()
