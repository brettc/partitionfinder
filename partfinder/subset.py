#Copyright (C) 2012 Robert Lanfear and Brett Calcott
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
#program, the RAxML program, and the PyParsing library,
#all of which are protected by their own licenses and conditions, using
#PartitionFinder implies that you agree with those licences and conditions as well.

import logging
log = logging.getLogger("subset")
import os
import weakref

from hashlib import md5

# import base64
# from zlib import compress

import cPickle as pickle
from math import log as logarithm
from alignment import Alignment, SubsetAlignment
from util import PartitionFinderError, remove_runID_files

FRESH, PREPARED, DONE = range(3)


class SubsetError(PartitionFinderError):
    pass


def count_subsets():
    return 1
    # return len(Subset._cache)
    #


def clear_subsets():
    pass
    # Subset._cache.clear()


class Subset(object):
    """A Subset of Partitions
    """
    # TODO: Move this to the config -- once we have a global one
    _cache = weakref.WeakValueDictionary()

    def __new__(cls, *parts):
        """Return the SAME subset if the partitions are identical. This is
        basically a pythonized factory. See here:
        http://codesnipers.com/?q=python-flyweights
        """

        cacheid = frozenset(parts)
        obj = Subset._cache.get(cacheid, None)
        # TODO Flush cache? USE MRU? functools.lrucache
        if not obj:
            obj = object.__new__(cls)
            Subset._cache[cacheid] = obj
            obj.init(cacheid, *parts)

        # obj = object.__new__(cls)
        # cacheid = frozenset(parts)
        # obj.init(cacheid, *parts)
        return obj

    def init(self, cacheid, *parts):
        # Error checking....
        self.status = FRESH

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

        self.partitions = cacheid

        # a list of columns in the subset
        self.columns = []
        self.columnset = set()
        for p in parts:
            self.columns += p.columns
            self.columnset |= p.columnset
        self.columns.sort()

        self.results = {}
        self.best_info_score = None  # e.g. AIC, BIC, AICc
        self.best_model = None
        self.best_params = None
        self.best_lnl = None
        self.alignment_path = None
        log.debug("Created %s", self)

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

    def add_result(self, cfg, model, result):
        result.model = model
        result.params = cfg.processor.models.get_num_params(model)

        K = float(result.params)
        n = float(len(self.columnset))
        lnL = float(result.lnl)
        #here we put in a catch for small subsets, where n<K+2
        #if this happens, the AICc actually starts rewarding very small datasets, which is wrong
        #a simple but crude catch for this is just to never allow n to go below k+2
        result.aic = (-2.0 * lnL) + (2.0 * K)
        result.bic = (-2.0 * lnL) + (K * logarithm(n))


        if n < (K + 2):
            if cfg.model_selection.lower() == "aicc":
                log.warning("The subset containing the following data_blocks: %s, has a very small"
                        " number of sites (%d) compared to the number of parameters"
                        " in the model being estimated (the %s model which has %d parameters)."
                        " This may give misleading AICc results, so please check carefully"
                        " if you are using the AICc for your analyses."
                        " The model selection results for this subset are in the following file:"
                        " /analysis/subsets/%s.txt\n" % (self, n, model, K, self.name))
            n = K + 2

        result.aicc = (-2.0 * lnL) + ((2.0 * K) * (n / (n - K - 1.0)))

        #this is the rate per site of the model - used in some clustering analyses
        result.site_rate = float(result.tree_size)

        log.debug("Adding model to subset. Model: %s, params %d, site_rate %f" % (model, K, result.site_rate))

        if model in self.results:
            log.error("Can't add model result %s, it already exists in %s",
                      model, self)
        self.results[model] = result

    def model_selection(self, cfg):
        # Model selection is done after we've added all the models
        # Note: we may have more models than we want if there is old data lying
        # around
        self.best_info_score = None  # Reset this before model selection
        meth = cfg.model_selection.lower()

        for model in cfg.models:
            result = self.results[model]
            try:
                info_score = getattr(result, meth)
            except AttributeError:
                log.error(
                    "Model selection option %s not recognised, "
                    "please check" % cfg.model_selection)
                raise SubsetError

            if self.best_info_score is None or info_score < self.best_info_score:
                self.best_lnl = result.lnl
                self.best_info_score = info_score
                self.best_model = result.model
                self.best_params = result.params
                self.best_site_rate = result.site_rate
                self.best_alpha = result.alpha
                self.best_freqs = result.freqs
                self.best_modelparams = result.rates

        log.debug("Model Selection. best model: %s, params: %d, site_rate: %f"
                  % (self.best_model, self.best_params, self.best_site_rate))

    def get_param_values(self):
        param_values = {}

        param_values["rate"] = self.best_site_rate
        param_values["alpha"] = self.best_alpha

        #not sure if this sorting is necessary, but it's here in case it's needed
        #to make sure that freqs and model parameters are always in the same order
        #can't hurt... (i hope).
        keys_f = self.best_freqs.keys()
        keys_f.sort()
        param_values["freqs"] = [self.best_freqs[key] for key in keys_f]

        keys_m = self.best_modelparams.keys()
        keys_m.sort()
        param_values["model"] = [self.best_modelparams[key] for key in keys_m]

        return param_values

    def finalise(self, cfg):
        if self.models_not_done:
            return False

        # We might already have done everything
        if self.status == DONE:
            return True

        # Do all the final cleanup
        cfg.reporter.write_subset_summary(self)
        self.save_results(cfg)
        self.model_selection(cfg)
        if not cfg.save_phylofiles:
            remove_runID_files(self.alignment_path)

        self.models_to_process = []
        self.status = DONE
        cfg.progress.subset_done(self)
        return True

    def prepare(self, cfg, alignment):
        """Get everything ready for running the analysis
        """
        # cfg.progress.update_subsets(self)
        cfg.progress.subset_begin(self)

        # Load the cached results
        self.load_results(cfg)

        # First, see if we've already got the results loaded. Then we can
        # shortcut all the other checks
        models_done = set(self.results.keys())
        self.models_not_done = cfg.models - models_done
        if self.finalise(cfg):
            return

        # Make an Alignment from the source, using this subset
        self.make_alignment(cfg, alignment)

        # Try and read in some previous analyses
        self.parse_results(cfg)
        if self.finalise(cfg):
            return

        self.models_to_process = list(self.models_not_done)
        # Now order them by difficulty
        self.models_to_process.sort(
            key=cfg.processor.models.get_model_difficulty,
            reverse=True)

        self.status = PREPARED

    def parse_results(self, cfg):
        """Read in the results and parse them"""
        for m in list(self.models_not_done):
            self.parse_model_result(cfg, m)

    def parse_model_result(self, cfg, model):
        pth, tree_path = cfg.processor.make_output_path(
            self.alignment_path, model)

        if not os.path.exists(pth):
            # If it ain't there, we can't do it
            return

        output = open(pth, 'rb').read()
        try:
            result = cfg.processor.parse(output, cfg.datatype)
            self.add_result(cfg, model, result)
            # Remove the current model from remaining ones
            self.models_not_done.remove(model)

            # Just used for below
            if not cfg.save_phylofiles:
                # We remove all files that have the specified RUN ID
                cfg.processor.remove_files(self.alignment_path, model)

        except cfg.processor.PhylogenyProgramError:
            # If we're loading old files, this is fine
            if self.status == FRESH:
                log.warning("Failed loading parse output from %s."
                            "Output maybe corrupted. I'll run it again.",
                            pth)
                cfg.processor.remove_files(self.alignment_path, model)
            else:
                # But if we're prepared, then we've just run this. And we're
                # screwed
                log.error(
                    "Failed to run models %s; not sure why",
                    ", ".join(list(self.models_not_done)))
                raise

    def make_alignment(self, cfg, alignment):
        # Make an Alignment from the source, using this subset
        sub_alignment = SubsetAlignment(alignment, self)
        sub_path = os.path.join(cfg.phylofiles_path, self.name + '.phy')
        # Add it into the sub, so we keep it around
        self.alignment_path = sub_path

        # Maybe it is there already?
        if os.path.exists(sub_path):
            log.debug("Found existing alignment file %s", sub_path)
            old_align = Alignment()
            old_align.read(sub_path)

            # It had better be the same!
            if not old_align.same_as(sub_alignment):
                log.error(
                    "It looks like you have changed one or more of the "
                    "data_blocks in the configuration file, "
                    "so the new subset alignments "
                    "don't match the ones stored for this analysis. "
                    "You'll need to run the program with --force-restart")
                raise SubsetError
        else:
            # We need to write it
            sub_alignment.write(sub_path)

    def get_subset_cache_path(self, cfg):
        return os.path.join(cfg.subsets_path, self.name + '.bin')

    def load_results(self, cfg):
        # We might have already saved a bunch of results, try there first
        if not self.results:
            log.debug("Reading in cached data from the subsets file")
            self.read_cache(self.get_subset_cache_path(cfg))

    def save_results(self, cfg):
        self.write_cache(self.get_subset_cache_path(cfg))

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
