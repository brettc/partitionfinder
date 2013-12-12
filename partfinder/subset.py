# Copyright (C) 2012-2013 Robert Lanfear and Brett Calcott
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

import logtools
log = logtools.get_logger(__file__)

import os
import weakref
import numpy

from alignment import Alignment, SubsetAlignment
from util import PartitionFinderError, remove_runID_files, get_aic, get_aicc, get_bic
import subset_ops
import database


FRESH, PREPARED, DONE = range(3)


class SubsetError(PartitionFinderError):
    pass


def count_subsets():
    return len(Subset._cache)


def clear_subsets():
    Subset._cache.clear()


class Subset(object):
    """Contains a set of columns in the Alignment
    """
    _cache = weakref.WeakValueDictionary()

    def __new__(cls, cfg, column_set):
        """Returns the identical subset if the columns are identical.

        This is basically a pythonized factory. See here:
        http://codesnipers.com/?q=python-flyweights
        """
        column_set = list(column_set)
        column_set.sort()
        column_set = frozenset(column_set)
        obj = Subset._cache.get(column_set, None)
        if not obj:
            obj = object.__new__(cls)
            Subset._cache[column_set] = obj
            obj.init(cfg, column_set)

        return obj

    def init(self, cfg, column_set):
        self.cfg = cfg
        self.column_set = column_set
        self.columns = list(column_set)
        self.columns.sort()
        self.status = FRESH

        # We put all results into this array, which is sized to the number of
        # models that we are analysing
        self.result_array = numpy.zeros(
            cfg.model_count, cfg.data_layout.data_type)

        # This points to the current empty array entry that we will fill up
        # next. When we're done it will equal the size of the array (and
        # accessing will cause an error)
        self.result_current = 0

        # This will get set to the best entry, once we've done all the analysis
        self.result_best = None

        self.models_not_done = set(cfg.models)

        self.fabricated = False
        self.analysis_error = None
        self.centroid = None

        # Site likelihoods calculated using GTR+G from the
        # processor.gen_per_site_stats()
        self.site_lnls_GTRG = []

        self.alignment_path = None
        log.debug("Created %s" % self)

    def add_description(self, name, description):
        """User created subsets can get some extra info"""
        self.full_name = name
        self.description = description

    @property
    def long_name(self):
        try:
            l = self.full_name[:]
            long_name = ', '.join(l)
            return long_name
        except:
            return "NA"

    @property
    def site_description(self):
        try:
            s = []
            for d in self.description:
                d = d[0]
                if d == 1:
                    text = "%s-%s" % (d[0], d[1])
                else:
                    text = "%s-%s\\%s" % tuple(d)
                s.append(text)
            site_description = ', '.join(s)
            return site_description 
        except:
            return ', '.join(map(str, self.columns))

    def __repr__(self):
        return "Subset(%s..)" % self.name[:5]

    @property
    def name(self):
        # Cache this
        if hasattr(self, '_name'):
            nm = self._name
        else:
            nm = subset_ops.subset_unique_name(self)
            self._name = nm
        return nm

    def load_results(self, cfg):
        matching = cfg.database.get_results_for_subset(self)
        # We might get models that we don't want, so we need to filter them
        for i, mod in enumerate(matching['model_id']):
            if mod in self.models_not_done:
                self.result_array[self.result_current] = matching[i]
                self.result_current += 1
                self.models_not_done.remove(mod)

    def add_result(self, cfg, model, result):
        """
        We get the result class from raxml or phyml. We need to transform this
        into a numpy record, and then store it locally, and in the database
        """
        K = float(cfg.processor.models.get_num_params(model))
        n = float(len(self.column_set))
        lnL = float(result.lnl)
        aic = get_aic(lnL, K)
        bic = get_bic(lnL, K, n)
        aicc = get_aicc(lnL, K, n)

        result.subset_id = self.name
        result.model_id = model
        result.params = K
        result.aic = aic
        result.aicc = aicc
        result.bic = bic

        # Now assign the data record in the result into the current model
        # result and save it to the database
        self.result_array[self.result_current] = result._data
        cfg.database.save_result(self, self.result_current)
        self.result_current += 1

        # We need to get the keys from the dictionary and convert them to
        # indexes into our fixed sized arrays
        # row_freqs = row['freqs']
        # getter = database.Freqs.get
        # for k, v in result.freqs.items():
        #     i = getter(k)
        #     if i is None:
        #         log.error("Unrecognised Frequency type %s", k)
        #         raise SubsetError
        #     row_freqs[i] = v

        # row_rates = row['rates']
        # getter = database.Rates.get
        # for k, v in result.rates.items():
        #     i = getter(k)
        #     if i is None:
        #         log.error("Unrecognised Rate type %s", k)
        #         raise SubsetError
        #     row_rates[i] = v


        log.debug("Added model to subset. Model: %s, params: %d, sites:%d, lnL:%.2f, site_rate %f"
                  % (model, K, n, lnL, result.site_rate))

    def model_selection(self, cfg):
        # We want the index of the smallest value
        method = cfg.model_selection
        self.result_best = numpy.argmin(self.result_array[method])
        best = self.result_array[self.result_best]

        # TODO: this is crappy. Anyone who wants this stuff should just access
        # the entire "best" item
        self.best_info_score = best[method]
        self.best_lnl = best['lnl']
        self.best_model = best['model_id']
        self.best_site_rate = best['site_rate']
        self.best_params = best['params']
        self.best_alpha = best['alpha']
        self.best_freqs = best['freqs']
        self.best_rates = best['rates']

        log.debug("Model Selection. best model: %s, params: %d, site_rate: %f"
                  % (self.best_model, self.best_params, self.best_site_rate))

    def get_param_values(self):
        param_values = {
            "rate": self.best_site_rate,
            "model": self.best_model,
            "alpha": self.best_alpha,
            "freqs": self.best_freqs,
            "rates": self.best_rates,
        }
        return param_values

    def finalise(self, cfg):
        if self.models_not_done:
            return False

        # We might already have done everything
        if self.status == DONE:
            return True

        self.model_selection(cfg)

        # Do all the final cleanup
        if cfg.save_phylofiles:
            # Write out a summary of the subsets
            cfg.reporter.write_subset_summary(self)
        else:
            # Otherwise, clean up files generated by the programs as well
            if self.alignment_path:
                remove_runID_files(self.alignment_path)

        self.models_to_process = []
        self.status = DONE
        cfg.progress.subset_done(self)
        return True

    def prepare(self, cfg, alignment):
        """Get everything ready for running the analysis
        """
        cfg.progress.subset_begin(self)

        # Load the cached results
        self.load_results(cfg)
        if self.finalise(cfg):
            return

        # Make an Alignment from the source, using this subset
        self.make_alignment(cfg, alignment)
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
            result = cfg.processor.parse(output, cfg)
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
                # screwed. Reraise the message
                log.error(
                    "Failed to run models %s; not sure why",
                    ", ".join(list(self.models_not_done)))
                raise

    def add_per_site_statistics(self, per_site_stats):
        self.site_lnls = per_site_stats[0]
        self.site_rates = per_site_stats[2]
        self.lnls_rates = per_site_stats[3]
        self.lnls_rate_cats = per_site_stats[1]

    def fabricate_result(self, cfg, model):
        '''If the subset fails to be analyzed, we throw some "fabricated"
        results'''
        processor = cfg.processor
        self.fabricated = True

        lnl = sum(self.site_lnls_GTRG)
        result = processor.fabricate(lnl)

        self.add_result(cfg, model, result)
        self.best_params = cfg.processor.models.get_num_params(model)
        self.best_lnl = result.lnl
        self.models_not_done.remove(model)

    def add_centroid(self, centroid):
        self.centroid = centroid

    FORCE_RESTART_MESSAGE = """
    It looks like you have changed one or more of the data_blocks in the
    configuration file, so the new subset alignments don't match the ones
    stored for this analysis.  You'll need to run the program with
    --force-restart
    """

    def make_alignment(self, cfg, alignment):
        # Make an Alignment from the source, using this subset
        sub_alignment = SubsetAlignment(alignment, self)
        sub_path = os.path.join(cfg.phylofiles_path, self.name + '.phy')
        # Add it into the sub, so we keep it around
        self.alignment_path = sub_path

        # Maybe it is there already?
        if os.path.exists(sub_path):
            log.debug("Found existing alignment file %s" % sub_path)
            old_align = Alignment()
            old_align.read(sub_path)

            # It had better be the same!
            if not old_align.same_as(sub_alignment):
                log.error(self.FORCE_RESTART_MESSAGE)
                raise SubsetError
        else:
            # We need to write it
            sub_alignment.write(sub_path)

    @property
    def is_done(self):
        return self.status == DONE

    @property
    def is_prepared(self):
        return self.status == PREPARED
