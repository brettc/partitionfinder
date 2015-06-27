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

import logtools
log = logtools.get_logger()

import subset_ops
import submodels

from util import PartitionFinderError, get_aic, get_aicc, get_bic


class SchemeError(PartitionFinderError):
    pass


class SchemeResult(object):
    def __init__(self, sch, nseq, branchlengths, model_selection):
        self.scheme_name = sch.name
        self.scheme = sch
        self.model_selection = model_selection

        # Calculate AIC, BIC, AICc for each scheme.
        # How you do this depends on whether brlens are linked or not.
        self.nsubs = len(sch.subsets)  # number of subsets
        sum_subset_k = sum([s.best_params for s in sch])  # sum of number of parameters in the best model of each subset

        log.debug("""Calculating number of parameters in scheme.
                  Total parameters from subset models: %d""" % (sum_subset_k))

        if branchlengths == 'linked':  # linked brlens - only one extra parameter per subset
            self.sum_k = sum_subset_k + (self.nsubs - 1) + (
                (2 * nseq) - 3)  # number of parameters in a scheme
            log.debug("Total parameters from brlens: %d" %
                      ((2 * nseq) - 3))
            log.debug("Parameters from subset multipliers: %d" %
                      (self.nsubs - 1))

        elif branchlengths == 'unlinked':  # unlinked brlens - every subset has its own set of brlens
            self.sum_k = sum_subset_k + (self.nsubs * (
                (2 * nseq) - 3))  # number of parameters in a scheme
            log.debug("Total parameters from brlens: %d" % ((
                2 * nseq) - 3) * self.nsubs)

        else:
            # WTF?
            log.error("Unknown option for branchlengths: %s", branchlengths)
            raise PartitionFinderError

        log.debug("Grand total parameters: %d" % (self.sum_k))

        self.lnl = sum([s.best_lnl for s in sch])
        self.nsites = sum([len(s.column_set) for s in sch])

        K = float(self.sum_k)
        n = float(self.nsites)
        lnL = float(self.lnl)

        log.debug("n: %d\tK: %d\tlnL: %d" % (n, K, lnL))

        self.aic = get_aic(lnL, K)
        self.bic = get_bic(lnL, K, n)
        self.aicc = get_aicc(lnL, K, n)

    @property
    def score(self):
        return getattr(self, self.model_selection)

    def __repr__(self):
        return "SchemeResult<score({0.model_selection}):{0.score}>".format(self)


class Scheme(object):
    def __init__(self, cfg, name, subsets, description=None):
        """A set of subsets of partitions"""
        self.name = name
        self.subsets = set(subsets)
        self.description = description

        # TODO: Fix this!
        if subset_ops.subsets_overlap(subsets):
           log.error("Scheme '%s' contains overlapping subsets", name)
           raise SchemeError
        #
        #if subset_ops.has_missing(subsets):
        #    log.error("Scheme '%s' has missing subsets", name)
        #    raise SchemeError
        #
        #log.debug("Created %s" % self)

    def __iter__(self):
        return iter(self.subsets)

    def __str__(self):
        ss = ', '.join([str(s) for s in self.subsets])
        return "Scheme(%s, %s)" % (self.name, ss)

    def get_fabricated_subsets(self):
        fabricated_subsets = []
        for sub in self.subsets:
            if sub.fabricated:
                fabricated_subsets.append(sub)
        return fabricated_subsets


class SchemeSet(object):
    """All the schemes added, and also a list of all unique subsets"""
    def __init__(self):
        """A collection of schemes"""
        self.clear_schemes()

    def clear_schemes(self):
        self.schemes_by_name = {}
        # self.schemes_by_subsets = {}

    def add_scheme(self, scheme):
        if scheme.name in self.schemes_by_name:
            log.error("Cannot add two schemes with same name: '%s'" %
                      scheme.name)
            raise SchemeError

        # TODO: Recheck schemes to make sure they're ok...
        # if scheme.part_subsets in self.schemes_by_subsets:
            # existing_scheme = \
                # self.schemes_by_subsets[scheme.part_subsets]
            # log.warning(
                # "Scheme named %s being added is identical to existing %s",
                # scheme.name, existing_scheme)
            # # raise SchemeError

        self.schemes_by_name[scheme.name] = scheme
        # self.schemes_by_subsets[scheme.part_subsets] = scheme

    def __len__(self):
        return len(self.schemes_by_name)

    def __iter__(self):
        return iter(self.schemes_by_name.itervalues())


def create_scheme(cfg, scheme_name, scheme_description):
    """
    Generate a single scheme given a list of numbers that represent the
    indexes of the partitions e.g. [0,1,2,3,4,5,6,7]
    """

    subset_count = len(cfg.user_subsets)

    # Check that the correct number of items are in the list
    if len(scheme_description) != subset_count:
        log.error("There's a problem with the description of scheme %s" %
                  scheme_name)
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
        subs_to_merge = [cfg.user_subsets[i] for i in sub_indexes]
        sub = subset_ops.merge_subsets(subs_to_merge)
        created_subsets.append(sub)

    return Scheme(cfg, str(scheme_name), created_subsets, description=scheme_description)


def model_to_scheme(model, scheme_name, cfg):
    """Turn a model definition e.g. [0, 1, 2, 3, 4] into a scheme"""
    subs = {}
    # We use the numbers returned to group the different subsets
    for sub_index, grouping in enumerate(model):
        insub = subs.setdefault(grouping, [])
        insub.append(sub_index)

    # We now have what we need to create a subset. Each entry will have a
    # set of values which are the index for the partition
    created_subsets = []
    for sub_indexes in subs.values():
        subs_to_merge = [cfg.user_subsets[i] for i in sub_indexes]
        sub = subset_ops.merge_subsets(subs_to_merge)
        created_subsets.append(sub)

    return Scheme(cfg, str(scheme_name), created_subsets)


def generate_all_schemes(cfg):
    """
    Convert the abstract schema given by the algorithm into subsets
    """

    log.info("Generating all possible schemes for the partitions...")

    subset_count = len(cfg.user_subsets)

    # Now generate the pattern for this many partitions
    all_schemes = submodels.get_submodels(subset_count)
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
            sub = subset_ops.merge_subsets(
                [cfg.user_subsets[i] for i in sub_indexes])
            created_subsets.append(sub)

        scheme_list.append(Scheme(cfg, str(scheme_name), created_subsets))

        log.debug("Created scheme %d of %d" % (scheme_name, len(all_schemes)))

        scheme_name += 1

    return scheme_list
